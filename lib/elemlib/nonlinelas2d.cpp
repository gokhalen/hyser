/*======================================================================
  Nachiket Gokhale gokhalen@bu.edu
  Two Dimensional Non-Linear Elasticity Element
  Builds stiffness matrix and right hand side
=======================================================================*/

#ifndef HAVE_NONLINELAS2D_HPP
#include "nonlinelas2d.hpp"
#endif

#include <cassert>
#include <cstdlib>

//=========================================================================
//implement constructor
nonlinelas2d::nonlinelas2d(int* cinte,int cnprop, int cmregparm, long int celemno){
  // have all information necessary. call constrctor of shapestrct

  // why is nprop set to a constant in here ?
  // localize nprop from inmesh!
  //nprop = 3;
  
  this->nprop = cnprop;
  this->mregparm = cmregparm;
  mdofn = 2;
  nonlinflag = true;

  shpstrct = new shapestruct(2,2,cinte,4,NONLINELAS2D,celemno);
  
  // allocate the rest after number of nodes has been set
  prop        = allocmat(shpstrct->numnodes,nprop,prop); 
  estiff      = allocmat(shpstrct->numnodes*mdofn,shpstrct->numnodes*mdofn,estiff);
  dofnatnode  = allocvec(shpstrct->numnodes,dofnatnode);

  // set the dofnatnode array, one degree of freedom at each node
  for(int inode = 0;inode<shpstrct->numnodes;inode++){
    dofnatnode[inode] = mdofn;
  }


  CIJ = allocmat(4,4,CIJ);
  krchstrs = allocmat(2,2,krchstrs);
  sptan    = alloc4D(2,2,2,2,sptan);
  // allocate primal
  primal  = allocmat(shpstrct->numnodes,mdofn,primal);
  zeromat(shpstrct->numnodes,mdofn,primal);

  lastupdate = allocmat(shpstrct->numnodes,mdofn,lastupdate);
  zeromat(shpstrct->numnodes,mdofn,lastupdate);

  // allocate dirichlet data and traction data vectors
  dirdata  = allocmat(shpstrct->numnodes,mdofn,dirdata);
  zeromat(shpstrct->numnodes,mdofn,dirdata);

  tracdata = allocmat(shpstrct->numnodes,mdofn,tracdata);
  zeromat(shpstrct->numnodes,mdofn,tracdata);

  // allocate right hand side vector
  erhs = allocvec(shpstrct->numnodes*mdofn,erhs);
  

  // allocate finite element arrays
  ideqn  = allocmat(shpstrct->numnodes,mdofn,ideqn); // (node num, dofn) -> glb eqn no
  idnode = allocvec(shpstrct->numnodes,idnode); // local -> global node num.
  isbc   = allocmat(shpstrct->numnodes,mdofn,isbc); //
  zeromat(shpstrct->numnodes,mdofn,isbc);


  // allocate space for inversion data
  measdata  = allocmat(shpstrct->numnodes,mdofn,measdata);
  zeromat(shpstrct->numnodes,mdofn,measdata);

  // allocate dual field data
  dualfield = allocmat(shpstrct->numnodes,mdofn,dualfield);
  zeromat(shpstrct->numnodes,mdofn,dualfield);

  // allocate the whole data - used in gradient calculation only
  // HACK for now - change with a constructor argument 
  primalwhole = alloc3D(shpstrct->numnodes,mdofn,6,primalwhole);
  zero3D(shpstrct->numnodes,mdofn,6,primalwhole);

  dualwhole = alloc3D(shpstrct->numnodes,mdofn,6,dualwhole);
  zero3D(shpstrct->numnodes,mdofn,6,dualwhole);

  // allocate gradient matrix
  
  gradient = allocmat(shpstrct->numnodes,nprop,gradient);
  zeromat(shpstrct->numnodes,nprop,gradient);

  // allocate the regparm array
  // can be done in the base class once the virtual function stuff
  // is figured out
  //regparm  = allocvec(this->nprop,regparm);

  functens = allocmat(shpstrct->ndime,shpstrct->ndime,functens);

  Tprim    = allocmat(shpstrct->numnodes,shpstrct->ndime,Tprim);
  Tmeas    = allocmat(shpstrct->numnodes,shpstrct->ndime,Tmeas);

  zeromat(shpstrct->numnodes,shpstrct->ndime,Tprim);
  zeromat(shpstrct->numnodes,shpstrct->ndime,Tmeas);

  // call memalloc for the common allocations 
  memalloc();
  
  }
//==========================================================================
nonlinelas2d::~nonlinelas2d(){

  deallocmat(4,CIJ);
  deallocmat(2,krchstrs);
  dealloc4D(2,2,2,2,sptan);

  deallocmat(shpstrct->numnodes,primal);
  deallocmat(shpstrct->numnodes,lastupdate);
  deallocmat(shpstrct->numnodes,prop);
  deallocmat(shpstrct->numnodes*mdofn,estiff);
  deallocvec(dofnatnode);
  deallocmat(shpstrct->numnodes,dirdata);
  deallocmat(shpstrct->numnodes,tracdata);
  deallocvec(erhs);

  deallocmat(shpstrct->numnodes,ideqn);
  deallocvec(idnode);
  deallocmat(shpstrct->numnodes,isbc);
  //shpstrct->~shapestruct();
  
   //deallocate measdata
  deallocmat(shpstrct->numnodes,measdata);

  // deallocate dualfield
  deallocmat(shpstrct->numnodes,dualfield);

  // deallocate the whole fields
  dealloc3D(shpstrct->numnodes,mdofn,primalwhole);
  dealloc3D(shpstrct->numnodes,mdofn,dualwhole);

  // deallocate gradient 
  deallocmat(shpstrct->numnodes,gradient);

  //deallocvec(regparm);
  deallocmat(shpstrct->ndime,functens);
  deallocmat(shpstrct->numnodes,Tprim);
  deallocmat(shpstrct->numnodes,Tmeas);

  memdealloc();
  
  delete shpstrct;

}
//========================================================================
// implement element stiffness routine
void nonlinelas2d::elemstiff(int which_stiff){
  
  // zero element stiffness K_{A,B} = \int_{\Omega}(N_{A,x}k(x)N_{B,X})
  zeromat(shpstrct->numnodes*mdofn,shpstrct->numnodes*mdofn,estiff);
  

  double *KK;
  // material displacement vector (at an integration point)
  double *disp;  
  double **defgrad;
  double **gradu;
  double **cauchy;
  
  double **firstpk;
  double **secondpiola;
  double **rightcauchygreen;
  double **rightcauchygreeninvtpose;
  double **linstrn;
  double **AIJ;
  double ****CIJKL;
  double ****AIJKL;
  double ****spatmod;

  
  // declaration of the conslib structure
  conselas23d *elasptr;

  // allocations, deallocated at the end of elemstiff()
  KK                  = allocvec(this->nprop,KK);
  disp                = allocvec(shpstrct->ndime,disp);
  defgrad             = allocmat(shpstrct->ndime,shpstrct->ndime,defgrad);
  gradu               = allocmat(shpstrct->ndime,shpstrct->ndime,gradu);
  linstrn             = allocmat(shpstrct->ndime,shpstrct->ndime,linstrn);
  AIJ                 = allocmat(shpstrct->ndime*shpstrct->ndime,shpstrct->ndime*shpstrct->ndime,AIJ);
  firstpk             = allocmat(shpstrct->ndime,shpstrct->ndime,firstpk);
  secondpiola         = allocmat(shpstrct->ndime,shpstrct->ndime,secondpiola);
  rightcauchygreen    = allocmat(shpstrct->ndime,shpstrct->ndime,rightcauchygreen);
  rightcauchygreeninvtpose = allocmat(shpstrct->ndime,shpstrct->ndime,rightcauchygreeninvtpose);
  CIJKL               = alloc4D(shpstrct->ndime,shpstrct->ndime,shpstrct->ndime,shpstrct->ndime,CIJKL);
   AIJKL               = alloc4D(shpstrct->ndime,shpstrct->ndime,shpstrct->ndime,shpstrct->ndime,AIJKL);

  // allocate memory for the elasptr
  elasptr             = (conselas23d *)(malloc(sizeof(conselas23d)));

  zeromat(shpstrct->ndime,shpstrct->ndime,secondpiola);

  // loop over integration points
  for(int iinte = 0;iinte<shpstrct->tinte;iinte++){

    // get jacobian determinant and weight
    
    double jacodet = shpstrct->jdet[iinte];
    double weight  = shpstrct->wts[iinte];

    // interpolate stiffness at this integration point
    
    // first zero property vector
    zerovec(this->nprop,KK);
    for(int iprop = 0;iprop<this->nprop;iprop++){
      for(int knode = 0;knode<shpstrct->numnodes;knode++){
	KK[iprop] +=  shpstrct->shape[iinte][knode]*prop[knode][iprop];      
      }
    }

    // interpolate the displacement at this point:
    // note: using dirdata is a redundant. primal has all the information
    zerovec(shpstrct->ndime,disp);

    for(int idime =0;idime<shpstrct->ndime;idime++){
      for(int knode =0;knode<shpstrct->numnodes;knode++){
	disp[idime] +=this->primal[knode][idime]*shpstrct->shape[iinte][knode];
      }
    }

    if(this->nprop==2){
      double lambda = KK[0];
      double mu     = KK[1];
    }

    // compute the gradient of the primal (displacement), GRADU, U(X),X 
    // we use the notation 
    // \gradu[i][j] = \partial{U_{i}}{X_{j}}
    // see TJRH ME 235 Stanford 
    
    zeromat(shpstrct->ndime,shpstrct->ndime,gradu);

    for(int idime=0;idime<shpstrct->ndime;idime++){
      for(int jdime=0;jdime<shpstrct->ndime;jdime++){
	
	for(int inode=0;inode<shpstrct->numnodes;inode++){
	  gradu[idime][jdime] += shpstrct->glbder[iinte][inode][jdime]*this->primal[inode][idime];
	  
	}// inode 
      } // jdime
    } // idime


    // compute the deformation gradient defgrad, \bf{1} + GRADU
    // calling the delta function is not the best (fastest) way to do this 
    // but by far the most clear 

    zeromat(shpstrct->ndime,shpstrct->ndime,defgrad);
    zeromat(shpstrct->ndime,shpstrct->ndime,linstrn);
    for(int idime = 0;idime<shpstrct->ndime;idime++){
      for(int jdime =0;jdime<shpstrct->ndime;jdime++){
	defgrad[idime][jdime] = (double)delta(idime,jdime) + gradu[idime][jdime];
	linstrn[idime][jdime] = (0.5)*(gradu[idime][jdime] + gradu[jdime][idime]);
      }
    }
    
    double dettdefgrad   = determinant(shpstrct->ndime, defgrad);
    if(dettdefgrad <= 0){
      cout <<__FILE__<<__LINE__<<" determinant of deformation gradient "<<dettdefgrad<<endl;
      abort();
    }
    
    // calculate the right cauchy green tensor F^{T}F
    
    tposetimesself(shpstrct->ndime,defgrad,rightcauchygreen);
    
    
      

    // we set elasptr data here               
    elasptr->grdflag  = 0;
    elasptr->ndime    =  shpstrct->ndime;
    elasptr->disp     = &disp[0];
    elasptr->defgrad  = defgrad;
    elasptr->rcgreen  = rightcauchygreen;
    elasptr->stfparam = &KK[0];
    elasptr->krchstrs = krchstrs;
    elasptr->secpk    = secondpiola;
    elasptr->cijklptr = CIJKL;
    elasptr->firstpk  = firstpk;
    elasptr->AIJ      = AIJ;
    elasptr->CIJ      = CIJ;
    elasptr->sptan    = sptan;
    elasptr->aijklptr = AIJKL;

    // we call elas2d -- this gives us the secondpk, C_IJKL, consistent tangent contrib
    elas2d(this->conseqnno,elasptr);
    
    // finally estiff assembly
    // refer equation 2.29 lengthy prospectus

    // first the CIJKL part
    
    for(int AA = 0;AA<shpstrct->numnodes;AA++){
      for(int ii = 0;ii<mdofn;ii++){
	int ieqn = mdofn*(AA-1+1) + ii;
	  
	  for(int BB=0;BB<shpstrct->numnodes;BB++){
	    for(int jj=0;jj<mdofn;jj++){
	      int jeqn = mdofn*(BB-1+1) + jj;
	      
	      for(int II=0;II<shpstrct->ndime;II++){
		double der1 = shpstrct->glbder[iinte][AA][II];

		for(int JJ=0;JJ<shpstrct->ndime;JJ++){
		  for(int KK=0;KK<shpstrct->ndime;KK++){
		    double der2 = shpstrct->glbder[iinte][BB][KK];

		    for(int LL=0;LL<shpstrct->ndime;LL++){
		      //cout<<__FILE__<<__LINE__<<" gets elem stiff core "<<endl;;
		      estiff[ieqn][jeqn] += der1*defgrad[ii][JJ]*CIJKL[JJ][II][KK][LL]*defgrad[jj][LL]*der2*jacodet*weight;
		    } //LL
		  }//KK
		}//JJ
	      }//BB
	    }//jj
	  }//BB
      }//ii
    }//AA
	    
    // second PK part
    // equation 2.42 expanded prospectus
    
    for(int AA = 0;AA<shpstrct->numnodes;AA++){
      for(int II=0;II<shpstrct->ndime;II++){

	double der1 = shpstrct->glbder[iinte][AA][II];

	for(int BB = 0;BB<shpstrct->numnodes;BB++){
	  for(int JJ = 0;JJ<shpstrct->ndime;JJ++){

	    double der2 = shpstrct->glbder[iinte][BB][JJ];

	    for(int ii = 0;ii<shpstrct->ndime;ii++){
	      int ieqn = mdofn*(AA-1+1) + ii;

	      for(int kk =0;kk<shpstrct->ndime;kk++){
		int jeqn = mdofn*(BB-1+1)+ kk;
		
		estiff[ieqn][jeqn] += der1*(secondpiola[JJ][II]*delta(ii,kk))*der2*jacodet*weight;
		
	      }
	    }

	  }
	}
      }
    }
    
  

    
    
  }// iinte

  deallocvec(KK);
  deallocvec(disp);
  deallocmat(shpstrct->ndime,gradu);
  deallocmat(shpstrct->ndime,defgrad);
  deallocmat(shpstrct->ndime,rightcauchygreen);
  deallocmat(shpstrct->ndime,firstpk);
  deallocmat(shpstrct->ndime,secondpiola);
  deallocmat(shpstrct->ndime,rightcauchygreeninvtpose);
  deallocmat(shpstrct->ndime,linstrn);
  deallocmat(shpstrct->ndime*shpstrct->ndime,AIJ);
  dealloc4D(shpstrct->ndime,shpstrct->ndime,shpstrct->ndime,shpstrct->ndime,CIJKL);
  dealloc4D(shpstrct->ndime,shpstrct->ndime,shpstrct->ndime,shpstrct->ndime,AIJKL);
  // deallocate memory for the elasptr
  free(elasptr);
}

//=======================================================================
//void nonlinelas2d::calcfunc(){
//  
///
//}

//=======================================================================
void nonlinelas2d:: elemrhs(){
  // this is very similar to building of stiffness matrix
  zerovec(shpstrct->numnodes*mdofn,erhs);


  double *KK;
  // material displacement vector (at an integration point)
  double *disp;  
  double **defgrad;
  double **gradu;
  double **firstpk;
  double **secondpiola;
  double **rightcauchygreen;
  double **rightcauchygreeninvtpose;
  double **linstrn;
  double ****CIJKL;
  double **AIJ;
  double ****AIJKL;
  
  // declaration of the conslib structure
  conselas23d *elasptr;

  // allocations, deallocated at the end of elemstiff()
  KK                  = allocvec(this->nprop,KK);
  disp                = allocvec(shpstrct->ndime,disp);
  defgrad             = allocmat(shpstrct->ndime,shpstrct->ndime,defgrad);
  gradu               = allocmat(shpstrct->ndime,shpstrct->ndime,gradu);
  linstrn             = allocmat(shpstrct->ndime,shpstrct->ndime,linstrn);
  AIJ                 = allocmat(shpstrct->ndime*shpstrct->ndime,shpstrct->ndime*shpstrct->ndime,AIJ);
  firstpk             = allocmat(shpstrct->ndime,shpstrct->ndime,firstpk);
  secondpiola         = allocmat(shpstrct->ndime,shpstrct->ndime,secondpiola);
  rightcauchygreen    = allocmat(shpstrct->ndime,shpstrct->ndime,rightcauchygreen);
  rightcauchygreeninvtpose = allocmat(shpstrct->ndime,shpstrct->ndime,rightcauchygreeninvtpose);
  CIJKL               = alloc4D(shpstrct->ndime,shpstrct->ndime,shpstrct->ndime,shpstrct->ndime,CIJKL);
  AIJKL               = alloc4D(shpstrct->ndime,shpstrct->ndime,shpstrct->ndime,shpstrct->ndime,AIJKL);
  // allocate memory for the elasptr
  elasptr             = (conselas23d *)(malloc(sizeof(conselas23d)));

  zeromat(shpstrct->ndime,shpstrct->ndime,secondpiola);

  // loop over integration points
  for(int iinte = 0;iinte<shpstrct->tinte;iinte++){

    // get jacobian determinant and weight
    
    double jacodet = shpstrct->jdet[iinte];
    double weight  = shpstrct->wts[iinte];

    // interpolate stiffness at this integration point
    
    // first zero property vector
    zerovec(this->nprop,KK);
    for(int iprop = 0;iprop<this->nprop;iprop++){
      for(int knode = 0;knode<shpstrct->numnodes;knode++){
	KK[iprop] +=  shpstrct->shape[iinte][knode]*prop[knode][iprop];      
      }
    }

    // interpolate the displacement at this point:
    // note: using dirdata is a redundant. primal has all the information
    zerovec(shpstrct->ndime,disp);

    for(int idime =0;idime<shpstrct->ndime;idime++){
      for(int knode =0;knode<shpstrct->numnodes;knode++){
	disp[idime] +=this->primal[knode][idime]*shpstrct->shape[iinte][knode];
      }
    }
    
    if(this->nprop==2){
      double lambda = KK[0];
      double mu     = KK[1];
    }

    // compute the gradient of the primal (displacement), GRADU, U(X),X 
    // we use the notation 
    // \gradu[i][j] = \partial{U_{i}}{X_{j}}
    // see TJRH ME 235 Stanford 
    
    zeromat(shpstrct->ndime,shpstrct->ndime,gradu);

    for(int idime=0;idime<shpstrct->ndime;idime++){
      for(int jdime=0;jdime<shpstrct->ndime;jdime++){
	
	for(int inode=0;inode<shpstrct->numnodes;inode++){
	  gradu[idime][jdime] += shpstrct->glbder[iinte][inode][jdime]*this->primal[inode][idime];
	  
	}// inode 
      } // jdime
    } // idime


    // compute the deformation gradient defgrad, \bf{1} + GRADU
    // calling the delta function is not the best (fastest) way to do this 
    // but by far the most clear 

    zeromat(shpstrct->ndime,shpstrct->ndime,defgrad);
    zeromat(shpstrct->ndime,shpstrct->ndime,linstrn);
    for(int idime = 0;idime<shpstrct->ndime;idime++){
      for(int jdime =0;jdime<shpstrct->ndime;jdime++){
	defgrad[idime][jdime] = (double)delta(idime,jdime) + gradu[idime][jdime];
	linstrn[idime][jdime] = (0.5)*(gradu[idime][jdime] + gradu[jdime][idime]);
      }
    }
    
    double dettdefgrad   = determinant(shpstrct->ndime, defgrad);
    if(dettdefgrad <= 0){
      cout <<__FILE__<<__LINE__<<" determinant of deformation gradient "<<dettdefgrad<<endl;
      abort();
    }
    
    // calculate the right cauchy green tensor F^{T}F
    
    tposetimesself(shpstrct->ndime,defgrad,rightcauchygreen);
    
    
      

    // we set elasptr data here
    elasptr->grdflag  = 0;
    elasptr->ndime    =  shpstrct->ndime;
    elasptr->disp     = &disp[0];
    elasptr->defgrad  = defgrad;
    elasptr->rcgreen  = rightcauchygreen;
    elasptr->stfparam = &KK[0];
    elasptr->krchstrs = krchstrs;
    elasptr->secpk    = secondpiola;
    elasptr->cijklptr = CIJKL;
    elasptr->firstpk  = firstpk;
    elasptr->AIJ      = AIJ;
    elasptr->CIJ      = CIJ;
    elasptr->secpk    = secondpiola;
    elasptr->sptan    = sptan;
    elasptr->cijklptr = CIJKL;
    elasptr->aijklptr = AIJKL;



    // we call elas2d -- this gives us the secondpk, C_IJKL, consistent tangent contrib
    elas2d(this->conseqnno,elasptr);
    
    // finally estiff assembly
    for(int AA = 0;AA<shpstrct->numnodes;AA++){
      for(int ii =0;ii<shpstrct->ndime;ii++){
	int ieqn  = mdofn*(AA-1+1)+ii;

	for(int II = 0;II<shpstrct->ndime;II++){
	  double der1 = shpstrct->glbder[iinte][AA][II];
	
	  for(int JJ = 0;JJ<shpstrct->ndime;JJ++){
	    erhs[ieqn] = erhs[ieqn]-der1*secondpiola[JJ][II]*defgrad[ii][JJ]*jacodet*weight;
	  }
	  
	}
      }
    }
  }

  //for(int ii = 0; ii<8;ii++){
  //  cout<<__FILE__<<" "<<__LINE__<<" "<<erhs[ii]<<endl;
  //}

  deallocvec(KK);
  deallocvec(disp);
  deallocmat(shpstrct->ndime,gradu);
  deallocmat(shpstrct->ndime,defgrad);
  deallocmat(shpstrct->ndime,rightcauchygreen);
  deallocmat(shpstrct->ndime,firstpk);
  deallocmat(shpstrct->ndime,secondpiola);
  deallocmat(shpstrct->ndime,rightcauchygreeninvtpose);
  deallocmat(shpstrct->ndime,linstrn);
  deallocmat(shpstrct->ndime*shpstrct->ndime,AIJ);
  dealloc4D(shpstrct->ndime,shpstrct->ndime,shpstrct->ndime,shpstrct->ndime,CIJKL);
  dealloc4D(shpstrct->ndime,shpstrct->ndime,shpstrct->ndime,shpstrct->ndime,AIJKL);

  // deallocate memory for the elasptr
  free(elasptr);
}

//=======================================================================
void nonlinelas2d:: calcgradient(){
  // we make use of the constitutive equation library.
  // again heavily based on the nonlinelas1d::calcgradient()

  // this is the data which will go into conslib
  // declarations specific to nonlinear problems;
  double *KK;
  double *propgrad; //store gradient of the i^th property
  // material displacement vector (at an integration point)
  double *disp;  
  double **defgrad;
  double **gradu;
  double **firstpk;
  double **secondpiola;
  double **rightcauchygreen;
  double **rightcauchygreeninvtpose;
  double **linstrn;
  double **AIJ;
  double ****CIJKL;
  double ****AIJKL;

  double **dualgrad; // material gradient of the dual field
  double **pk2grad;  // gradient of the second piola kirchoff wrt material properties

  
  // declare and allocate a pointer to a conselas1d structure
  conselas23d *elasptr; elasptr = (conselas23d *)malloc(sizeof(conselas23d));

  KK                  = allocvec(this->nprop,KK);
  propgrad            = allocvec(shpstrct->ndime,propgrad);
  disp                = allocvec(shpstrct->ndime,disp);
  defgrad             = allocmat(shpstrct->ndime,shpstrct->ndime,defgrad);
  gradu               = allocmat(shpstrct->ndime,shpstrct->ndime,gradu);
  linstrn             = allocmat(shpstrct->ndime,shpstrct->ndime,linstrn);
  firstpk             = allocmat(shpstrct->ndime,shpstrct->ndime,firstpk);
  secondpiola         = allocmat(shpstrct->ndime,shpstrct->ndime,secondpiola);
  rightcauchygreen    = allocmat(shpstrct->ndime,shpstrct->ndime,rightcauchygreen);
  AIJ                 = allocmat(shpstrct->ndime*shpstrct->ndime,shpstrct->ndime*shpstrct->ndime,AIJ);
  CIJKL               = alloc4D(shpstrct->ndime,shpstrct->ndime,shpstrct->ndime,shpstrct->ndime,CIJKL);
  AIJKL               = alloc4D(shpstrct->ndime,shpstrct->ndime,shpstrct->ndime,shpstrct->ndime,AIJKL);
  dualgrad            = allocmat(shpstrct->ndime,shpstrct->ndime,dualgrad);
  pk2grad             = allocmat(shpstrct->ndime,shpstrct->ndime,pk2grad);

  rightcauchygreeninvtpose = allocmat(shpstrct->ndime,shpstrct->ndime,rightcauchygreeninvtpose);  
 

  // loop over the integration points
  for(int iinte = 0; iinte<shpstrct->tinte; iinte++){
    double jacodet = shpstrct->jdet[iinte];
    double weight  = shpstrct->wts[iinte];
    
    // interpolate stiffness at this integration point
    // first zero property vector
    zerovec(this->nprop,KK);
    for(int iprop = 0;iprop<this->nprop;iprop++){
      for(int knode = 0;knode<shpstrct->numnodes;knode++){
	KK[iprop] +=  shpstrct->shape[iinte][knode]*prop[knode][iprop];      
      }
    }

    // interpolate the displacement at this point:
    // note: using dirdata is a redundant. primal has all the information
    zerovec(shpstrct->ndime,disp);

    for(int idime =0;idime<shpstrct->ndime;idime++){
      for(int knode =0;knode<shpstrct->numnodes;knode++){
	disp[idime] +=this->primal[knode][idime]*shpstrct->shape[iinte][knode];
      }
    }

    // compute the gradient of the primal (displacement), GRADU, U(X),X 
    // also compute the gradient of the dualfield W(X).
    // we use the notation 
    // \gradu[i][j] = \partial{U_{i}}{X_{j}}
    // see TJRH ME 235 Stanford, the unauthorized notes.
    
    zeromat(shpstrct->ndime,shpstrct->ndime,gradu);
    zeromat(shpstrct->ndime,shpstrct->ndime,dualgrad);
    for(int idime=0;idime<shpstrct->ndime;idime++){
      for(int jdime=0;jdime<shpstrct->ndime;jdime++){
	
	for(int inode=0;inode<shpstrct->numnodes;inode++){
	  gradu[idime][jdime] += shpstrct->glbder[iinte][inode][jdime]*this->primal[inode][idime];
	  dualgrad[idime][jdime] +=shpstrct->glbder[iinte][inode][jdime]*this->dualfield[inode][idime];
	}// inode 
      } // jdime
    } // idime

    
    // compute the deformation gradient defgrad, \bf{1} + GRADU
    // calling the delta function is not the best (fastest) way to do this 
    // but by far the most clear 
    
    zeromat(shpstrct->ndime,shpstrct->ndime,defgrad);
    zeromat(shpstrct->ndime,shpstrct->ndime,linstrn);
    for(int idime = 0;idime<shpstrct->ndime;idime++){
      for(int jdime =0;jdime<shpstrct->ndime;jdime++){
	defgrad[idime][jdime] = (double)delta(idime,jdime) + gradu[idime][jdime];
	linstrn[idime][jdime] = (0.5)*(gradu[idime][jdime] + gradu[jdime][idime]);
      }
    }
    
    // test if determinant of deformation gradient is positive
    double dettdefgrad   = determinant(shpstrct->ndime, defgrad);
    if(dettdefgrad <= 0){
      cout <<__FILE__<<__LINE__<<" determinant of deformation gradient "<<dettdefgrad<<endl;
      
      abort();
    }
    
    // calculate the right cauchy green tensor F^{T}F
    
    tposetimesself(shpstrct->ndime,defgrad,rightcauchygreen);

    // have all data required to initialize elasptr structure
    // we set elasptr data here
    elasptr->grdflag  = 1;
    elasptr->ndime    =  shpstrct->ndime;
    elasptr->disp     = &disp[0];
    elasptr->defgrad  = defgrad;
    elasptr->rcgreen  = rightcauchygreen;
    elasptr->stfparam = &KK[0];
    elasptr->krchstrs = krchstrs;
    elasptr->firstpk  = firstpk;
    elasptr->secpk    = secondpiola;
    elasptr->cijklptr = CIJKL;
    
    elasptr->AIJ      = AIJ;
    elasptr->CIJ      = CIJ;
    elasptr->sptan    = sptan;
    elasptr->aijklptr = AIJKL;
    elasptr->pk2grad  = pk2grad;
    

    // in the next loop we set which "stiffness parameter" we are taking the gradient wrt.

    for(int iprop = 0;iprop <nprop;iprop++){

      elasptr->grdparm = iprop;
      // fills pk2grad according to constitutive law number. see conslib:elas2d.cpp
      elas2d(this->conseqnno,elasptr); 

      for(int inode = 0;inode<shpstrct->numnodes;inode++){
	double gradshp = shpstrct->shape[iinte][inode];
	// refer expression from prospectus: equation 3.24 in expanded prospectus
	
	for(int II = 0; II<shpstrct->ndime; II++){
	  for(int JJ = 0; JJ<shpstrct->ndime; JJ++){

	    // si --> small i
	    for(int si = 0; si <shpstrct->ndime; si++){
	      this->gradient[inode][iprop] += dualgrad[si][II]*defgrad[si][JJ]*elasptr->pk2grad[JJ][II]*gradshp*jacodet*weight;
	      
	    }//si
	  }//JJ
	} // II

	// add the regularization contribution
	// H1-seminorm we need the gradient of \mu
	// The property to regularize is picked by iprop
	// the shape function is picked by inode

	// first calculate the gradient of this (iprop^th) property
	zerovec(shpstrct->ndime,propgrad);

	for (int jdime = 0; jdime < shpstrct->ndime; jdime++){
	  for (int jnode =0; jnode <shpstrct->numnodes; jnode++){
	    propgrad[jdime] += shpstrct->glbder[iinte][jnode][jdime]*this->prop[jnode][iprop];
	  }
	}
	// calculate dot product of gradient of properties
	double dotprdprop = pow(vecnorm(shpstrct->ndime,propgrad),2.0);
	

	
	for (int idime = 0; idime < shpstrct->ndime; idime++){
	  // calculate \mu_{idime}
	  double mugrad = 0.0;
	  for (int jnode = 0;jnode <shpstrct->numnodes;jnode++){
	    mugrad += shpstrct->glbder[iinte][jnode][idime]*this->prop[jnode][iprop];
	  }
	  for (int jdime = 0; jdime < shpstrct->ndime; jdime++){
	    if(jdime == idime){
	      //double first = this->regparm[iprop][0];
	      //double last  = jacodet*weight;
	      
	      switch(this->which_regularization){
	      case H1SEMINORM:
		//this->gradient[inode][iprop] += first*mugrad*shpstrct->glbder[iinte][inode][jdime]*last;
		break;

	      case TOTALVARIATION:
		//double oosqrt = (1.0/sqrt(dotprdprop + 0.01));
		//this->gradient[inode][iprop] += 0.5*first*oosqrt*propgrad[jdime]*shpstrct->glbder[iinte][inode][jdime]*jacodet*weight;
		break;
		
	      }//switch

	    }
	  }
	}
	

      }// inode
    } // iprop
  }//iinte
  
  deallocvec(KK);
  deallocvec(propgrad);
  deallocvec(disp);
  deallocmat(shpstrct->ndime,gradu);
  deallocmat(shpstrct->ndime,defgrad);
  deallocmat(shpstrct->ndime,rightcauchygreen);
  deallocmat(shpstrct->ndime,firstpk);
  deallocmat(shpstrct->ndime,secondpiola);
  deallocmat(shpstrct->ndime,rightcauchygreeninvtpose);
  deallocmat(shpstrct->ndime,linstrn);
  deallocmat(shpstrct->ndime*shpstrct->ndime,AIJ);
  dealloc4D(shpstrct->ndime,shpstrct->ndime,shpstrct->ndime,shpstrct->ndime,CIJKL);
  dealloc4D(shpstrct->ndime,shpstrct->ndime,shpstrct->ndime,shpstrct->ndime,AIJKL);
  deallocmat(shpstrct->ndime,dualgrad);
  deallocmat(shpstrct->ndime,pk2grad);

  // free the elasptr structure
  free(elasptr);
}

//======================================================================
void nonlinelas2d::getinteg(){
  // get integration points
  gauss2d(shpstrct);
}
//====================================================================
void nonlinelas2d::getshape(){
  //cout<<"getshape in nonlinelas2d()"<<endl;
  shape2d(shpstrct);
}
//======================================================================
void nonlinelas2d::printeltype(){
  cout <<endl<<"This is NONLINELAS2D"<<endl;
}
//=====================================================================
//void nonlinelas2d::memalloc(){
  // material properties, connectivity, stiffnessmatrix, righthand side
  // ideqn,idnode,isbc
  // The best way is to inherit this..
  
//}
//=====================================================================
//void nonlinelas2d::memdealloc(){
  // dealloc memory

//}

//==============================================================
//void nonlinelas2d::calcdual(){
//  
//  
//}
//==============================================================
//gvoid nonlinelas2d::setdualfield(double ddata,int inode,int idofn){}
//===============================================================
//void nonlinelas2d::calculate_Tmeas_Tprim(){

//}
//========================================================
// Old code for the dual forcing

// pick dimension to force
    //for(int jdime = 0;jdime<shpstrct->ndime;jdime++){
    //for(int jnode =0;jnode <shpstrct->numnodes;jnode++){

    //	for(int idime = 0;idime<shpstrct->ndime;idime++){
    //	  double prim = 0;
    //	  double meas = 0;
    //	  // interpolate primal and dual at this point
    //	  for(int inode = 0;inode<shpstrct->numnodes;inode++){
    //	    prim += Tprim[inode][idime]*shpstrct->shape[iinte][inode];
    //	    meas += Tmeas[inode][idime]*shpstrct->shape[iinte][inode];
    //	  } // inode
    //	  // interpolation complete, now force the dual
    //	  if(jdime == idime){
    //	    int ieqn = mdofn*(jnode - 1 + 1) + idime; 
    //	       double Tcontrib   = functens[][];
    // 	    this->erhs[ieqn] += -(prim-meas)*shpstrct->shape[iinte][jnode]*jacodet*weight;
    //	  }
    // 	}//idime
	
    //} // jnode
    //} // jdime

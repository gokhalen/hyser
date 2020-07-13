 /*======================================================================
  Nachiket Gokhale gokhalen@bu.edu
  Three Dimensional Non-Linear Elasticity Element
  Builds stiffness matrix and right hand side
=======================================================================*/

#ifndef HAVE_NONLINELAS3D_HPP
#include "nonlinelas3d.hpp"
#endif

#include <cassert>
#include <cstdlib>

//=========================================================================
//implement constructor
nonlinelas3d::nonlinelas3d(int* cinte, int cnprop, int cmregparm, long int celemno){
  // have all information necessary. call constrctor of shapestrct
  // we put the number of properties in here, it is variable for non-linear
  nprop = cnprop;
  this->mregparm = cmregparm;
  mdofn = 3;
  nonlinflag = true;

  int tempdime = 3;
  int tempnumnodes = 8;
  
  shpstrct = new shapestruct(tempdime,tempdime,cinte,tempnumnodes,NONLINELAS3D,celemno);

  // allocate the rest after number of nodes has been set
  prop        = allocmat(shpstrct->numnodes,nprop,prop); 
  estiff      = allocmat(shpstrct->numnodes*mdofn,shpstrct->numnodes*mdofn,estiff);
  dofnatnode  = allocvec(shpstrct->numnodes,dofnatnode);

  // set the dofnatnode array, one degree of freedom at each node
  for(int inode = 0;inode<shpstrct->numnodes;inode++){
    dofnatnode[inode] = mdofn;
  }

  // allocate primal
  primal  = allocmat(shpstrct->numnodes,mdofn,primal);
  zeromat(shpstrct->numnodes,mdofn,primal);

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

  // allocate gradient matrix
  
  gradient = allocmat(shpstrct->numnodes,nprop,gradient);
  zeromat(shpstrct->numnodes,nprop,gradient);

  }
//==========================================================================
nonlinelas3d::~nonlinelas3d(){

  deallocmat(shpstrct->numnodes,prop);
  deallocmat(shpstrct->numnodes*mdofn,estiff);
  deallocvec(dofnatnode);
  deallocmat(shpstrct->numnodes,primal);
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

  // deallocate gradient 
  deallocmat(shpstrct->numnodes,gradient);

  delete shpstrct;
}
//========================================================================
// implement element stiffness routine
void nonlinelas3d::elemstiff(int which_stiffness){
  
  // declare \& allocate memory for elasptr
  conselas23d *elasptr;
  elasptr = (conselas23d *)malloc(sizeof(conselas23d));

  // zero element stiffness K_{A,B} = \int_{\Omega}(N_{A,x}k(x)N_{B,X})
  zeromat(shpstrct->numnodes*mdofn,shpstrct->numnodes*mdofn,estiff);
  zerovec(shpstrct->numnodes*mdofn,erhs);

  // declarations specific to nonlinear problems;
  double *KK;
  double *disp;
  double **defgrad;
  double **gradu;
  double **secondpiola;
  double **rightcauchygreen;
  double **rightcauchygreeninvtpose;
  double **linstrn;
  double ****CIJKL;

  // allocations, deallocated at the end of elemstiff()
  KK                       = allocvec(this->nprop,KK);
  disp                     = allocvec(shpstrct->ndime,disp);
  defgrad                  = allocmat(shpstrct->ndime,shpstrct->ndime,defgrad);
  gradu                    = allocmat(shpstrct->ndime,shpstrct->ndime,gradu);
  linstrn                  = allocmat(shpstrct->ndime,shpstrct->ndime,linstrn);
  secondpiola              = allocmat(shpstrct->ndime,shpstrct->ndime,secondpiola);
  rightcauchygreen         = allocmat(shpstrct->ndime,shpstrct->ndime,rightcauchygreen);
  rightcauchygreeninvtpose = allocmat(shpstrct->ndime,shpstrct->ndime,rightcauchygreeninvtpose);
  CIJKL                    = alloc4D(shpstrct->ndime,shpstrct->ndime,shpstrct->ndime,shpstrct->ndime,CIJKL);
  
  zeromat(shpstrct->ndime,shpstrct->ndime,secondpiola);

  // loop over integration points
  for(int iinte = 0;iinte<shpstrct->tinte;iinte++){

    // get jacobian determinant and weight
    
    double jacodet = shpstrct->jdet[iinte];
    double weight  = shpstrct->wts[iinte];
    double wtjac  = jacodet*weight;

    // interpolate stiffness at this integration point
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

    double lambda = KK[0];
    double mu     = KK[1];

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
    
    // have rightcauchygreen now calculate the second PK  and CIJKL
    // first the second PK
    /*{
      // invert 
      invertsmall(shpstrct->ndime, rightcauchygreen, rightcauchygreeninvtpose);

      // take transpose
      tposeinplace(shpstrct->ndime,rightcauchygreeninvtpose);

      for(int irow = 0; irow<shpstrct->ndime;irow++){
	for(int icol = 0;icol<shpstrct->ndime;icol++){
	  double dett   = determinant(shpstrct->ndime, rightcauchygreen);
	  double right  = rightcauchygreeninvtpose[irow][icol];
	  double secpiola  = lambda*(dett - 1)*0.5*right + mu*(delta(irow,icol) -right);
	  
	  // secondpiola[irow][icol] = secpiola;
	  
	}
      }
      
    } // close PK scope


    for(int ii  =0;ii<shpstrct->ndime;ii++){
      for(int jj =0;jj<shpstrct->ndime;jj++){
	assert(secondpiola[ii][jj]==secondpiola[jj][ii]);
	assert(rightcauchygreeninvtpose[ii][jj] == rightcauchygreeninvtpose[jj][ii]);
	assert(rightcauchygreen[ii][jj] == rightcauchygreen[jj][ii]);
      }
    }

    // evaluate CIJKL
    // Oden material 
    {
      
      for(int ii = 0;ii<shpstrct->ndime;ii++){
	for(int jj = 0;jj<shpstrct->ndime;jj++){
	  for(int kk = 0;kk<shpstrct->ndime;kk++){
	    for(int ll = 0;ll<shpstrct->ndime;ll++){

	      double dett   = determinant(shpstrct->ndime, rightcauchygreen);

	      double symmfourth1  = rightcauchygreeninvtpose[ii][kk]*rightcauchygreeninvtpose[jj][ll];
	      double symmfourth2  = rightcauchygreeninvtpose[jj][kk]*rightcauchygreeninvtpose[ii][ll];
	      double symmfourth   = 0.5*(symmfourth1+symmfourth2);

	      double Cinv_ijCinv_kl  = rightcauchygreeninvtpose[ii][jj]*rightcauchygreeninvtpose[kk][ll];
	      
	      
	      double first  = (mu-0.5*lambda*(dett-1))*symmfourth;
	      double second = 0.5*lambda*dett*Cinv_ijCinv_kl;

	      double lintrace = tracemat(shpstrct->ndime, linstrn);
	      
		

	      //CIJKL[ii][jj][kk][ll] = 2*(first) + 2*(second);

	     
	     
	    }
	  }
	}
      }
      } // close CIJKL scope*/
     
    // we set elasptr data here
    elasptr->grdflag  = 0;
    elasptr->ndime    =  shpstrct->ndime;
    elasptr->disp     = &disp[0];
    elasptr->rcgreen  = rightcauchygreen;
    elasptr->stfparam = &KK[0];
    elasptr->secpk    = secondpiola;
    elasptr->cijklptr = CIJKL;

    
    
    // finally estiff assembly
    // refer equation 2.29 lengthy prospectus

    // first the CIJKL part
    elas2d(this->conseqnno,elasptr);
  
    for(int AA = 0;AA<shpstrct->numnodes;AA++){
      for(int ii = 0;ii<mdofn;ii++){
	int ieqn = mdofn*(AA-1+1) + ii;
	  
	  for(int BB=0;BB<shpstrct->numnodes;BB++){
	    for(int jj=0;jj<mdofn;jj++){
	      int jeqn = mdofn*(BB-1+1) + jj;
	      
	      double contrib=0;
	      for(int II=0;II<shpstrct->ndime;II++){
		double der1 = shpstrct->glbder[iinte][AA][II];

		for(int JJ=0;JJ<shpstrct->ndime;JJ++){
		  for(int KK=0;KK<shpstrct->ndime;KK++){
		    double der2 = shpstrct->glbder[iinte][BB][KK];

		    for(int LL=0;LL<shpstrct->ndime;LL++){
		      //estiff[ieqn][jeqn] += der1*defgrad[ii][JJ]*CIJKL[JJ][II][KK][LL]*defgrad[jj][LL]*der2*wtjac;
		      //contrib += defgrad[ii][JJ]*CIJKL[JJ][II][KK][LL]*defgrad[jj][LL];
		      contrib  += der1*defgrad[ii][JJ]*CIJKL[JJ][II][KK][LL]*defgrad[jj][LL]*der2;
		    } //LL
		    
		    // contrib = contrib*der2;
		  }//KK
		}//JJ

		//contrib = der1*contrib;
	      }//II
	      
	      estiff[ieqn][jeqn] += contrib*wtjac;

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
		
		estiff[ieqn][jeqn] += der1*(secondpiola[JJ][II]*delta(ii,kk))*der2*wtjac;
		
	      }
	    }

	  }
	}
      }
    }
    
    // finally the rhs
    // Hughes ME235 pg 28 Chapter 4.102

    for(int AA = 0;AA<shpstrct->numnodes;AA++){
      for(int ii =0;ii<shpstrct->ndime;ii++){
	int ieqn  = mdofn*(AA-1+1)+ii;

	for(int II = 0;II<shpstrct->ndime;II++){
	  double der1 = shpstrct->glbder[iinte][AA][II];
	
	  for(int JJ = 0;JJ<shpstrct->ndime;JJ++){
	    erhs[ieqn] = erhs[ieqn]-der1*secondpiola[JJ][II]*defgrad[ii][JJ]*wtjac;
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
  deallocmat(shpstrct->ndime,secondpiola);
  deallocmat(shpstrct->ndime,rightcauchygreeninvtpose);
  deallocmat(shpstrct->ndime,linstrn);
  dealloc4D(shpstrct->ndime,shpstrct->ndime,shpstrct->ndime,shpstrct->ndime,CIJKL);

  // free elasptr
  free(elasptr);
}

//=======================================================================
void nonlinelas3d::calcfunc(){
  // calculate the functional 
  // done with nonlinelas1d as a model
  
  // begin by zeroing the functional
  this->functional = 0;
  for(int iinte =0; iinte<shpstrct->tinte;iinte++){
    double jacodet = shpstrct->jdet[iinte];
    double weight  = shpstrct->wts[iinte];

    for(int idime = 0;idime < shpstrct->ndime ; idime++){
      double prim = 0;
      double meas = 0;
      for(int inode  = 0;inode < shpstrct->numnodes;inode++){
	prim += this->primal[inode][idime]*shpstrct->shape[iinte][inode];
	meas += measdata[inode][idime]*shpstrct->shape[iinte][inode];
      }
      this->functional += 0.5*(meas-prim)*(meas-prim)*jacodet*weight;
    }
    
  }
}

//=======================================================================
void nonlinelas3d:: elemrhs(){

}

//=======================================================================
void nonlinelas3d:: calcgradient(){
    // we make use of the constitutive equation library.
  // again heavily based on the nonlinelas1d::calcgradient()

  // this is the data which will go into conslib
  // declarations specific to nonlinear problems;
  double *KK;
  // material displacement vector (at an integration point)
  double *disp;  
  double **defgrad;
  double **gradu;
  double **secondpiola;
  double **rightcauchygreen;
  double **rightcauchygreeninvtpose;
  double **linstrn;
  double ****CIJKL;
  double **dualgrad; // material gradient of the dual field
  double **pk2grad;  // gradient of the second piola kirchoff wrt material properties

  // declare and allocate a pointer to a conselas1d structure
  conselas23d *elasptr; elasptr = (conselas23d *)malloc(sizeof(conselas23d));

  KK                  = allocvec(this->nprop,KK);
  disp                = allocvec(shpstrct->ndime,disp);
  defgrad             = allocmat(shpstrct->ndime,shpstrct->ndime,defgrad);
  gradu               = allocmat(shpstrct->ndime,shpstrct->ndime,gradu);
  linstrn             = allocmat(shpstrct->ndime,shpstrct->ndime,linstrn);
  secondpiola         = allocmat(shpstrct->ndime,shpstrct->ndime,secondpiola);
  rightcauchygreen    = allocmat(shpstrct->ndime,shpstrct->ndime,rightcauchygreen);

  CIJKL               = alloc4D(shpstrct->ndime,shpstrct->ndime,shpstrct->ndime,shpstrct->ndime,CIJKL);
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
    // see TJRH ME 235 Stanford, unauthorized notes.
    
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
    elasptr->rcgreen  = rightcauchygreen;
    elasptr->stfparam = &KK[0];
    elasptr->secpk    = secondpiola;
    elasptr->cijklptr = CIJKL;
    elasptr->pk2grad  = pk2grad;

    // in the next loop we set which "stiffness parameter" we are taking the gradient wrt.

    for(int iprop = 0;iprop <nprop;iprop++){

      elasptr->grdparm = iprop;
      elas2d(this->conseqnno,elasptr); // fills pk2grad accordingly.

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
	
      }// inode
    } // iprop

  }//iinte
  
  deallocvec(KK);
  deallocvec(disp);
  deallocmat(shpstrct->ndime,gradu);
  deallocmat(shpstrct->ndime,defgrad);
  deallocmat(shpstrct->ndime,rightcauchygreen);
  deallocmat(shpstrct->ndime,secondpiola);
  deallocmat(shpstrct->ndime,rightcauchygreeninvtpose);
  deallocmat(shpstrct->ndime,linstrn);
  dealloc4D(shpstrct->ndime,shpstrct->ndime,shpstrct->ndime,shpstrct->ndime,CIJKL);
  deallocmat(shpstrct->ndime,dualgrad);
  deallocmat(shpstrct->ndime,pk2grad);

  // free the elasptr structure
  free(elasptr);
}

//======================================================================
void nonlinelas3d::getinteg(){
  // get integration points
  gauss3d(shpstrct);
}
//====================================================================
void nonlinelas3d::getshape(){
  //cout<<"getshape in nonlinelas3d()"<<endl;
  shape3d(shpstrct);
}
//======================================================================
void nonlinelas3d::printeltype(){
  cout <<endl<<"This is NONLINELAS3D"<<endl;
}
//=====================================================================
void nonlinelas3d::memalloc(){
  // material properties, connectivity, stiffnessmatrix, righthand side
  // ideqn,idnode,isbc
  // The best way is to inherit this..
  
}
//=====================================================================
void nonlinelas3d::memdealloc(){
  // dealloc memory

}
//=====================================================================
void nonlinelas3d::setcoord(double cval,int idime,int inode){
  // cval = coordinate value
  // inode = local nodenumber
  // idime = dimension
  // Again should, be in base class
  shpstrct->coord[inode][idime] = cval;
}
//====================================================================
void nonlinelas3d::setprop(double propval,int iprop,int inode){
  // cval = coordinate value
  // inode = local nodenumber
  // idime = dimension
  // Again should, be in base class
  prop[inode][iprop] = propval;
}
//===================================================================
void nonlinelas3d::setideqn(int inode,int idofn,long int eqnno){

  assert(inode<shpstrct->numnodes);  // 
 
 
  assert(idofn<dofnatnode[inode]);   // 
  


  ideqn[inode][idofn] = eqnno;
}
//==================================================================
void nonlinelas3d::setidnode(int inode,long int globalnode){
  assert(inode<shpstrct->numnodes);
  idnode[inode]= globalnode;
}
//=================================================================
void nonlinelas3d::setisbc(int inode,int idofn,int value){
  assert(inode<shpstrct->numnodes);
  assert(idofn<dofnatnode[inode]);
  
  isbc[inode][idofn] = value;

}
//==============================================================
void nonlinelas3d::setdirdata(double ddata,int inode,int idofn){
   assert(inode<shpstrct->numnodes);
   assert(idofn<dofnatnode[inode]);
   
   dirdata[inode][idofn] = ddata;
}
//==============================================================
void nonlinelas3d::settracdata(double tdata,int inode,int idofn){
  assert(inode<shpstrct->numnodes);
  assert(idofn<dofnatnode[inode]);

  tracdata[inode][idofn] = tdata;
}
//=============================================================
//=============================================================
void nonlinelas3d::setmeasdata(double mdata,int inode,int idofn){
  assert(inode<shpstrct->numnodes);
  assert(idofn<dofnatnode[inode]);
  
  measdata[inode][idofn] = mdata;

}
//==============================================================
void nonlinelas3d::calcdual(){
  // the whole point is to calculate the forcing for the dualfield
  // again, closely modelled on nonlinelas1d::calcdual;

  // begin by zeroing rhsvector
  zerovec(shpstrct->numnodes*mdofn,erhs);
  
  // loop over each integration point
  
  for(int iinte = 0;iinte<shpstrct->tinte;iinte++){
    double jacodet = shpstrct->jdet[iinte];
    double weight  = shpstrct->wts[iinte];

    // pick dimension to force
    for(int jdime = 0;jdime<shpstrct->ndime;jdime++){
      for(int jnode =0;jnode <shpstrct->numnodes;jnode++){

	for(int idime = 0;idime<shpstrct->ndime;idime++){
	  double prim = 0;
	  double meas = 0;
	  // interpolate primal and dual at this point
	  for(int inode = 0;inode<shpstrct->numnodes;inode++){
	    prim += primal[inode][idime]*shpstrct->shape[iinte][inode];
	    meas += measdata[inode][idime]*shpstrct->shape[iinte][inode];
	  } // inode
	  // interpolation complete, now force the dual
	  if(jdime == idime){
	    int ieqn = mdofn*(jnode - 1 + 1) + idime;

	  this->erhs[ieqn] += -(prim-meas)*shpstrct->shape[iinte][jnode]*jacodet*weight;
	  }
	}//idime
	
      } // jnode
    } // jdime
    
  }
  
}
//==============================================================
void nonlinelas3d::setdualfield(double ddata,int inode,int idofn){
  assert(inode<shpstrct->numnodes);
  assert(idofn<dofnatnode[inode]);
  
  dualfield[inode][idofn] = ddata;
}

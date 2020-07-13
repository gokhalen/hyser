/*======================================================================
 *  Nachiket Gokhale gokhalen@bu.edu
 *
 *  Two dimensional Linear Elasticity Element
 *  Builds stiffness matrix and right hand side
 *
 **========================================================================*/

/* selective integration -- used in the rhs for initial strain, standard stiffness matrix, and dirichlet rhs ( automatic). If a change is made, then it should be made in all places for consistency.*/


#ifndef HAVE_PIEZOREDUCED2D_HPP
#include "piezoreduced2d.hpp"
#endif

#include <cassert>
#include <cstdlib>

using namespace std;


// implement constrctor 

piezoreduced2d::piezoreduced2d(int *cinte, int cmregparm, long int celemno){
  nprop = 8;
  nregprop = 2;
  mdofn = 3;
  nonlinflag = false;
  this->mregparm = cmregparm;
  
  uscal2 = 1.0;
  vscal2 = 1.0;
  uscal = sqrt(uscal2);
  vscal = sqrt(vscal2);
  uvscal = uscal*vscal;

  shpstrct = new shapestruct(2,2,cinte,4,PIEZOREDUCED2D,celemno);
  
  // allocate the rest after number of nodes has been set
  prop        = allocmat(shpstrct->numnodes,nprop,prop); 
  estiff      = allocmat(shpstrct->numnodes*mdofn,shpstrct->numnodes*mdofn,estiff);
  estiffimag   = allocmat(shpstrct->numnodes*mdofn,shpstrct->numnodes*mdofn,estiffimag);
  ecmass      = allocmat(shpstrct->numnodes*mdofn,shpstrct->numnodes*mdofn,ecmass);
  dofnatnode  = allocvec(shpstrct->numnodes,dofnatnode);
  density     = allocvec(shpstrct->numnodes,density);

  // set the dofnatnode array, one degree of freedom at each node
  for(int inode = 0;inode<shpstrct->numnodes;inode++){
    dofnatnode[inode] = mdofn;
  }

  // allocate primal
  primal  = allocmat(shpstrct->numnodes,mdofn,primal);
  zeromat(shpstrct->numnodes,mdofn,primal);

  lastupdate  = allocmat(shpstrct->numnodes,mdofn,lastupdate);
  zeromat(shpstrct->numnodes,mdofn,lastupdate);

   // allocate an array for initial strain vector
  initstrn = allocvec(3*(shpstrct->ndime-1) , initstrn);

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

  // allocate primal field data
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
  
  memalloc();
}
//==================================================================
piezoreduced2d::~piezoreduced2d(){
  deallocmat(shpstrct->numnodes,primal);
  deallocmat(shpstrct->numnodes,lastupdate);
  deallocvec(initstrn);
  deallocmat(shpstrct->numnodes,prop);
  deallocmat(shpstrct->numnodes*mdofn,estiff);
  deallocmat(shpstrct->numnodes*mdofn,estiffimag);
  deallocmat(shpstrct->numnodes*mdofn,ecmass);
  deallocvec(dofnatnode);
  deallocvec(density);
  deallocmat(shpstrct->numnodes,dirdata);
  deallocmat(shpstrct->numnodes,tracdata);
  deallocvec(erhs);

  deallocmat(shpstrct->numnodes,ideqn);
  deallocvec(idnode);
  deallocmat(shpstrct->numnodes,isbc);
  //shpstrct->~shapestruct();
  //deallocate measdata
  deallocmat(shpstrct->numnodes,measdata);


  // deallocate the whole fields
  dealloc3D(shpstrct->numnodes,mdofn,primalwhole);
  dealloc3D(shpstrct->numnodes,mdofn,dualwhole);

  // deallocate dualfield
  deallocmat(shpstrct->numnodes,dualfield);
  // deallocate gradient 
  deallocmat(shpstrct->numnodes,gradient);
  //deallocvec(regparm);
  deallocmat(shpstrct->ndime,functens);
  deallocmat(shpstrct->numnodes,Tprim);
  deallocmat(shpstrct->numnodes,Tmeas);
  
  memdealloc();
  delete shpstrct;
}
//================================================================

//===============================================================
void piezoreduced2d::elemstiff(int which_stiffness){
  // zero element stiffness K_{A,B} = \int_{\Omega}(N_{A,x}k(x)N_{B,X})
  zeromat(shpstrct->numnodes*mdofn,shpstrct->numnodes*mdofn,estiff);
  zeromat(shpstrct->numnodes*mdofn,shpstrct->numnodes*mdofn,estiffimag);

  int tempninte[2],temptinte;
  tempninte[0] =shpstrct->ninte[0];
  tempninte[1] =shpstrct->ninte[1];
  temptinte    =shpstrct->tinte;
  //pick property to implement reduced integration
  // if iprop == 1 then we do selective integration on lambda 
  for (int iprop = 0; iprop < 2; iprop++){

   if(iprop == 0){
      shpstrct->tinte    = 1;
      shpstrct->ninte[0] = 1;
      gauss2d(shpstrct);
      shape2d(shpstrct);
    }
   else
    {
      //cout << "reducecd integration stopped"<<endl;
      shpstrct->tinte    = temptinte;
      shpstrct->ninte[0] = tempninte[0];
      shpstrct->ninte[1] = tempninte[1];
      gauss2d(shpstrct);
      shape2d(shpstrct);
    }
    
    //loop over integration points
    for(int iinte=0;iinte<shpstrct->tinte;iinte++){
      double jacodet = shpstrct->jdet[iinte];
      double weight  = shpstrct->wts[iinte];
      //cout<<__FILE__<<"WTJAC =="<<shpstrct->jdet[0]*shpstrct->wts[0];
      double lambda = 0; double mu = 0;
      double ee31   = 0.0; double ee15 = 0.0; double ee33=0.0;
      double kk11   = 0.0; double kk22 = 0.0;
						
      
      // interpolate \lambda and \mu at this point
      for(int inode = 0;inode<shpstrct->numnodes;inode++){
	lambda += shpstrct->shape[iinte][inode]*prop[inode][0];
	mu     += shpstrct->shape[iinte][inode]*prop[inode][1];
	ee31   += shpstrct->shape[iinte][inode]*prop[inode][3];
	ee33   += shpstrct->shape[iinte][inode]*prop[inode][4];
	ee15   += shpstrct->shape[iinte][inode]*prop[inode][5];
	kk11   += shpstrct->shape[iinte][inode]*prop[inode][6];
	kk22   += shpstrct->shape[iinte][inode]*prop[inode][7];
      }// close interpolation of \lambda & \mu


    
      // cout<<__FILE__<<" lambda mu "<<lambda<<" "<<mu;
      // set up D matrix
      double DD[3][3];   // elastic
      double EET[2][3];   // piezoelectric properties
      double EE[3][2];
      double KK[2][2];   // dielectric 
      
    
      if(iprop == 0){
	DD[0][0] = lambda/(uscal2);        DD[0][1] = lambda/uscal2;        DD[0][2] = 0;
	DD[1][0] = lambda/(uscal2);        DD[1][1] = lambda/uscal2;        DD[1][2] = 0;
	DD[2][0] = 0;             DD[2][1] = 0;             DD[2][2] = 0;
      }
      else{
	DD[0][0] = 2*mu/uscal2;          DD[0][1] = 0;         DD[0][2] = 0;
	DD[1][0] = 0;             DD[1][1] = 2*mu/uscal2;      DD[1][2] = 0;
	DD[2][0] = 0;             DD[2][1] = 0;                DD[2][2] = mu/uscal2;
      }

      // pz 
      {
	EET[0][0] = 0.0;            EET[0][1] = 0.0;         EET[0][2] = ee15/uvscal; 
	EET[1][0] = ee31/(uvscal); EET[1][1] = ee33/uvscal; EET[1][2] = 0.0; 
	
	for (int itmp = 0; itmp < 3; itmp++){
	  for (int jtmp = 0; jtmp < 2; jtmp++){
	    EE[itmp][jtmp] = EET[jtmp][itmp];
	  }
	}

      }

      // dielect
      {
	KK[0][0] = kk11/vscal2;  KK[0][1] = 0.0;
	KK[1][0] = 0.0;          KK[1][1] = kk22/vscal2;
      }
    
      for(int inode = 0; inode<shpstrct->numnodes; inode++){
	// set up B_{a} and initialize it
	double BA[3][2];
	double BAKK[2];  // dielectric;

	double kvecx = this->kvec[0];
	double kvecy = this->kvec[1];
      
	// standard B-matrices
	BA[0][0]   = shpstrct->glbder[iinte][inode][0]; 
	BA[1][0]   = 0;
	BA[2][0]   = shpstrct->glbder[iinte][inode][1];
      
	BA[0][1]   = 0;
	BA[1][1]   = shpstrct->glbder[iinte][inode][1];
	BA[2][1]   = shpstrct->glbder[iinte][inode][0];

	
	// dielectric matrix
	BAKK[0] = shpstrct->glbder[iinte][inode][0];
	BAKK[1] = shpstrct->glbder[iinte][inode][1];
      
      
	for(int jnode =0; jnode<shpstrct->numnodes; jnode++){
	  double BB[3][2];
	  double BBKK[2];
	  
	  
	  // standard B mats
	  BB[0][0]   = shpstrct->glbder[iinte][jnode][0];
	  BB[1][0]   = 0;
	  BB[2][0]   = shpstrct->glbder[iinte][jnode][1];
	
	  BB[0][1]   = 0;
	  BB[1][1]   = shpstrct->glbder[iinte][jnode][1];
	  BB[2][1]   = shpstrct->glbder[iinte][jnode][0];

	  // dielectric matrix
	  BBKK[0] = shpstrct->glbder[iinte][jnode][0];
	  BBKK[1] = shpstrct->glbder[iinte][jnode][1];

	
	  for(int idofn = 0;idofn<dofnatnode[inode];idofn++){
	    // need the (inode  - 1 + 1) because of zero based indexing
	    // not correct: assumes dofnatnode = mdofn 
	    // use mdofn instead
	    for(int jdofn = 0;jdofn<dofnatnode[jnode];jdofn++){
	    
	      // elasticity matrix 
	      if (( idofn < 2 ) && ( jdofn < 2 )){
		int ieqn = (mdofn*(inode-1 + 1) + idofn);
		int jeqn = (mdofn*(jnode-1 + 1) + jdofn);

		for(int ii=2;ii>=0;--ii){
		  for(int jj=2;jj>=0;--jj){
		    estiff[ieqn][jeqn] +=   BA[ii][idofn]*DD[ii][jj]*BB[jj][jdofn]*jacodet*weight;
		  }//jj
		}//ii
		
	      }// closes elasticity

	      if ( iprop == 1 ){
		// dielectric
		if ( (idofn == 2 ) && ( jdofn == 2 )){
		  int ieqn = mdofn*(inode-1 + 1) + idofn; 
		  int jeqn = mdofn*(jnode-1 + 1) + jdofn;
		  
		  for ( int ii = 0; ii < 2 ; ii++){
		    for ( int jj = 0; jj < 2; jj++){
		      estiff[ieqn][jeqn] -=BAKK[ii]*KK[ii][jj]*BBKK[jj]*jacodet*weight;  
		    }//jj
		  }//ii
		  
		} // closes dielectric
		
		// piezo-electric 
		if ( ( idofn != 2 )  && ( jdofn == 2 ) ){
		  int ieqn = mdofn*(inode-1 + 1) + idofn; 
		  int jeqn = mdofn*(jnode-1 + 1) + jdofn;
		  
		  for ( int ii = 0; ii < 3; ii++){
		    for ( int jj = 0; jj < 2; jj++){
		      estiff[ieqn][jeqn] += BA[ii][idofn]*EE[ii][jj]*BBKK[jj]*jacodet*weight;
		      // transpose
		      estiff[jeqn][ieqn] += estiff[ieqn][jeqn];
		    }
		  }
		  
		}//closes pz
	      }// closes iprop
	      
		 
	    
	    }// jdofn
	  } // idofn
	

	} // closes jnode
      
      }//closes inode
    
    } // closes iinte
  }//iprop

  //cout<<endl;

  //cout <<"estiff"<<endl;
  //for(int ii = 0;ii<12;ii++){
  //   for(int jj =0;jj<12;jj++){
  //     cout<<" "<<estiff[ii][jj]<<" ";
  //   }
  //  cout<<endl;
  //}
  

}//closes elemstiff
//====================================================================
void piezoreduced2d::elemcmass(int which_mass){
  zeromat(shpstrct->numnodes*mdofn,shpstrct->numnodes*mdofn,ecmass);
    //loop over integration points

    for(int iinte=0;iinte<shpstrct->tinte;iinte++){
      double jacodet = shpstrct->jdet[iinte];
      double weight  = shpstrct->wts[iinte];
      
      // interpolate - density
      double rhotmp = 0.0;
      for(int inode = 0;inode<shpstrct->numnodes;inode++){
	rhotmp +=  this->density[inode]*shpstrct->shape[iinte][inode];
      }

      for(int inode = 0;inode<shpstrct->numnodes;inode++){
	for(int idofn = 0;idofn<dofnatnode[inode]; idofn++){
	    int ieqn = (dofnatnode[inode]*(inode-1 + 1) + idofn);

	    for(int jnode = 0;jnode<shpstrct->numnodes;jnode++){
	      for(int jdofn = 0;jdofn<dofnatnode[inode]; jdofn++){
		int jeqn = (dofnatnode[jnode]*(jnode-1 + 1) + jdofn);
		//cout <<"ieqn = "<<ieqn<<" jeqn= "<<jeqn<<endl;
		ecmass[ieqn][jeqn] += delta(idofn,jdofn)*rhotmp*shpstrct->shape[iinte][inode]*shpstrct->shape[iinte][jnode]*jacodet*weight;
		//ecmass[ieqn][jeqn] = delta(ieqn,jeqn);
	      } //jdofn
	    } // jnode

	}//idofn
      } // inode
      
    }//close integ loop

    // lump if required
    if ( which_mass == BUILD_LUMPED_MASS ){

      for(int inode = 0;inode<shpstrct->numnodes;inode++){
	for(int idofn = 0;idofn<dofnatnode[inode]; idofn++){
	  int ieqn = (dofnatnode[inode]*(inode-1 + 1) + idofn);

	  for(int jnode = 0;jnode<shpstrct->numnodes;jnode++){
	    for(int jdofn = 0;jdofn<dofnatnode[inode]; jdofn++){
	      int jeqn = (dofnatnode[jnode]*(jnode-1 + 1) + jdofn);
	      // put off diagonal mass into the diagonal
	      ecmass[ieqn][ieqn] += (1-delta(ieqn,jeqn))*ecmass[ieqn][jeqn];
	      // zero out non-diagonal entries
	      ecmass[ieqn][jeqn]  = delta(ieqn,jeqn)*ecmass[ieqn][jeqn];
	    }//jdofn
	  }//jnode

	}//idofn
      }//inode 

    }

} // closes elemcmass
      
      

//void piezoreduced2d::calcfunc(){
//
//}
//====================================================================
void piezoreduced2d::elemrhs(){
  // zero right hand side vector
  zerovec(shpstrct->numnodes*mdofn,erhs);
   
  // elemstiff should be call from outside before calling elemrhs
  // dirichlet data: take product of matrix with right hand side
  
  for(int inode = 0;inode<shpstrct->numnodes;inode++){
    for(int idofn = 0;idofn<dofnatnode[inode];idofn++){
      
      // zero based indexing. therefore, -1 + 1 
      int ieqn = (inode-1 + 1)*mdofn + idofn;
      //take matrix-vector product

      for(int jnode=0;jnode<shpstrct->numnodes;jnode++){
	for(int jdofn=0;jdofn<dofnatnode[jnode];jdofn++){
	  // zero based indexing. therefore, -1 + 1
	  int jeqn    = (jnode-1+1)*dofnatnode[jnode] + jdofn;
	  if ( jdofn < 2 ) {
	    erhs[ieqn] -=  estiff[ieqn][jeqn]*dirdata[jnode][jdofn];
	  }
	  else{
	    erhs[ieqn] -=  estiff[ieqn][jeqn]*dirdata[jnode][jdofn];
	  }
	} // jnode
      }   // jdofn
    }     // idofn
  }       // inode

  //for(int jj =0;jj<8;jj++){
  // cout<<__FILE__<<__LINE__<<"erhs "<<erhs[jj]<<" "<<endl;
  //}
 
    
}//elemrhs()
//====================================================================
void piezoreduced2d::calcgradient(){

  // the implementation follows the \b{w}^{T}\b{B}^{T}{D}\b{B}^{T}\b{\epsilon} 
  // representation of the weak form of the elasticity equations

  // see Linear Hughes ()
  // this treats \lambda and \mu as two sep
  // optimization variables
  

  // loop over integration points


  double dispvec[8];
  double dualvec[8];

  // initialize the displacement and dual vectors

  {
    int ictr = 0;
    for (int inode = 0; inode < shpstrct->numnodes; inode++){
      for (int idofn = 0; idofn <dofnatnode[inode]; idofn++){
	//cout <<__FILE__<<__LINE__<<"ictr== "<<ictr<<endl;
	dispvec[ictr] = this->primal[inode][idofn];
	dualvec[ictr] = this->dualfield[inode][idofn];
	ictr++;
      }
    }
  }


  for(int iinte = 0; iinte <shpstrct->tinte; iinte++){
    double jacodet  = shpstrct->jdet[iinte];
    double weight   = shpstrct->wts[iinte];
    
    for(int iprop = 0; iprop <this->nprop; iprop++){
      // this is the property to take the gradient of 
      
      double DD[3][3];
      
      for(int knode = 0; knode <shpstrct->numnodes; knode++){
	// this is the node at which_we take the gradient 
	double gradshp = shpstrct->shape[iinte][knode];

	if(iprop ==0){
	  DD[0][0] = gradshp;   DD[0][1] = gradshp;        DD[0][2] = 0;
	  DD[1][0] = gradshp;   DD[1][1] = gradshp;        DD[1][2] = 0;
	  DD[2][0] = 0;         DD[2][1] = 0;              DD[2][2] = 0;
	}
	else{
	  DD[0][0] = 2*gradshp;   DD[0][1] = 0;              DD[0][2] = 0;
	  DD[1][0] = 0;           DD[1][1] = 2*gradshp;      DD[1][2] = 0;
	  DD[2][0] = 0;           DD[2][1] = 0;              DD[2][2] = gradshp;
	}

	for (int inode = 0; inode < shpstrct->numnodes; inode++){
	 
	  double BA[3][2];
	  
	  
	  BA[0][0]   = shpstrct->glbder[iinte][inode][0]; 
	  BA[1][0]   = 0;
	  BA[2][0]   = shpstrct->glbder[iinte][inode][1];
	  
	  BA[0][1]   = 0;
	  BA[1][1]   = shpstrct->glbder[iinte][inode][1];
	  BA[2][1]   = shpstrct->glbder[iinte][inode][0];

	  for (int jnode = 0; jnode <shpstrct->numnodes; jnode++){
	    double BB[3][2];

	    BB[0][0]   = shpstrct->glbder[iinte][jnode][0];
	    BB[1][0]   = 0;
	    BB[2][0]   = shpstrct->glbder[iinte][jnode][1];
	    
	    BB[0][1]   = 0;
	    BB[1][1]   = shpstrct->glbder[iinte][jnode][1];
	    BB[2][1]   = shpstrct->glbder[iinte][jnode][0];
	    
	    for (int idofn = 0; idofn < mdofn ; idofn++){
	      // the index for  the dual field
	      int ictr = mdofn*(inode-1+1) + idofn;
	      for (int jdofn = 0; jdofn < mdofn; jdofn++){
		// the index for the primal field
		int jctr = mdofn*(jnode-1+1) + jdofn;
		
		for(int ii = 0; ii<3;ii++){
		  for(int jj = 0; jj<3;jj++){
		    this->gradient[knode][iprop] += dualvec[ictr]*BA[ii][idofn]*DD[ii][jj]*BB[jj][jdofn]*dispvec[jctr]*weight*jacodet;
		  }
		}
		
	      }
	    }
	  }
	}

	// add the effects of regularization
	// the regularization code is the same for elements
	// what should happen here is that this element
	// this code should be placed in the base class and 
	// the easiest way this can be done is by creating a function 
	// called add regularization in the base class


      }//knode
    } //iprop
  }//iinte
  
}
//====================================================================
void piezoreduced2d::getinteg(){
  gauss2d(shpstrct);
}
//====================================================================
void piezoreduced2d::getshape(){
  shape2d(shpstrct);
}
//====================================================================
void piezoreduced2d::printeltype(){
  //cout<<endl<<"This is PIEZOREDUCED2D"<<endl;
}
//====================================================================
//void piezoreduced2d::memalloc(){

  // inherited
//}
//====================================================================
//void piezoreduced2d::memdealloc(){
//  // dealloc memory
//}
//===================================================================
void piezoreduced2d::eleminitstrn(){
  // zero the rhs
  zerovec(shpstrct->numnodes*mdofn,erhs);

  if ( initstrntype > 0 ){
    int pp, qq;
    getkl(&pp,&qq);
  
    for ( int iinte = 0; iinte < shpstrct->tinte; iinte++){
      
      double jacodet = shpstrct->jdet[iinte];
      double weight  = shpstrct->wts[iinte];
      double lambda  = 0.0;
      double mu      = 0.0;
      double ee31    = 0.0;
      double ee33    = 0.0;
      double ee15    = 0.0;
      double kk11    = 0.0;
      double kk22    = 0.0;

      // interpolation
      for(int inode = 0;inode<shpstrct->numnodes;inode++){
	lambda += shpstrct->shape[iinte][inode]*prop[inode][0];
	mu     += shpstrct->shape[iinte][inode]*prop[inode][1];
	ee31   += shpstrct->shape[iinte][inode]*prop[inode][3];
	ee33   += shpstrct->shape[iinte][inode]*prop[inode][4];
	ee15   += shpstrct->shape[iinte][inode]*prop[inode][5];
	kk11   += shpstrct->shape[iinte][inode]*prop[inode][6];
	kk22   += shpstrct->shape[iinte][inode]*prop[inode][7];
      }
    
      for ( int inode = 0; inode < shpstrct->numnodes; inode++){
	for ( int idofn = 0; idofn < mdofn; idofn++){
	  int ieqn = (mdofn*(inode-1 + 1) + idofn);
	
	  // first problem 
	  if  ( initstrntype <= 4 )  {
	  
	    if ( idofn < 2 ) {
	      double cijpqj = 0.0;

	      for ( int jdofn = 0; jdofn < 2; jdofn++){
		double cijpq = getelas(lambda,mu,idofn,jdofn,pp,qq);
		erhs[ieqn] -= shpstrct->glbder[iinte][inode][jdofn]*(cijpq/uscal)*jacodet*weight;
	      }

	  
	    }//idofn < 2
	    if ( idofn == 2 ){
	      // pz forc
	      for ( int jdofn = 0; jdofn < 2; jdofn++){
		for ( int kk = 0 ; kk < 2 ; kk++){
		  for ( int ll = 0; ll < 2; ll++){
		    double ejkl = getpz(jdofn,kk,ll,ee15,ee33,ee31);
		    erhs[ieqn] -= shpstrct->glbder[iinte][inode][jdofn]*(ejkl/vscal)*delta(kk,pp)*delta(ll,qq)*jacodet*weight;
		  }
		}
	      }
	    }//idofn ==2 

	  }//end first problem
	  // second problem
	  if ( initstrntype > 4 ) {
	    if ( idofn < 2 ){
	      // elas
	      for ( int jdofn = 0 ; jdofn < 2; jdofn++){
		for ( int kk = 0; kk < 2; kk++){
		  double ekij = getpz(kk,idofn,jdofn,ee15,ee33,ee31);
		  erhs[ieqn] -= shpstrct->glbder[iinte][inode][jdofn]*(ekij/uscal)*delta(kk,pp)*jacodet*weight;
		}
	      }
	    }
	    if ( idofn == 2 ){
	      // pz
	      for ( int jdofn = 0 ; jdofn < 2; jdofn++){
		for ( int mm = 0; mm < 2; mm++){
		  double kjm = getkk(jdofn,mm,kk11,kk22);
		  erhs[ieqn] += shpstrct->glbder[iinte][inode][jdofn]*(kjm/vscal)*delta(pp,mm)*jacodet*weight;
		}
	      }

	    }
	  }

	}//idofn
      }//inode
    }//iinte
  }//initstrn
  //if ( initstrntype == 4 ){
  //  for (int ii = 0; ii< 12; ii++){
  //    cout<<"ii= "<<ii<<" rhs= "<<erhs[ii]<<endl;
  //  }
  // }
  //cout<<"------------------------------------";
}//elemrhs()
// ====================================================================
void piezoreduced2d::homogenizehassanihinton(){
   zero4D(2,2,2,2,this->BASE_EHIJKL);
   zero3D(2,2,2,this->BASE_PZHKIJ);
   zeromat(2,2,this->BASE_KHIJ);

   if ( initstrntype !=0 ) {
   if ( initstrntype > 6 ) {
     cout<<endl<<"INVALID VALUE FOR initstrntype in elem_base.cpp"<<endl;
     exit(1);
   }
   // select field, which defines the indices over which we sum
   int kk; int ll;
   getkl(&kk,&ll);

   // rescale primal 
   for(int inode = 0 ; inode < 4; inode++){
     this->primal[inode][0] = this->primal[inode][0]/uscal;
     this->primal[inode][1] = this->primal[inode][1]/uscal;
     this->primal[inode][2] = this->primal[inode][2]/vscal;
   }

   
   for ( int iinte = 0; iinte < shpstrct->tinte; iinte++){

     double jacodet = shpstrct->jdet[iinte];
     double weight  = shpstrct->wts[iinte];
     double lambda  = 0.0;
     double mu      = 0.0;
     double ee31   = 0.0;
     double ee33   = 0.0;
     double ee15   = 0.0;
     double kk11   = 0.0;
     double kk22   = 0.0;

     double chikl[2][2]; chikl[0][0] = 0.0;chikl[0][1] = 0.0;chikl[1][0] = 0.0;chikl[1][1] = 0.0;
     double psikl[2]; psikl[0]=0.0; psikl[1]=0.0; 
     double phik[2][2]; phik[0][0] = 0.0;phik[0][1] = 0.0;phik[1][0] = 0.0;phik[1][1] = 0.0;
     double Rk[2]; Rk[0] = 0.0; Rk[1] = 0.0;

     // interpolation
     for(int inode = 0;inode<shpstrct->numnodes;inode++){
       lambda += shpstrct->shape[iinte][inode]*prop[inode][0];
       mu     += shpstrct->shape[iinte][inode]*prop[inode][1];
       ee31   += shpstrct->shape[iinte][inode]*prop[inode][3];
       ee33   += shpstrct->shape[iinte][inode]*prop[inode][4];
       ee15   += shpstrct->shape[iinte][inode]*prop[inode][5];
       kk11   += shpstrct->shape[iinte][inode]*prop[inode][6];
       kk22   += shpstrct->shape[iinte][inode]*prop[inode][7];
     }


     // effective elastic tensor
     if ( initstrntype <= 4 ) {
       // get grads of inf. func. the kl index is chosen implicitly thru localization of primal
       for ( int pp = 0; pp < 2; pp++){
	 for ( int qq = 0; qq < 2; qq++){
	   for (int inode = 0; inode < shpstrct->numnodes; inode++){
	     chikl[pp][qq] += 0.5*( shpstrct->glbder[iinte][inode][qq]*this->primal[inode][pp] + shpstrct->glbder[iinte][inode][pp]*this->primal[inode][qq]); 
	   }
	   
	 }//qq

	 for (int inode = 0; inode < shpstrct->numnodes; inode++){
	   psikl[pp] += shpstrct->glbder[iinte][inode][pp]*this->primal[inode][2];
	 }

       }//pp
       
       for ( int ii = 0; ii < 2 ; ii++){
	 for ( int jj = 0; jj < 2; jj++){

	   double cav    = getelas(lambda,mu,ii,jj,kk,ll);

	   double corr1 = 0.0;
	   double corr2 = 0.0;
	   // compute corrections
	   
	   for ( int pp = 0; pp < 2; pp++){
	     for ( int qq = 0; qq < 2; qq++){
	       double cijpq = getelas(lambda,mu,ii,jj,pp,qq);
	       corr1 += chikl[pp][qq]*cijpq;
	     }
	   }	     
		     
	   for ( int pp = 0; pp < 2; pp++){
	     double epij = getpz(pp,ii,jj,ee15,ee33,ee31);
	     corr2 += epij*psikl[pp];
	   }
	   
	   double cijkl = cav + corr1 + corr2;

	   this->BASE_EHIJKL[ii][jj][kk][ll] += (1.0/this->unitcellvolume)*cijkl*jacodet*weight;
	 } //jj
       } //ii
     }//end elastic
     
     if ( initstrntype > 4 ){
       // compute influence functions
       for(int pp = 0 ; pp < 2 ; pp++){
	 for(int qq = 0 ; qq < 2 ; qq++){
	   for(int inode=0; inode<shpstrct->numnodes; inode++){
	     phik[pp][qq] += 0.5*( shpstrct->glbder[iinte][inode][qq]*this->primal[inode][pp] + shpstrct->glbder[iinte][inode][pp]*this->primal[inode][qq]); 
	   }//inode
	 }//qq
	 for(int inode=0; inode<shpstrct->numnodes; inode++){
	   Rk[pp] += shpstrct->glbder[iinte][inode][pp]*this->primal[inode][2];
	 }//inode
       }//pp

       // effective piezo tensor
       for ( int ii = 0; ii < shpstrct->ndime ; ii++){
	 for ( int jj = 0; jj < shpstrct->ndime; jj++){
	   double ekijav  = getpz(kk,ii,jj,ee15,ee33,ee31);
	   double corr1 = 0.0;
	   double corr2 = 0.0;
	   for (int pp = 0; pp < shpstrct->ndime; pp++){
	     for (int qq = 0; qq < shpstrct->ndime; qq++){
	        corr1 += getelas(lambda,mu,ii,jj,pp,qq)*phik[pp][qq] ;
	     }
	   }
	   for (int pp = 0; pp < shpstrct->ndime; pp++){
	     double epij = getpz(pp,ii,jj,ee15,ee33,ee31);
	     corr2 +=epij*Rk[pp];
	   }
	   
	   double ekij = ekijav + corr1 + corr2;
	   this->BASE_PZHKIJ[kk][ii][jj] += (1.0/this->unitcellvolume)*ekij*jacodet*weight;

	 } // ii
       } //jj

       // effective dielectric tensor
	 for ( int ii = 0; ii < shpstrct->ndime; ii++){
	   double kikav = getkk(ii,kk,kk11,kk22);
	   double corr1 = 0.0;
	   double corr2 = 0.0;

	   for ( int pp = 0; pp < 2; pp++){
	     for ( int qq = 0; qq < 2; qq++){
	       double eipq = getpz(ii,pp,qq,ee15,ee33,ee31);
	       corr1 += eipq*phik[pp][qq];
	     }
	   }

	   for ( int pp = 0; pp < 2; pp++){
	     //	     corr2 +=kk11*delta(ii,pp)*Rk[pp];
	     corr2 += getkk(ii,pp,kk11,kk22)*Rk[pp];
	   }
	   
	   double kik = kikav -corr1 + corr2;
	   this->BASE_KHIJ[ii][kk] += (1.0/this->unitcellvolume)*kik*jacodet*weight;
	 }
     
       
     }// end initstrntype > 4 (i.e. pz & props )
   }//iinte
   } //initstrntype !=0
}
// =========================================================================================
void piezoreduced2d::computeaveprops(){
  
  zero4D(shpstrct->ndime,shpstrct->ndime,shpstrct->ndime,shpstrct->ndime,this->BASE_EHAVIJKL);
  zero3D(shpstrct->ndime,shpstrct->ndime,shpstrct->ndime,this->BASE_PZHAVKIJ);
  zeromat(shpstrct->ndime,shpstrct->ndime,this->BASE_KHAVIJ);

  for(int iinte=0; iinte < shpstrct->tinte; iinte++){
    double jacodet = shpstrct->jdet[iinte];
    double weight  = shpstrct->wts[iinte];
      
    double lambda = 0.0;
    double mu     = 0.0;
    double ee31   = 0.0;
    double ee33   = 0.0;
    double ee15   = 0.0;
    double kk11   = 0.0;
    double kk22   = 0.0;

    for(int inode = 0;inode<shpstrct->numnodes;inode++){
      lambda += shpstrct->shape[iinte][inode]*prop[inode][0]*jacodet*weight;
      mu     += shpstrct->shape[iinte][inode]*prop[inode][1]*jacodet*weight;
      ee31   += shpstrct->shape[iinte][inode]*prop[inode][3]*jacodet*weight;
      ee33   += shpstrct->shape[iinte][inode]*prop[inode][4]*jacodet*weight;
      ee15   += shpstrct->shape[iinte][inode]*prop[inode][5]*jacodet*weight;
      kk11   += shpstrct->shape[iinte][inode]*prop[inode][6]*jacodet*weight;
      kk22   += shpstrct->shape[iinte][inode]*prop[inode][7]*jacodet*weight;
    }// close interpolation of \lambda & \mu
    
    // set CIJKL and PZKIJ ave
    for (int ii = 0; ii < shpstrct->ndime; ii++){
      for (int  jj = 0; jj < shpstrct->ndime; jj++){
	for (int kk = 0; kk < shpstrct->ndime; kk++){
	  for (int ll = 0; ll < shpstrct->ndime; ll++){    
	    //double first = delta(ii,jj)*delta(kk,ll);
	    //double second = delta(ii,kk)*delta(jj,ll) + delta(ii,ll)*delta(jj,kk);
	    //double temp = lambda*first + mu*second;
	    this->BASE_EHAVIJKL[ii][jj][kk][ll] += getelas(lambda,mu,ii,jj,kk,ll);
	  }
	}
      }
    }
    
    // e11 = e111 = ee[1][1][1]
    // e31 = e311 = ee[2][1][1]
    // e13 = e133 = ee[1][2][2]
    // e33 = e333 = ee[2][2][2]
    // e15 = e113 = ee[1][1][2] = ee[1][2][1]   
    // e35 = e313 = ee[2][1][2] = ee[2][2][1]

    //this->BASE_PZHAVKIJ[1][1][1] += ee11;
    this->BASE_PZHAVKIJ[1][0][0] += ee31;
    //this->BASE_PZHAVKIJ[1][2][2] += ee13;
    this->BASE_PZHAVKIJ[1][1][1] += ee33;
    this->BASE_PZHAVKIJ[0][0][1] += ee15;
    this->BASE_PZHAVKIJ[0][1][0] += ee15;
    //this->BASE_PZHAVKIJ[2][1][2] += ee35;
    //this->BASE_PZHAVKIJ[2][2][1] += ee35;
    
    
    // compute the average piezoelectric tensor
    this->BASE_KHAVIJ[0][0] += kk11;
    this->BASE_KHAVIJ[1][1] += kk22;
    
  } //iinte
  
}


//==================================================================
double piezoreduced2d::getpz(int kk,int ii,int jj,double ee15,double ee33,double ee31){
  double ekij = 0.0;
  if ( ( kk == 0 ) && ( ii != jj )){
    ekij =  ee15;
  }
  if ( (kk == 1) && (ii == jj )){
    if ( ii == 0 ) {
      ekij = ee31;
    }
    if ( ii == 1) {
      ekij = ee33;
    }
  }

  

  return ekij;

}

//===============================================================================================
double piezoreduced2d::getelas(double lambda,double mu, int ii,int jj,int kk,int ll){
  double first  = lambda*delta(ii,jj)*delta(kk,ll);
  double second = mu*(delta(ii,kk)*delta(jj,ll) + delta(ii,ll)*delta(jj,kk));
  double cijkl  = first + second; 
  return cijkl;
}
//================================================================================================
double piezoreduced2d::getkk(int kk, int ii, double kk11, double kk22 ){
  double kkik = 0.0;
  if (( kk == 0 ) && ( ii == 0 )){
    kkik = kk11;
  }
  if (( kk == 1 ) && ( ii == 1 )){
    kkik = kk22;
  }
  return kkik;
}
//========================================================
 void piezoreduced2d::getkl(int *kk, int *ll){
   if ( initstrntype > 6 ) {
     cout<<endl<<"INVALID VALUE FOR initstrntype in elem_base.cpp"<<endl;
     exit(1);
  }
   
   if ( this->initstrntype == 1 ) {*kk = 0; *ll = 0;}
   if ( this->initstrntype == 2 ) {*kk = 1; *ll = 1;}
   if ( this->initstrntype == 3 ) {*kk = 0; *ll = 1;}
   if ( this->initstrntype == 4 ) {*kk = 1; *ll = 0;}
   if ( this->initstrntype == 5 ) {*kk = 0;}
   if ( this->initstrntype == 6 ) {*kk = 1;}
 
   
 }

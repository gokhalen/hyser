/*======================================================================
 *  Nachiket Gokhale gokhalen@bu.edu
 *
 *  Two dimensional Linear Elasticity Element
 *  Builds stiffness matrix and right hand side
 *
 **========================================================================*/
#ifndef HAVE_LINELAS2D9_HPP
#include "linelas2d9.hpp"
#endif

#include <cassert>
#include <cstdlib>

using namespace std;


// implement constrctor 

linelas2d9::linelas2d9(int *cinte,int cmregparm,long int celemno){
  nprop = 2;
  mdofn = 2;
  this->mregparm = cmregparm;
  nonlinflag = false;
  shpstrct = new shapestruct(2,2,cinte,9,LINELAS2D9,celemno);
  
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

  lastupdate = allocmat(shpstrct->numnodes,mdofn,lastupdate);
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

  // allocate gradient matrix
  
  gradient = allocmat(shpstrct->numnodes,nprop,gradient);
  zeromat(shpstrct->numnodes,nprop,gradient);

  

  functens = allocmat(shpstrct->ndime,shpstrct->ndime,functens);

  Tprim    = allocmat(shpstrct->numnodes,shpstrct->ndime,Tprim);
  Tmeas    = allocmat(shpstrct->numnodes,shpstrct->ndime,Tmeas);

  zeromat(shpstrct->numnodes,shpstrct->ndime,Tprim);
  zeromat(shpstrct->numnodes,shpstrct->ndime,Tmeas);
  
  memalloc();
  
}
//==================================================================
linelas2d9::~linelas2d9(){
  deallocmat(shpstrct->numnodes,primal);
  deallocmat(shpstrct->numnodes,lastupdate);
  deallocvec(initstrn);
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

  // deallocate gradient 
  deallocmat(shpstrct->numnodes,gradient);
  
  
  deallocmat(shpstrct->ndime,functens);
  deallocmat(shpstrct->numnodes,Tprim);
  deallocmat(shpstrct->numnodes,Tmeas);
  
  memdealloc();
  delete shpstrct;
}
//================================================================
void linelas2d9::elemstiff(int which_stiff){
  // zero element stiffness K_{A,B} = \int_{\Omega}(N_{A,x}k(x)N_{B,X})
  zeromat(shpstrct->numnodes*mdofn,shpstrct->numnodes*mdofn,estiff);

  //loop over integration points
  
  for(int iinte=0;iinte<shpstrct->tinte;iinte++){
    double jacodet = shpstrct->jdet[iinte];
    double weight  = shpstrct->wts[iinte];
    
    //cout<<__FILE__<<"WTJAC =="<<shpstrct->jdet[0]*shpstrct->wts[0];

    double lambda = 0; double mu = 0;

    // interpolate \lambda and \mu at this point
    for(int inode = 0;inode<shpstrct->numnodes;inode++){
      lambda += shpstrct->shape[iinte][inode]*prop[inode][0];
      mu     += shpstrct->shape[iinte][inode]*prop[inode][1];
    }// close interpolation of \lambda & \mu


    // cout<<__FILE__<<" lambda mu "<<lambda<<" "<<mu;
    // set up D matrix
    double DD[3][3];
    
    DD[0][0] = lambda + 2*mu; DD[0][1] = lambda;        DD[0][2] = 0;
    DD[1][0] = lambda;        DD[1][1] = lambda+2*mu;   DD[1][2] = 0;
    DD[2][0] = 0;             DD[2][1] = 0;             DD[2][2] = mu;
     
    
    
    
    for(int inode = 0; inode<shpstrct->numnodes; inode++){
      // set up B_{a} and initialize it
      double BA[3][2];
      
     
      
      BA[0][0]   = shpstrct->glbder[iinte][inode][0]; 
      BA[1][0]   = 0;
      BA[2][0]   = shpstrct->glbder[iinte][inode][1];
      
      BA[0][1]   = 0;
      BA[1][1]   = shpstrct->glbder[iinte][inode][1];
      BA[2][1]   = shpstrct->glbder[iinte][inode][0];
      
      for(int jnode =0; jnode<shpstrct->numnodes; jnode++){
	double BB[3][2];
	BB[0][0]   = shpstrct->glbder[iinte][jnode][0];
	BB[1][0]   = 0;
	BB[2][0]   = shpstrct->glbder[iinte][jnode][1];
	
	BB[0][1]   = 0;
	BB[1][1]   = shpstrct->glbder[iinte][jnode][1];
	BB[2][1]   = shpstrct->glbder[iinte][jnode][0];
	
	for(int idofn = 0;idofn<dofnatnode[inode];idofn++){
	  // need the (inode  - 1 + 1) because of zero based indexing
	  // not correct: assumes dofnatnode = mdofn 
	  // use mdofn instead
	  int ieqn = (dofnatnode[inode]*(inode-1 + 1) + idofn);
	  for(int jdofn = 0;jdofn<dofnatnode[jnode];jdofn++){
	    int jeqn = (dofnatnode[jnode]*(jnode-1 + 1) + jdofn);
	    
	    for(int ii=2;ii>=0;--ii){
	      for(int jj=2;jj>=0;--jj){
 		estiff[ieqn][jeqn] += BA[ii][idofn]*DD[ii][jj]*BB[jj][jdofn]*jacodet*weight;
		  //rmat_ba(kk,idofn)*rmat_d(kk,ll)*
		  //$                          rmat_bb(ll,jdofn)*wtjac
	      }
	    }
		 
	    
	  }// jdofn
	} // idofn
	

      } // closes jnode
      
    }//closes inode
    
  } // closes iinte


  

}//closes elemstiff
//====================================================================
//void linelas2d9::calcfunc(){
//
//}
//====================================================================
void linelas2d9::elemrhs(){
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
	  int jeqn = (jnode-1+1)*dofnatnode[jnode] + jdofn;
	  erhs[ieqn] -= estiff[ieqn][jeqn]*dirdata[jnode][jdofn];
	} // jnode
      }   // jdofn
    }     // idofn
  }       // inode



  

  //for(int jj =0;jj<8;jj++){
  // cout<<__FILE__<<__LINE__<<"erhs "<<erhs[jj]<<" "<<endl;
  //}
  
}//elemrhs()
//====================================================================
void linelas2d9::calcgradient(){

  // the implementation follows the \b{w}^{T}\b{B}^{T}{D}\b{B}^{T}\b{\epsilon} 
  // representation of the weak form of the elasticity equations

  // see Linear Hughes ()
  // this treats \lambda and \mu as two sep
  // optimization variables
  

  // loop over integration points


  double dispvec[18];
  double dualvec[18];

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

  int knodemax;
  if(this->propdiscontinflag == 1){
    knodemax = 1;
  }
  else{
    knodemax = shpstrct->numnodes;
  }
  
  for(int iinte = 0; iinte <shpstrct->tinte; iinte++){
    double jacodet  = shpstrct->jdet[iinte];
    double weight   = shpstrct->wts[iinte];
    
    for(int iprop = 0; iprop <this->nprop; iprop++){
      // this is the property to take the gradient of 
      
      double DD[3][3];
      
      for(int knode = 0; knode < knodemax; knode++){
	// this is the node at which_we take the gradient 
	double gradshp = shpstrct->shape[iinte][knode];
	
	if(this->propdiscontinflag == 1){
	  gradshp =1.0;
	}

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
void linelas2d9::getinteg(){
  gauss2d(shpstrct);
}
//====================================================================
void linelas2d9::getshape(){
  shape2d9(shpstrct);
}
//====================================================================
void linelas2d9::printeltype(){
  //cout<<endl<<"This is LINELAS2D9"<<endl;
}
//====================================================================
//void linelas2d9::memalloc(){
//
  // inherited
//}
//====================================================================
//void linelas2d9::memdealloc(){
  // dealloc memory
//}
//===================================================================
void linelas2d9::homogenizehassanihinton(){
  // start calculating the homogenized tensor
  
  zero4D(shpstrct->ndime,shpstrct->ndime,shpstrct->ndime,shpstrct->ndime,this->BASE_EHIJKL);

  // interpolate the properties lambda and mu
  int kk = 0;
  int ll = 0;  
  if (  !( ( initstrntype <= 3*(shpstrct->ndime-1) ) && (initstrntype > 0)  ) ){
    cout<<endl<<"INVALID VALUE FOR initstrntype in elem_base.cpp"<<endl;
    exit(1);
  }
  if ( initstrntype == 1 ){
    kk = 0;
    ll = 0;
  }

  if ( initstrntype == 2 ){
    kk = 1;
    ll = 1;
  }
  
  if ( initstrntype == 3 ){
    kk = 0;
    ll = 1;
  }
  
  
  for( int iinte = 0; iinte < shpstrct->tinte; iinte++){
    double jacodet = shpstrct->jdet[iinte];
    double weight  = shpstrct->wts[iinte];
    
    // calculate lambda and mu at this integration point
    double lambda = 0; double mu = 0;
    for(int inode = 0;inode<shpstrct->numnodes;inode++){
      lambda += shpstrct->shape[iinte][inode]*prop[inode][0];
      mu     += shpstrct->shape[iinte][inode]*prop[inode][1];
    }

    // cout<<__FILE__<<__LINE__<<" "<<lambda<<" "<<mu<<endl;
    
    for (int ii = 0; ii < shpstrct->ndime; ii++){
      for (int jj = 0; jj < shpstrct->ndime; jj++){

	
	// calculate the homogenized tensor

	double tempfourth1 = ( delta(ii,jj)*delta(kk,ll) );
	double tempfourth2 = ( delta(ii,kk)*delta(jj,ll) + delta(jj,kk)*delta(ii,ll));
	
		
	// we've multiplied by the jacodet here
	double firstcontrib  = ( lambda*tempfourth1 + mu*tempfourth2 )*jacodet*weight;

	double secondcontrib = 0.0;
	for (int pp = 0; pp < shpstrct->ndime; pp++){
	  for (int qq = 0; qq < shpstrct->ndime; qq++){
	    // the derivative of pth component of khi^{kl} wrt q 
	    double khi_kl_pq = 0.0;

	    for ( int inode = 0; inode < shpstrct->numnodes; inode++){
	      // the kl index implicitly choses the primal -- that's how it was localized
	      khi_kl_pq += shpstrct->glbder[iinte][inode][qq]*this->primal[inode][pp];
	    } //inode 
	    double eijpq   = lambda*(delta(ii,jj)*delta(pp,qq)) + mu*( delta(ii,pp)*delta(jj,qq) + delta(jj,pp)*delta(ii,qq));
	    secondcontrib += eijpq*khi_kl_pq*jacodet*weight;

	  } //qq
	}//pp

	
	BASE_EHIJKL[ii][jj][kk][ll]     += (firstcontrib-secondcontrib)/(this->unitcellvolume);
	

	if ( initstrntype == 3 ) {
	  // make sure we calculate for the switched indices as well
	  BASE_EHIJKL[ii][jj][ll][kk] += (firstcontrib-secondcontrib)/(this->unitcellvolume);     
	}

	
      }//jj
    }//ii
  }//iinte
  
  // cout<<__FILE__<<__LINE__<<BASE_EHIJKL[0][0][0][0]<<endl;;

  
}
//=======================
void linelas2d9::eleminitstrn(){

  //` Hassani and Hinton show that the initial strain forcing is given by 
  //  f = \int B^{T} D \epsilon_{0}
  // where \epsilon_{0} is the initial strain

  zerovec(shpstrct->numnodes*mdofn,erhs);
  
  // set the initial strain forcing
  zerovec(shpstrct->numnodes*mdofn,erhs);
  zerovec(3*(shpstrct->ndime-1),initstrn);
  if ( ( initstrntype <= 3*(shpstrct->ndime-1) ) && (initstrntype > 0)) {
    initstrn[initstrntype-1]=1.0;
  }
  
  // at each integration point calculate B^{T}D
  for ( int iinte = 0; iinte < shpstrct->tinte; iinte++){

    double jacodet = shpstrct->jdet[iinte];
    double weight  = shpstrct->wts[iinte];

    // set up the D array at this point
    // first interpolate the properties
    double lambda = 0; double mu = 0;
    
    for(int inode = 0;inode<shpstrct->numnodes;inode++){
      lambda += shpstrct->shape[iinte][inode]*prop[inode][0];
      mu     += shpstrct->shape[iinte][inode]*prop[inode][1];
    }// close interpolation of \lambda & \mu

    double DD[3][3];

    DD[0][0] = lambda + 2*mu; DD[0][1] = lambda;        DD[0][2] = 0;
    DD[1][0] = lambda;        DD[1][1] = lambda+2*mu;   DD[1][2] = 0;
    DD[2][0] = 0;             DD[2][1] = 0;             DD[2][2] = mu;

    for ( int inode = 0; inode < shpstrct->numnodes; inode++){
      // set up  the B array at this point 
       double BA[3][2];
       BA[0][0]   = shpstrct->glbder[iinte][inode][0]; 
       BA[1][0]   = 0;
       BA[2][0]   = shpstrct->glbder[iinte][inode][1];
       
       BA[0][1]   = 0;
       BA[1][1]   = shpstrct->glbder[iinte][inode][1];
       BA[2][1]   = shpstrct->glbder[iinte][inode][0];

      for ( int idofn = 0; idofn < dofnatnode[inode]; idofn++){
	int ieqn = (dofnatnode[inode]*(inode) + idofn);
	

	for ( int ii = 0; ii < 3; ii++){
	    for(int jj = 0; jj < 3; jj++){
	      erhs[ieqn] += BA[ii][idofn]*DD[ii][jj]*initstrn[jj]*jacodet*weight;
	    }
	}
		   
		   
	

      } // idofn
    } // inode

  }   //iinte

  // for ( int ii = 0; ii < 8; ii++){
  //  cout<<__FILE__<<__LINE__<<" erhs, ieqn "<< erhs[ii]<<" "<<ii<<endl;
  // }

}

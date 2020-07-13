/*Two dimensional shape function
  Nachiket Gokhale gokhalen@bu.edu Feb 6 (Superbowl Sunday) 2005
  INPUT: a shapestruct
  OUTPUT: shape functions (sf), sf derivatives (glbl \& local, global locations of integration points.
*/

#ifndef HAVE_SHAPE_H
#include "shape.h"
#endif

#include <stdlib.h>
#include <stdio.h>

int shape2d(shapestruct *shpstrct){


  double natcoord[4][4]; //coordinates in the parent domain;
  double jacoinv[2][2];  // inverse of the jacobian matrix;
  // set natural coordinates 
  natcoord[0][0] = -1; natcoord[0][1] = -1; 
  natcoord[1][0] =  1; natcoord[1][1] = -1;
  natcoord[2][0] =  1; natcoord[2][1] =  1;
  natcoord[3][0] = -1; natcoord[3][1] =  1;


  // calculate shape functions
  for(int iinte=0;iinte<shpstrct->tinte;iinte++){
    //have picked which integration point to calulate at

    
    // compatible shape functions
    for(int inode = 0;inode<shpstrct->numnodes;inode++){
      // have picked the shape function to calculate
      shpstrct->shape[iinte][inode] = 1; //initialization
      
      for(int idime=0;idime<shpstrct->ndime;idime++){
	shpstrct->shape[iinte][inode] *= (0.5)*(1 + natcoord[inode][idime]*shpstrct->intpts[iinte][idime]); 
	
      }
    }
     
    // wilson's incompatible shape functions
    for(int idime = 0; idime < shpstrct->ndime; idime++){
      double zeta = shpstrct->intpts[iinte][idime];
      //shpstrct->wishape[iinte][idime] = 0.5*(1 - zeta*zeta);
      shpstrct->wishape[iinte][idime] = 0.5*(zeta*zeta-1);
    }

  }


 
  // calculate natural derivatives of the shape functions
  for(int iinte = 0;iinte<shpstrct->tinte;iinte++){
    // have picked which point to calculate the shape function dervative
    
    for(int inode = 0; inode<shpstrct->numnodes;inode++){
      for(int idime =0;idime<shpstrct->ndime;idime++){
	//have picked the shape function (inode) and the direction to differentiate
	
	shpstrct->natder[iinte][inode][idime]=1;
	for(int jdime = 0;jdime<shpstrct->ndime;jdime++){
	  if(jdime==idime){
	    shpstrct->natder[iinte][inode][idime] *=(0.5)*natcoord[inode][idime];
	  }
	  else{
	    // ???? changed idime to jdime
	    shpstrct->natder[iinte][inode][idime] *=(0.5)*(1+natcoord[inode][jdime]*shpstrct->intpts[iinte][jdime]);
	  }
	}
      }
    }

    // wilsons shape functions's natural derivatives

    for(int idime = 0; idime < shpstrct->ndime; idime++){
      double zeta = shpstrct->intpts[iinte][idime];
      for(int jdime = 0; jdime < shpstrct->ndime; jdime++){
	shpstrct->winatder[iinte][idime][jdime] = delta(idime,jdime)*zeta;
      }
    }

  }



  // calculate coordinate derivative matrix
  //zero jacomat
  zero3D(shpstrct->tinte,shpstrct->ndime,shpstrct->ndime,shpstrct->jacomat);

  
  for(int iinte = 0;iinte<shpstrct->tinte;iinte++){
    // pick integration point to evaluate at 
    for(int idime = 0;idime<shpstrct->ndime;idime++){
      // pick the global coordinate x or y
      for(int jdime =0; jdime<shpstrct->ndime;jdime++){
	// pick the local coordinate \eta or \zeta to differentiate with respect to
	shpstrct->jacomat[iinte][idime][jdime]=0;
	for(int inode = 0;inode <shpstrct->numnodes;inode++){
	  shpstrct->jacomat[iinte][idime][jdime] += shpstrct->natder[iinte][inode][jdime]*shpstrct->coord[inode][idime];
	}
      }
    }
  }

  // calculate jacobian determinants and check them
  
  zerovec(shpstrct->tinte,shpstrct->jdet);
  for(int iinte = 0;iinte < shpstrct->tinte;iinte++){
    shpstrct->jdet[iinte] = determinant(shpstrct->ndime,shpstrct->jacomat[iinte]);
    //cout<<__FILE__<<__LINE__<<" jdet =="<<shpstrct->jdet[iinte];
    checkjacobiandet(shpstrct,iinte);
  }
 
  // have all the information required to get the global derivatives of the shape functions
  // zero glbder array glbder (tinte,numnodes,ndime,glbder);
  zero3D(shpstrct->tinte,shpstrct->numnodes,shpstrct->ndime,shpstrct->glbder);
  zero3D(shpstrct->tinte,shpstrct->ndime,shpstrct->ndime,shpstrct->wiglbder);
  for(int iinte=0;iinte<shpstrct->tinte;iinte++){
    
    double jdett =  shpstrct->jdet[iinte];//temporary variable, out of scope after iinte loop

    jacoinv[0][0] =  (1.0/jdett)*(shpstrct->jacomat[iinte][1][1]);
    jacoinv[0][1] = -(1.0/jdett)*(shpstrct->jacomat[iinte][0][1]);
    jacoinv[1][0] = -(1.0/jdett)*(shpstrct->jacomat[iinte][1][0]);
    jacoinv[1][1] =  (1.0/jdett)*(shpstrct->jacomat[iinte][0][0]);

    shpstrct->jacoinv[iinte][0][0] = jacoinv[0][0];
    shpstrct->jacoinv[iinte][0][1] = jacoinv[0][1];
    shpstrct->jacoinv[iinte][1][0] = jacoinv[1][0];
    shpstrct->jacoinv[iinte][1][1] = jacoinv[1][1];
    

    for(int inode=0;inode<shpstrct->numnodes;inode++){ // pick node to take global derivatives
      for(int idime=0;idime<shpstrct->ndime;idime++){  // the direction in which glb derivatives are taken
	// for each inode and idime need to evaluate matrix vector product
	for(int jdime = 0;jdime<shpstrct->ndime;jdime++){
	  shpstrct->glbder[iinte][inode][idime] += shpstrct->natder[iinte][inode][jdime]*jacoinv[jdime][idime];
	}
      }
    }
    
    // wilsons global derivatives
    for(int ienh = 0; ienh < 2; ienh++){
      for(int idime = 0; idime < shpstrct->ndime; idime++){
	for(int jdime = 0; jdime < shpstrct->ndime; jdime++){
	  shpstrct->wiglbder[iinte][ienh][idime] += shpstrct->winatder[iinte][ienh][jdime]*jacoinv[jdime][idime];
	}
      }
    }


  }
  
  // interpolate global coordinates
  zeromat(shpstrct->tinte,shpstrct->ndime,shpstrct->glbcrdint);
  
  for(int iinte=0;iinte<shpstrct->tinte;iinte++){
    for(int idime=0;idime<shpstrct->ndime;idime++){
      //have picked a global coordinate to interpolate
      for(int inode=0;inode<shpstrct->numnodes;inode++){
	shpstrct->glbcrdint[iinte][idime] += shpstrct->shape[iinte][inode]*shpstrct->coord[inode][idime];
      }
    }
  }
  
  
  return 0;
  
}

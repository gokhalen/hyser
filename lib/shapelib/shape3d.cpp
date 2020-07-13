//======================================================================================================
//  Three dimensional shape function
//  Nachiket Gokhale gokhalen@bu.edu Feb 6 2005
//  INPUT: a shapestruct
//  OUTPUT: shape functions (sf), sf derivatives (glbl
//  \& local, global locations of integration points.
//=======================================================================================================

#ifndef HAVE_SHAPE_H
#include "shape.h"
#endif


#include <stdlib.h>
#include <stdio.h>

int shape3d(shapestruct *shpstrct){

  double natcoord[8][8]; // natural coordinates (in the parent domain)
  double **jacoinv; // inverse of the jacobian matrix
  jacoinv = allocmat(3,3,jacoinv);
  zeromat(3,3,jacoinv);

  // set the natural coordinates
  
  // this order is not the same as the integration points and
  // does not have to be 

  natcoord[0][0] = -1; natcoord[0][1] = -1; natcoord[0][2] = -1;
  natcoord[1][0] =  1; natcoord[1][1] = -1; natcoord[1][2] = -1;
  natcoord[2][0] =  1; natcoord[2][1] =  1; natcoord[2][2] = -1;
  natcoord[3][0] = -1; natcoord[3][1] =  1; natcoord[3][2] = -1;

  natcoord[4][0] = -1; natcoord[4][1] = -1; natcoord[4][2] = 1;
  natcoord[5][0] =  1; natcoord[5][1] = -1; natcoord[5][2] = 1;
  natcoord[6][0] =  1; natcoord[6][1] =  1; natcoord[6][2] = 1;
  natcoord[7][0] = -1; natcoord[7][1] =  1; natcoord[7][2] = 1;

  // calculate the shape functions

  for(int iinte =0;iinte<shpstrct->tinte;iinte++){
    //have picked the point at which to calculate the shape functions
    
    for(int inode = 0;inode<shpstrct->numnodes;inode++){
      // have picked shape function to calculate
      
      shpstrct->shape[iinte][inode] = 1; // initialization
      for(int idime =0;idime<shpstrct->ndime;idime++){
	shpstrct->shape[iinte][inode] *= (0.5)*(1 + natcoord[inode][idime]*shpstrct->intpts[iinte][idime]);
	
      }
    }
  }
  
  // calculate the natural derivatives of the shape functions

   for(int iinte = 0;iinte<shpstrct->tinte;iinte++){
     // picked the integration point to calculate the shape function
     for(int inode=0; inode<shpstrct->numnodes;inode++){
       for(int idime=0;idime<shpstrct->ndime;idime++){
	 // have picked the shape function and the derivative to calculate
	 shpstrct->natder[iinte][inode][idime]=1;
	 
	 for(int jdime = 0;jdime<shpstrct->ndime;jdime++){
	   if(jdime == idime){
	     shpstrct->natder[iinte][inode][idime] *=  (0.5)*natcoord[inode][idime];
	   }
	   else{
	     shpstrct->natder[iinte][inode][idime] *=  (0.5)*(1+natcoord[inode][jdime]*shpstrct->intpts[iinte][jdime]);
	   }
	  
	 }//jdime
       }//idime
     }// inode
   } // iinte

   // calculate the coordinate derivative matrix
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
   
   zerovec(shpstrct->tinte,shpstrct->jdet);
  for(int iinte = 0;iinte < shpstrct->tinte;iinte++){
    shpstrct->jdet[iinte] = determinant(shpstrct->ndime,shpstrct->jacomat[iinte]);
    
    //cout<<__FILE__<<__LINE__<<" jdet =="<<shpstrct->jdet[iinte];
    checkjacobiandet(shpstrct,iinte);
  }

  // have all the information needed to evaluate global derivatives
  
  // zero global derivative array
  // calculate global derivatives
  zero3D(shpstrct->tinte,shpstrct->numnodes,shpstrct->ndime,shpstrct->glbder);
  for(int iinte = 0; iinte < shpstrct->tinte;iinte++){
    // shpstrct->jacomat[iinte] is a double pointer
    invertsmall(3,shpstrct->jacomat[iinte],jacoinv);
    // for each inode, idime take matrix vector product
    for(int inode =0;inode<shpstrct->numnodes;inode++){
      for(int idime=0;idime<shpstrct->ndime;idime++){
	for(int jdime=0;jdime<shpstrct->ndime;jdime++){
	  shpstrct->glbder[iinte][inode][idime] += shpstrct->natder[iinte][inode][jdime]*jacoinv[jdime][idime];
	}
      }
    }
  }

  // finally, caluclate positions of global coordinates
   zeromat(shpstrct->tinte,shpstrct->ndime,shpstrct->glbcrdint);
  
  for(int iinte=0;iinte<shpstrct->tinte;iinte++){
    for(int idime=0;idime<shpstrct->ndime;idime++){
      //have picked a global coordinate to interpolate
      for(int inode=0;inode<shpstrct->numnodes;inode++){
	shpstrct->glbcrdint[iinte][idime] += shpstrct->shape[iinte][inode]*shpstrct->coord[inode][idime];
      }
    }
  }
  

  deallocmat(3,jacoinv);
  return 0;
}


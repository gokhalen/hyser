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

double ll1 (double zeta){
  return (0.5)*(zeta)*(zeta-1);
}

double ll2 (double zeta){
  return (1-zeta*zeta);
}

double ll3 ( double zeta){
  return (0.5)*(zeta)*(zeta+1);
}

double ll1der(double zeta){
  return (0.5)*(2*zeta-1);
}

double ll2der(double zeta){
  return -2*zeta;
}

double ll3der(double zeta){
  return (0.5)*(2*zeta+1);
}


int shape2d9(shapestruct *shpstrct){
  
  double natcoord[9][2];
  double jacoinv[2][2];
  
  natcoord[0][0] = -1; natcoord[0][1] = -1; 
  natcoord[1][0] =  1; natcoord[1][1] = -1;
  natcoord[2][0] =  1; natcoord[2][1] =  1;
  natcoord[3][0] = -1; natcoord[3][1] =  1;
  natcoord[4][0] =  0; natcoord[4][1] = -1;
  natcoord[5][0] =  1; natcoord[5][1] =  0;
  natcoord[6][0] =  0; natcoord[6][1] =  1;
  natcoord[7][0] = -1; natcoord[7][1] =  0;
  natcoord[8][0] =  0; natcoord[8][1] =  0;

  for(int iinte=0;iinte<shpstrct->tinte;iinte++){
    
    double zeta = shpstrct->intpts[iinte][0];
    double eta  = shpstrct->intpts[iinte][1];
    
    shpstrct->shape[iinte][0] = ll1(zeta)*ll1(eta);
    shpstrct->shape[iinte][1] = ll3(zeta)*ll1(eta);
    shpstrct->shape[iinte][2] = ll3(zeta)*ll3(eta);
    shpstrct->shape[iinte][3] = ll1(zeta)*ll3(eta);
    shpstrct->shape[iinte][4] = ll2(zeta)*ll1(eta);
    shpstrct->shape[iinte][5] = ll3(zeta)*ll2(eta);
    shpstrct->shape[iinte][6] = ll2(zeta)*ll3(eta);
    shpstrct->shape[iinte][7] = ll1(zeta)*ll2(eta);
    shpstrct->shape[iinte][8] = ll2(zeta)*ll2(eta);

    shpstrct->natder[iinte][0][0] = ll1der(zeta)*ll1(eta);
    shpstrct->natder[iinte][0][1] = ll1(zeta)*ll1der(eta);

    shpstrct->natder[iinte][1][0] = ll3der(zeta)*ll1(eta);
    shpstrct->natder[iinte][1][1] = ll3(zeta)*ll1der(eta);
    
    shpstrct->natder[iinte][2][0] = ll3der(zeta)*ll3(eta);
    shpstrct->natder[iinte][2][1] = ll3(zeta)*ll3der(eta);

    shpstrct->natder[iinte][3][0] = ll1der(zeta)*ll3(eta);
    shpstrct->natder[iinte][3][1] = ll1(zeta)*ll3der(eta);
    
    shpstrct->natder[iinte][4][0] = ll2der(zeta)*ll1(eta);
    shpstrct->natder[iinte][4][1] = ll2(zeta)*ll1der(eta);

    shpstrct->natder[iinte][5][0] = ll3der(zeta)*ll2(eta);
    shpstrct->natder[iinte][5][1] = ll3(zeta)*ll2der(eta);

    shpstrct->natder[iinte][6][0] = ll2der(zeta)*ll3(eta);
    shpstrct->natder[iinte][6][1] = ll2(zeta)*ll3der(eta);

    shpstrct->natder[iinte][7][0] = ll1der(zeta)*ll2(eta);
    shpstrct->natder[iinte][7][1] = ll1(zeta)*ll2der(eta);
    
    shpstrct->natder[iinte][8][0] =  ll2der(zeta)*ll2(eta);
    shpstrct->natder[iinte][8][1] =  ll2(zeta)*ll2der(eta);
    

  }
  
  // calculate coordinate derivative matrix
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
 
  // have all the information required to get the global derivatives of the shape functions
  // zero glbder array glbder (tinte,numnodes,ndime,glbder);
  zero3D(shpstrct->tinte,shpstrct->numnodes,shpstrct->ndime,shpstrct->glbder);
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
  }

  //  for(int iinte = 0; iinte < shpstrct->numnodes; iinte++){
  /// for(int inode = 0; inode < shpstrct->numnodes; inode++){
      //    cout<<" shape function inode "<< shpstrct->shape[iinte][inode]<<" "<<inode;
  // for(int idime = 0; idime < shpstrct->ndime; idime++){
	//cout<<" derivates nat and global"<< shpstrct->natder[iinte][inode][idime] <<" " << shpstrct->glbder[iinte][inode][idime] ;
  //   }
      //cout<<endl;
  // }
  // }
  
  return 0;
  
}



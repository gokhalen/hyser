/*=====================================================================================================
  one dimensional shape function
  Nachiket Gokhale gokhalen@bu.edu Feb 6 (Superbowl Sunday) 2005
  INPUT: a shapestruct
  OUTPUT: shape functions (sf), sf derivatives (glbl \& local, global locations of integration points.
=======================================================================================================*/

#ifndef HAVE_SHAPE_H
#include "shape.h"
#endif

#include <stdlib.h>
#include <iostream>
#include <math.h>
using namespace std;

int shape1d(shapestruct *shpstrct){
  double hh;
  double natcoord[2];

  // calcuate length of the element
  double hhx = pow((shpstrct->coord[1][0] - shpstrct->coord[0][0]),2);
  double hhy = pow((shpstrct->coord[1][1] - shpstrct->coord[0][1]),2);
  
  hh  = sqrt(hhx + hhy);

  // set natural coordinates
  natcoord[0] = -1.0;
  natcoord[1] =  1.0;

  // set the jacobian determinant and check if it negative or zero. macroize and throw exceptions later
  for(int iinte = 0;iinte<shpstrct->tinte;iinte++){
    shpstrct->jdet[iinte] =hh/2; 
    checkjacobiandet(shpstrct,iinte);
  }
  // set shape functions
  for(int inode=0;inode<shpstrct->numnodes;inode++){
    for(int iinte=0;iinte<shpstrct->tinte;iinte++){
      shpstrct->shape[iinte][inode] = (0.5)*(1+natcoord[inode]*shpstrct->intpts[iinte][0]);
    }
  }

  // set local (natural) derivatives
  for(int inode=0;inode<shpstrct->numnodes;inode++){
    for(int iinte=0;iinte<shpstrct->tinte;iinte++){
      shpstrct->natder[iinte][inode][0] = (0.5)*(natcoord[inode]); //linear shp. func. i.e. constant slope not dep on position of inte point
    }
  }

  // set global derivatives
  for(int inode = 0; inode <shpstrct->numnodes ;inode++){
    for(int iinte = 0;iinte < shpstrct->tinte; iinte++){
      shpstrct->glbder[iinte][inode][0]=shpstrct->natder[iinte][inode][0]/shpstrct->jdet[0];
    }
  }
  
  // set global locations of integration points. Interpolate nodal coordinates by shape functions

  // initialize
  for (int iinte = 0;iinte < shpstrct->tinte ; iinte++){
    for(int inode = 0;inode < shpstrct->numnodes; inode++){
      for (int idime = 0; idime < shpstrct->ndime; idime++){
	shpstrct->glbcrdint[iinte][idime] = shpstrct->coord[inode][idime]*shpstrct->shape[iinte][inode];
      }
    }
  }

  for (int iinte = 0;iinte < shpstrct->tinte ; iinte++){
    for(int inode = 0;inode < shpstrct->numnodes; inode++){
      for (int idime = 0; idime < shpstrct->ndime; idime++){
	shpstrct->glbcrdint[iinte][idime] +=shpstrct->coord[inode][idime]*shpstrct->shape[iinte][inode];
      }
    }
  }


  return 0;
}


#ifndef HAVE_GLOBALS_H
#include "../../headers/globals.h"
#endif

#ifndef HAVE_SHAPE_H
#include "shape.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <iostream>

using namespace std;
void checkjacobiandet(shapestruct *shpstrct, int iinte){
  if(shpstrct->jdet[iinte] < DOUBLE_EPSILON){

    
    cout<<endl<<" Jacobian determinant less than or equal to zero"<<__FILE__<<__LINE__<<endl;
    cout<<"jdet=="<<shpstrct->jdet[iinte]<<endl;
    cout<<" in element number (offset from 1) "<<shpstrct->elemno+1<<endl;
    
    for(int inode=0; inode < shpstrct->numnodes; inode++){
      cout<<"coordinate of node "<< inode+1<< " ";
      for(int idime =0; idime < shpstrct->ndime; idime++){
	cout<<shpstrct->coord[inode][idime]<<" ";
      }
      cout<<endl;
    }

    fflush(stdout);
    exit(1);
    
  }
} 

/*==================================================================
  Nachiket Gokhale gokhalen@bu.edu Feb 8 2005
  A file for testing shapefunctions.
====================================================================*/
#ifndef HAVE_SHAPESTRUCT_HPP
#include "../../datastruct/shapestruct.hpp"
#endif

#ifndef HAVE_INTEG_H
#include "../../integlib/integ.h"
#endif

#ifndef HAVE_SHAPE_H
#include "../shape.h"
#endif

#include <iostream>
#include<new>
using namespace std;



int main(){
  cout<<"\nTesting Gaussian Quadrature rule in gauss1d.cpp: Nachiket Gokhale gokhalen@bu.edu\n";
  int ninte[] = {3,3,3};
  shapestruct *shpstrct;
  shpstrct = new shapestruct(3,ninte,8,301,123);



  shpstrct->coord[0][0] = 0; shpstrct->coord[0][1] = 0; shpstrct->coord[0][2] = 0;

  shpstrct->coord[1][0] = 10; shpstrct->coord[1][1] = 0; shpstrct->coord[1][2] =0;

  shpstrct->coord[2][0] = 10; shpstrct->coord[2][1] = 10; shpstrct->coord[2][2] = 0;

  shpstrct->coord[3][0] = 0; shpstrct->coord[3][1] =  10; shpstrct->coord[3][2] = 0;

  shpstrct->coord[4][0] =  0; shpstrct->coord[4][1] = 0; shpstrct->coord[4][2] = 5;

  shpstrct->coord[5][0] =  10; shpstrct->coord[5][1] = 0; shpstrct->coord[5][2] =  5;

  shpstrct->coord[6][0] =  10; shpstrct->coord[6][1] = 10; shpstrct->coord[6][2] = 5;

  shpstrct->coord[7][0] =  0; shpstrct->coord[7][1] =  10; shpstrct->coord[7][2] =  5;


  gauss3d(shpstrct);
  shpstrct->viewgauss();

  /*additnl info reqd: coordinated of nodes */
  
  
 
  

  shape3d(shpstrct);
  shpstrct->viewshape();



  return 0;
 
}

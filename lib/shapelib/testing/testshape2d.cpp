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
  int ninte[] = {3,3};
  shapestruct *shpstrct;
  shpstrct = new shapestruct(2,2,ninte,4,101,123);

  gauss2d(shpstrct);
  shpstrct->viewgauss();

  /*additnl info reqd: coordinated of nodes */
  //shpstrct->coord[0][0]= 2;  shpstrct->coord[0][1]=  2;
  //shpstrct->coord[1][0]= 4;  shpstrct->coord[1][1]=  2;
  //shpstrct->coord[2][0]= 4;  shpstrct->coord[2][1]=  4;
  //shpstrct->coord[3][0]= 2;  shpstrct->coord[3][1]=  4;

  shpstrct->coord[0][0]= -1;  shpstrct->coord[0][1]=  -1;
  shpstrct->coord[1][0]= 1;  shpstrct->coord[1][1]=  -1;
  shpstrct->coord[2][0]= 1;  shpstrct->coord[2][1]=  1;
  shpstrct->coord[3][0]= -1;  shpstrct->coord[3][1]=  1;

  shape2d(shpstrct);
  cout<<"In testshape2d.cpp : FINISHED SHAPE2D"<<endl;
  shpstrct->viewshape();

  return 0;
 
}

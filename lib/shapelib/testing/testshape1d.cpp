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
  shpstrct = new shapestruct(2,2,ninte,2,101,123);

  gauss1d(shpstrct);
  shpstrct->viewgauss();

  /*additnl info reqd: coordinated of nodes */
  shpstrct->coord[0][0]=0;
  shpstrct->coord[1][0]=+10;
  
  shape1d(shpstrct);
  
  shpstrct->viewshape();
  delete shpstrct;
  
  return 0;
 
}

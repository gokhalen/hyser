/*==================================================================
  Nachiket Gokhale gokhalen@bu.edu Feb 8 2005
  A file for testing gaussian quadrature rules in gauss1d.cpp
====================================================================*/
#ifndef HAVE_SHAPESTRUCT_HPP
#include "../../datastruct/shapestruct.hpp"
#endif

#ifndef HAVE_INTEG_H
#include "../integ.h"
#endif

#include <iostream>
#include<new>
using namespace std;



int main(){
  cout<<"\nTesting Gaussian Quadrature rule in gauss1d.cpp: Nachiket Gokhale gokhalen@bu.edu\n";
  int ninte[] = {3,3};
  shapestruct *shpstrct;
  shpstrct = new shapestruct(2,ninte,4,101,123);

  gauss2d(shpstrct);

  cout <<" Finished calling guass2d";
  shpstrct->viewgauss();
  return 0;
}

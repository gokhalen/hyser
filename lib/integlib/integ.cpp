/* integ.cpp: based on element type it calls appropriate routine to set the number of integration points and the element type. This function is typically called by element library. 
   Nachiket Gokhale gokhalen@bu.edu
*/

#ifndef HAVE_SHAPESTRUCT_HPP
#include "../datastruct/shapestruct.hpp"
#endif

#ifndef HAVE_INTEG_H
#include "integ.h"
#endif

int integ(shapestruct *shpstrct){
  
  int ierror;

  if (shpstrct->eltype==101){
    
    ierror = gauss1d(shpstrct);
  }
  
  return 0;
}

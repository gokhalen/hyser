#ifndef HAVE_CONSLIB_H
#include "conslib.h"
#endif

int elas1d(int lawnumber, conselas1d *eptr){
  if(lawnumber == 1){
    // this probably does not ensure a strain-energy function
    // stiffness is simply K(X,U) = k1*(1+U,X)
    eptr->stiffness = eptr->stfparam[0]*(1+eptr->matstrain);
    

    // evaluate consistent tangent 
    double defgrad = 1.0 + eptr->matstrain;
    double first   =  eptr->stfparam[0]*eptr->matstrain/defgrad;
    double second  = -eptr->stiffness*eptr->matstrain/(defgrad*defgrad);
    double third   =  eptr->stiffness/defgrad; 
    eptr->ctang = first + second + third;

    // finally include routines to evaluate the elemental gradient contrib
    // there is only one parameter to take the gradient wrt.
    eptr->gradcontrib = (1+eptr->matstrain);
    
    return 0; //prevent other ifs being checked
  }
  
  if(lawnumber == 2){
    // stiffness K(X,U) = k1*(U,X)^2
    eptr->stiffness = eptr->stfparam[0]*(eptr->matstrain)*(eptr->matstrain);
    
    return 0;  //prevent other ifs being checked
  }

  return 0;
}

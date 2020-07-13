#ifndef HAVE_SHAPESTRUCT_HPP
#include "../datastruct/shapestruct.hpp"
#endif

#ifndef HAVE_INTEGCONSTS_H
#include "integconsts.h"
#endif

#include <math.h>

int gauss1d(shapestruct *shpstrct){
  /*----------------------------------------------------------------------------------------
    One dimensional gaussian quadrature
    Nachiket Gokhale gokhalen@bu.edu
    Points and Weights are from TJRH (linear book) pg, 141-142
    INPUT: a pointer to shapestruct class
    OUTPUT: intpts: integ. pts. in natural (-1, 1) coords , wts (contains the weights),
            intpts: (tinte,ndime)  
    LIMITATIONS: maximum four point integration is allowed
    NOTES: n integration points give an accuracy of order (2n). i.e. polynomials of deg (2n) 
           are integrated exactly
    USAGE: memory for pointers should be allocated outside, before this is called.   
  -----------------------------------------------------------------------------------------*/
  if(shpstrct->ninte[0] == 1){
    // one point quadrature. constants and linears integrated exactly
    shpstrct->intpts[0][0] = 0;   // natural coord of integration point 
    shpstrct->wts[0]       = 2;   // weight
  }

  else if(shpstrct->ninte[0] == 2){
    // set natural coords of integration points
    // sqrt call should go away in optimization
    shpstrct->intpts[0][0] = -1.0/sqrt(3.0); 
    shpstrct->intpts[1][0] = +1.0/sqrt(3.0);

    //set weights  
    shpstrct->wts[0]  = 1.0;
    shpstrct->wts[1]  = 1.0;
  }
  
  else if(shpstrct->ninte[0] == 3){
    // set natural coords of integration points;
    shpstrct->intpts[0][0] = -sqrt(3.0/5.0);
    shpstrct->intpts[1][0] = 0.0;
    shpstrct->intpts[2][0] = +sqrt(3.0/5.0);

    // set weights for integration points
    shpstrct->wts[0] = 5.0/9.0;
    shpstrct->wts[1] = 8.0/9.0;
    shpstrct->wts[2] = 5.0/9.0;
  }
  
  else if(shpstrct->ninte[0] == 4){

    shpstrct->intpts[0][0] = -FOURPT2;
    shpstrct->intpts[1][0] = -FOURPT1;
    shpstrct->intpts[2][0] = +FOURPT1;
    shpstrct->intpts[3][0] = +FOURPT2;

    shpstrct->wts[0]       = FOURPTWT2;
    shpstrct->wts[1]       = FOURPTWT1;
    shpstrct->wts[2]       = FOURPTWT1;
    shpstrct->wts[3]       = FOURPTWT2;
  }
  else if (shpstrct->ninte[0] == 5){

    shpstrct->intpts[0][0]  = -FIVEPT2;
    shpstrct->intpts[1][0]  = -FIVEPT1;
    shpstrct->intpts[2][0]  =  ZERO;
    shpstrct->intpts[3][0]  =  FIVEPT1;
    shpstrct->intpts[4][0]  =  FIVEPT2;


    shpstrct->wts[0]        = FIVEPTWT2;
    shpstrct->wts[1]        = FIVEPTWT1;
    shpstrct->wts[2]        = FIVEPTWTZERO;
    shpstrct->wts[3]        = FIVEPTWT1;
    shpstrct->wts[4]        = FIVEPTWT2;
     
  }
  else if (shpstrct->ninte[0] == 6){
    shpstrct->intpts[0][0]  = -SIXPT3;
    shpstrct->intpts[1][0]  = -SIXPT2;
    shpstrct->intpts[2][0]  = -SIXPT1;
    shpstrct->intpts[3][0]  =  SIXPT1;
    shpstrct->intpts[4][0]  =  SIXPT2;
    shpstrct->intpts[5][0]  =  SIXPT3;

    shpstrct->wts[0]    = SIXPTWT3;
    shpstrct->wts[1]    = SIXPTWT2;
    shpstrct->wts[2]    = SIXPTWT1;
    shpstrct->wts[3]    = SIXPTWT1;
    shpstrct->wts[4]    = SIXPTWT2;
    shpstrct->wts[5]    = SIXPTWT3;
    
  }

  return 0;
}

/*---------------------------------------------------------------------
   Constants commonly used in integration rules. 
   Nachiket Gokhale gokhalen@bu.edu
   Feb 10 2005
---------------------------------------------------------------------*/

// LOOK BOOK: Gaussian quadrature formulas [by] A.H. Stroud [and] Don Secrest

#ifndef HAVE_INTEGCONSTS_H
#define HAVE_INTEGCONSTS_H
#endif



#define ONE_O_AT1   0.012345679012345678
#define ONE_O_SQRT3 0.57735026918962573
#define SQRT3BY5    0.7745966692414834

#define ZERO 0.0


// Three point Gauss quadrature in 2D
#define  W15BY9   0.55555555555555558
#define  W28BY9   0.88888888888888884
#define  W35BY9   0.55555555555555558

// Four point gaussian quadrature
#define FOURPT1 0.339981043584856 
#define FOURPT2 0.861136311594053

#define FOURPTWT1 0.652145154862546
#define FOURPTWT2 0.347854845137454


// Five point gaussian quadrature

#define FIVEPT1   0.538469310105683 
#define FIVEPT2   0.906179845938664

#define FIVEPTWT1    0.478628670499366
#define FIVEPTWT2    0.236926885056189
#define FIVEPTWTZERO 0.568888888888889


// Six point gaussian quadrature



#define SIXPT1  0.238619186083197
#define SIXPT2  0.661209386466265
#define SIXPT3  0.932469514203152

#define SIXPTWT1 0.467913934572691
#define SIXPTWT2 0.360761573048139
#define SIXPTWT3 0.171324492379170

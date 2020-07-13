#ifndef HAVE_SHAPESTRUCT_HPP
#include "../datastruct/shapestruct.hpp"
#endif

#ifndef HAVE_INTEGCONSTS_H
#include "integconsts.h"
#endif

#include <math.h>

/*---------------------------------------------------------------
  Nachiket Gokhale gokhalen@bu.edu 
  Two Dimensional UNIFORM Gaussian Quadrature
  There are ninte[0] quadrature points in each direction
  These integration rules can be found in Abramowitz and Stegan 
  ----------------------------------------------------------------*/

int gauss2d(shapestruct *shpstrct){
  
  /*1X1 gaussian quadrature*/
  if(shpstrct->ninte[0] == 1){
    
    shpstrct->intpts[0][0] = 0;  //  natural(x) coord of only integration point
    shpstrct->intpts[0][1] = 0;  //  natural(y) coord of only integration point 
    shpstrct->wts[0] = 4;        //  set weight
  }

  /*2X2 gaussian quadrature*/
  if(shpstrct->ninte[0] == 2){
    /*set x and y coordinates */
    shpstrct->intpts[0][0] = -ONE_O_SQRT3;    shpstrct->intpts[0][1] = -ONE_O_SQRT3;
    shpstrct->intpts[1][0] =  ONE_O_SQRT3;    shpstrct->intpts[1][1] = -ONE_O_SQRT3;
    shpstrct->intpts[2][0] =  ONE_O_SQRT3;    shpstrct->intpts[2][1] =  ONE_O_SQRT3;
    shpstrct->intpts[3][0] = -ONE_O_SQRT3;    shpstrct->intpts[3][1] =  ONE_O_SQRT3;

    // all the weights are = 1
    shpstrct->wts[0] = shpstrct->wts[1] = shpstrct->wts[2] = shpstrct->wts[3] = 1;  

  }
  
  /*3X3 gaussian quadrature*/
  if(shpstrct->ninte[0] == 3){
    /*natural coordinates of integration points*/
    shpstrct->intpts[0][0] = -SQRT3BY5;    shpstrct->intpts[0][1] = -SQRT3BY5;
    shpstrct->intpts[1][0] =      ZERO;    shpstrct->intpts[1][1] = -SQRT3BY5;
    shpstrct->intpts[2][0] =  SQRT3BY5;    shpstrct->intpts[2][1] = -SQRT3BY5;
    shpstrct->intpts[3][0] = -SQRT3BY5;    shpstrct->intpts[3][1] =      ZERO;
    shpstrct->intpts[4][0] =      ZERO;    shpstrct->intpts[4][1] =      ZERO;
    shpstrct->intpts[5][0] =  SQRT3BY5;    shpstrct->intpts[5][1] =      ZERO;
    shpstrct->intpts[6][0] = -SQRT3BY5;    shpstrct->intpts[6][1] =  SQRT3BY5;
    shpstrct->intpts[7][0] =     ZERO ;    shpstrct->intpts[7][1] =  SQRT3BY5;
    shpstrct->intpts[8][0] =  SQRT3BY5;    shpstrct->intpts[8][1] =  SQRT3BY5;

    /*weights of integration points*/ 
    shpstrct->wts[0] =  W15BY9*W15BY9;

    shpstrct->wts[1] =  W15BY9*W28BY9;

    shpstrct->wts[2] =  W15BY9*W35BY9;

    shpstrct->wts[3] =  W28BY9*W15BY9;

    shpstrct->wts[4] =  W28BY9*W28BY9;

    shpstrct->wts[5] =  W28BY9*W35BY9;

    shpstrct->wts[6] =  W35BY9*W15BY9;

    shpstrct->wts[7] =  W35BY9*W28BY9;

    shpstrct->wts[8] =  W35BY9*W35BY9;

  }
  
  /*4X4 gaussian quadrature*/
  if(shpstrct->ninte[0] == 4){
    
    shpstrct->intpts[0][0]  = -FOURPT2;  shpstrct->intpts[0][1]  = -FOURPT2;
    shpstrct->intpts[1][0]  = -FOURPT2;  shpstrct->intpts[1][1]  = -FOURPT1;
    shpstrct->intpts[2][0]  = -FOURPT2;  shpstrct->intpts[2][1]  =  FOURPT1;
    shpstrct->intpts[3][0]  = -FOURPT2;  shpstrct->intpts[3][1]  =  FOURPT2;

    shpstrct->intpts[4][0]  = -FOURPT1;  shpstrct->intpts[4][1]  = -FOURPT2;
    shpstrct->intpts[5][0]  = -FOURPT1;  shpstrct->intpts[5][1]  = -FOURPT1;
    shpstrct->intpts[6][0]  = -FOURPT1;  shpstrct->intpts[6][1]  =  FOURPT1;
    shpstrct->intpts[7][0]  = -FOURPT1;  shpstrct->intpts[7][1]  =  FOURPT2;

    shpstrct->intpts[8][0]  =  FOURPT1;  shpstrct->intpts[8][1]  = -FOURPT2;
    shpstrct->intpts[9][0]  =  FOURPT1;  shpstrct->intpts[9][1]  = -FOURPT1;
    shpstrct->intpts[10][0] =  FOURPT1;  shpstrct->intpts[10][1] =  FOURPT1;
    shpstrct->intpts[11][0] =  FOURPT1;  shpstrct->intpts[11][1] =  FOURPT2;

    shpstrct->intpts[12][0] =  FOURPT2;  shpstrct->intpts[12][1] = -FOURPT2;
    shpstrct->intpts[13][0] =  FOURPT2;  shpstrct->intpts[13][1] = -FOURPT1;
    shpstrct->intpts[14][0] =  FOURPT2;  shpstrct->intpts[14][1] =  FOURPT1;
    shpstrct->intpts[15][0] =  FOURPT2;  shpstrct->intpts[15][1] =  FOURPT2;
    
    shpstrct->wts[0]   = FOURPTWT2*FOURPTWT2;
    shpstrct->wts[1]   = FOURPTWT2*FOURPTWT1;
    shpstrct->wts[2]   = FOURPTWT2*FOURPTWT1;
    shpstrct->wts[3]   = FOURPTWT2*FOURPTWT2;

    shpstrct->wts[4]   = FOURPTWT1*FOURPTWT2;
    shpstrct->wts[5]   = FOURPTWT1*FOURPTWT1;
    shpstrct->wts[6]   = FOURPTWT1*FOURPTWT1;
    shpstrct->wts[7]   = FOURPTWT1*FOURPTWT2;
    
    shpstrct->wts[8]   = FOURPTWT1*FOURPTWT2;
    shpstrct->wts[9]   = FOURPTWT1*FOURPTWT1;
    shpstrct->wts[10]  = FOURPTWT1*FOURPTWT1;
    shpstrct->wts[11]  = FOURPTWT1*FOURPTWT2;

    shpstrct->wts[12]  = FOURPTWT2*FOURPTWT2;
    shpstrct->wts[13]  = FOURPTWT2*FOURPTWT1;
    shpstrct->wts[14]  = FOURPTWT2*FOURPTWT1;
    shpstrct->wts[15]  = FOURPTWT2*FOURPTWT2;
    
    
  }
  
  if(shpstrct->ninte[0]==5){
    
    shpstrct->intpts[0][0]   = -FIVEPT2;  shpstrct->intpts[0][1]   = -FIVEPT2;
    shpstrct->intpts[1][0]   = -FIVEPT2;  shpstrct->intpts[1][1]   = -FIVEPT1;
    shpstrct->intpts[2][0]   = -FIVEPT2;  shpstrct->intpts[2][1]   =     ZERO;
    shpstrct->intpts[3][0]   = -FIVEPT2;  shpstrct->intpts[3][1]   =  FIVEPT1;
    shpstrct->intpts[4][0]   = -FIVEPT2;  shpstrct->intpts[4][1]   =  FIVEPT2;

    shpstrct->intpts[5][0]   = -FIVEPT1;  shpstrct->intpts[5][1]   = -FIVEPT2;
    shpstrct->intpts[6][0]   = -FIVEPT1;  shpstrct->intpts[6][1]   = -FIVEPT1;
    shpstrct->intpts[7][0]   = -FIVEPT1;  shpstrct->intpts[7][1]   =     ZERO;
    shpstrct->intpts[8][0]   = -FIVEPT1;  shpstrct->intpts[8][1]   =  FIVEPT1;
    shpstrct->intpts[9][0]   = -FIVEPT1;  shpstrct->intpts[9][1]   =  FIVEPT2;

    shpstrct->intpts[10][0]  =     ZERO;  shpstrct->intpts[10][1]  = -FIVEPT2;
    shpstrct->intpts[11][0]  =     ZERO;  shpstrct->intpts[11][1]  = -FIVEPT1;
    shpstrct->intpts[12][0]  =     ZERO;  shpstrct->intpts[12][1]  =     ZERO;
    shpstrct->intpts[13][0]  =     ZERO;  shpstrct->intpts[13][1]  =  FIVEPT1;
    shpstrct->intpts[14][0]  =     ZERO;  shpstrct->intpts[14][1]  =  FIVEPT2;

    shpstrct->intpts[15][0]  =  FIVEPT1;  shpstrct->intpts[15][1]  = -FIVEPT2;
    shpstrct->intpts[16][0]  =  FIVEPT1;  shpstrct->intpts[16][1]  = -FIVEPT1;
    shpstrct->intpts[17][0]  =  FIVEPT1;  shpstrct->intpts[17][1]  =     ZERO;
    shpstrct->intpts[18][0]  =  FIVEPT1;  shpstrct->intpts[18][1]  =  FIVEPT1;
    shpstrct->intpts[19][0]  =  FIVEPT1;  shpstrct->intpts[19][1]  =  FIVEPT2;

    shpstrct->intpts[20][0]  =  FIVEPT2;  shpstrct->intpts[20][1]  = -FIVEPT2;
    shpstrct->intpts[21][0]  =  FIVEPT2;  shpstrct->intpts[21][1]  = -FIVEPT1;
    shpstrct->intpts[22][0]  =  FIVEPT2;  shpstrct->intpts[22][1]  =     ZERO;
    shpstrct->intpts[23][0]  =  FIVEPT2;  shpstrct->intpts[23][1]  =  FIVEPT1;
    shpstrct->intpts[24][0]  =  FIVEPT2;  shpstrct->intpts[24][1]  =  FIVEPT2;

    shpstrct->wts[0]    = FIVEPTWT2*FIVEPTWT2;
    shpstrct->wts[1]    = FIVEPTWT2*FIVEPTWT1;
    shpstrct->wts[2]    = FIVEPTWT2*FIVEPTWTZERO;
    shpstrct->wts[3]    = FIVEPTWT2*FIVEPTWT1;
    shpstrct->wts[4]    = FIVEPTWT2*FIVEPTWT2;

    shpstrct->wts[5]    = FIVEPTWT1*FIVEPTWT2;
    shpstrct->wts[6]    = FIVEPTWT1*FIVEPTWT1;
    shpstrct->wts[7]    = FIVEPTWT1*FIVEPTWTZERO;
    shpstrct->wts[8]    = FIVEPTWT1*FIVEPTWT1;
    shpstrct->wts[9]    = FIVEPTWT1*FIVEPTWT2;

    shpstrct->wts[10]   = FIVEPTWTZERO*FIVEPTWT2;
    shpstrct->wts[11]   = FIVEPTWTZERO*FIVEPTWT1;
    shpstrct->wts[12]   = FIVEPTWTZERO*FIVEPTWTZERO;
    shpstrct->wts[13]   = FIVEPTWTZERO*FIVEPTWT1;
    shpstrct->wts[14]   = FIVEPTWTZERO*FIVEPTWT2;

    shpstrct->wts[15]   = FIVEPTWT1*FIVEPTWT2;
    shpstrct->wts[16]   = FIVEPTWT1*FIVEPTWT1;
    shpstrct->wts[17]   = FIVEPTWT1*FIVEPTWTZERO;
    shpstrct->wts[18]   = FIVEPTWT1*FIVEPTWT1;
    shpstrct->wts[19]   = FIVEPTWT1*FIVEPTWT2;

    shpstrct->wts[20]   = FIVEPTWT2*FIVEPTWT2;
    shpstrct->wts[21]   = FIVEPTWT2*FIVEPTWT1;
    shpstrct->wts[22]   = FIVEPTWT2*FIVEPTWTZERO;
    shpstrct->wts[23]   = FIVEPTWT2*FIVEPTWT1;
    shpstrct->wts[24]   = FIVEPTWT2*FIVEPTWT2;

  }
  

  if(shpstrct->ninte[0] == 6){
    
    shpstrct->intpts[0][0]   = -SIXPT3;  shpstrct->intpts[0][1]  = -SIXPT3;
    shpstrct->intpts[1][0]   = -SIXPT3;  shpstrct->intpts[1][1]  = -SIXPT2;
    shpstrct->intpts[2][0]   = -SIXPT3;  shpstrct->intpts[2][1]  = -SIXPT1;
    shpstrct->intpts[3][0]   = -SIXPT3;  shpstrct->intpts[3][1]  =  SIXPT1;
    shpstrct->intpts[4][0]   = -SIXPT3;  shpstrct->intpts[4][1]  =  SIXPT2;
    shpstrct->intpts[5][0]   = -SIXPT3;  shpstrct->intpts[5][1]  =  SIXPT3;

    shpstrct->intpts[6][0]   = -SIXPT2;  shpstrct->intpts[6][1]  = -SIXPT3;
    shpstrct->intpts[7][0]   = -SIXPT2;  shpstrct->intpts[7][1]  = -SIXPT2;
    shpstrct->intpts[8][0]   = -SIXPT2;  shpstrct->intpts[8][1]  = -SIXPT1;
    shpstrct->intpts[9][0]   = -SIXPT2;  shpstrct->intpts[9][1]  =  SIXPT1;
    shpstrct->intpts[10][0]  = -SIXPT2;  shpstrct->intpts[10][1]  =  SIXPT2;
    shpstrct->intpts[11][0]  = -SIXPT2;  shpstrct->intpts[11][1]  =  SIXPT3;

    shpstrct->intpts[12][0]  = -SIXPT1;  shpstrct->intpts[12][1]  = -SIXPT3;
    shpstrct->intpts[13][0]  = -SIXPT1;  shpstrct->intpts[13][1]  = -SIXPT2;
    shpstrct->intpts[14][0]  = -SIXPT1;  shpstrct->intpts[14][1]  = -SIXPT1;
    shpstrct->intpts[15][0]  = -SIXPT1;  shpstrct->intpts[15][1]  =  SIXPT1;
    shpstrct->intpts[16][0]  = -SIXPT1;  shpstrct->intpts[16][1]  =  SIXPT2;
    shpstrct->intpts[17][0]  = -SIXPT1;  shpstrct->intpts[17][1]  =  SIXPT3;

    shpstrct->intpts[18][0]  = SIXPT1;  shpstrct->intpts[18][1]   = -SIXPT3;
    shpstrct->intpts[19][0]  = SIXPT1;  shpstrct->intpts[19][1]   = -SIXPT2;
    shpstrct->intpts[20][0]  = SIXPT1;  shpstrct->intpts[20][1]   = -SIXPT1;
    shpstrct->intpts[21][0]  = SIXPT1;  shpstrct->intpts[21][1]   =  SIXPT1;
    shpstrct->intpts[22][0]  = SIXPT1;  shpstrct->intpts[22][1]   =  SIXPT2;
    shpstrct->intpts[23][0]  = SIXPT1;  shpstrct->intpts[23][1]   =  SIXPT3;

    shpstrct->intpts[24][0]  = SIXPT2;  shpstrct->intpts[24][1]   = -SIXPT3;
    shpstrct->intpts[25][0]  = SIXPT2;  shpstrct->intpts[25][1]   = -SIXPT2;
    shpstrct->intpts[26][0]  = SIXPT2;  shpstrct->intpts[26][1]   = -SIXPT1;
    shpstrct->intpts[27][0]  = SIXPT2;  shpstrct->intpts[27][1]   =  SIXPT1;
    shpstrct->intpts[28][0]  = SIXPT2;  shpstrct->intpts[28][1]   =  SIXPT2;
    shpstrct->intpts[29][0]  = SIXPT2;  shpstrct->intpts[29][1]   =  SIXPT3;

    shpstrct->intpts[30][0]  = SIXPT3;  shpstrct->intpts[30][1]   = -SIXPT3;
    shpstrct->intpts[31][0]  = SIXPT3;  shpstrct->intpts[31][1]   = -SIXPT2;
    shpstrct->intpts[32][0]  = SIXPT3;  shpstrct->intpts[32][1]   = -SIXPT1;
    shpstrct->intpts[33][0]  = SIXPT3;  shpstrct->intpts[33][1]   =  SIXPT1;
    shpstrct->intpts[34][0]  = SIXPT3;  shpstrct->intpts[34][1]   =  SIXPT2;
    shpstrct->intpts[35][0]  = SIXPT3;  shpstrct->intpts[35][1]   =  SIXPT3;

    shpstrct->wts[0]    = SIXPTWT3*SIXPTWT3;
    shpstrct->wts[1]    = SIXPTWT3*SIXPTWT2;
    shpstrct->wts[2]    = SIXPTWT3*SIXPTWT1;
    shpstrct->wts[3]    = SIXPTWT3*SIXPTWT1;
    shpstrct->wts[4]    = SIXPTWT3*SIXPTWT2;
    shpstrct->wts[5]    = SIXPTWT3*SIXPTWT3;

    shpstrct->wts[6]    = SIXPTWT2*SIXPTWT3;
    shpstrct->wts[7]    = SIXPTWT2*SIXPTWT2;
    shpstrct->wts[8]    = SIXPTWT2*SIXPTWT1;
    shpstrct->wts[9]    = SIXPTWT2*SIXPTWT1;
    shpstrct->wts[10]    = SIXPTWT2*SIXPTWT2;
    shpstrct->wts[11]    = SIXPTWT2*SIXPTWT3;

    shpstrct->wts[12]    = SIXPTWT1*SIXPTWT3;
    shpstrct->wts[13]    = SIXPTWT1*SIXPTWT2;
    shpstrct->wts[14]    = SIXPTWT1*SIXPTWT1;
    shpstrct->wts[15]    = SIXPTWT1*SIXPTWT1;
    shpstrct->wts[16]    = SIXPTWT1*SIXPTWT2;
    shpstrct->wts[17]    = SIXPTWT1*SIXPTWT3;

    shpstrct->wts[18]    = SIXPTWT1*SIXPTWT3;
    shpstrct->wts[19]    = SIXPTWT1*SIXPTWT2;
    shpstrct->wts[20]    = SIXPTWT1*SIXPTWT1;
    shpstrct->wts[21]    = SIXPTWT1*SIXPTWT1;
    shpstrct->wts[22]    = SIXPTWT1*SIXPTWT2;
    shpstrct->wts[23]    = SIXPTWT1*SIXPTWT3;

    shpstrct->wts[24]    = SIXPTWT2*SIXPTWT3;
    shpstrct->wts[25]    = SIXPTWT2*SIXPTWT2;
    shpstrct->wts[26]    = SIXPTWT2*SIXPTWT1;
    shpstrct->wts[27]    = SIXPTWT2*SIXPTWT1;
    shpstrct->wts[28]    = SIXPTWT2*SIXPTWT2;
    shpstrct->wts[29]    = SIXPTWT2*SIXPTWT3;

    shpstrct->wts[30]    = SIXPTWT3*SIXPTWT3;
    shpstrct->wts[31]    = SIXPTWT3*SIXPTWT2;
    shpstrct->wts[32]    = SIXPTWT3*SIXPTWT1;
    shpstrct->wts[33]    = SIXPTWT3*SIXPTWT1;
    shpstrct->wts[34]    = SIXPTWT3*SIXPTWT2;
    shpstrct->wts[35]    = SIXPTWT3*SIXPTWT3;
    
    
    
  }

  return 0;
}

#ifndef HAVE_INTEG_H
#include "integ.h"
#endif

#ifndef HAVE_SHAPESTRUCT_HPP
#include "../datastruct/shapestruct.hpp"
#endif
#include <math.h>

#ifndef HAVE_INTEGCONSTS_H
#include "integconsts.h"
#endif

int gauss3d(shapestruct *shpstrct){


  // 1X1X1 gaussian quadrature
  if(shpstrct->ninte[0] == 1){

    shpstrct->intpts[0][0] = 0; // x-coordinate
    shpstrct->intpts[0][1] = 0; // y-coordinate
    shpstrct->intpts[0][2] = 0; // z-coordinate
    
    shpstrct->wts[0]       = 8;       
  }


  if(shpstrct->ninte[0] == 2){
    
    // check this later. derived in about 2 mins 
    shpstrct->intpts[0][0] = -ONE_O_SQRT3; shpstrct->intpts[0][1] = -ONE_O_SQRT3; shpstrct->intpts[0][2] = -ONE_O_SQRT3;
    shpstrct->intpts[1][0] = -ONE_O_SQRT3; shpstrct->intpts[1][1] = -ONE_O_SQRT3; shpstrct->intpts[1][2] =  ONE_O_SQRT3;
    shpstrct->intpts[2][0] = -ONE_O_SQRT3; shpstrct->intpts[2][1] =  ONE_O_SQRT3; shpstrct->intpts[2][2] = -ONE_O_SQRT3;
    shpstrct->intpts[3][0] = -ONE_O_SQRT3; shpstrct->intpts[3][1] =  ONE_O_SQRT3; shpstrct->intpts[3][2] =  ONE_O_SQRT3;
    shpstrct->intpts[4][0] =  ONE_O_SQRT3; shpstrct->intpts[4][1] = -ONE_O_SQRT3; shpstrct->intpts[4][2] = -ONE_O_SQRT3;
    shpstrct->intpts[5][0] =  ONE_O_SQRT3; shpstrct->intpts[5][1] = -ONE_O_SQRT3; shpstrct->intpts[5][2] =  ONE_O_SQRT3;
    shpstrct->intpts[6][0] =  ONE_O_SQRT3; shpstrct->intpts[6][1] =  ONE_O_SQRT3; shpstrct->intpts[6][2] = -ONE_O_SQRT3;
    shpstrct->intpts[7][0] =  ONE_O_SQRT3; shpstrct->intpts[7][1] =  ONE_O_SQRT3; shpstrct->intpts[7][2] =  ONE_O_SQRT3;

    for(int iwts = 0;iwts<8;iwts++){
      shpstrct->wts[iwts] = 1;
    }
    
  }

  // three point integration
  
  if(shpstrct->ninte[0] == 3){

    double intcoord[3]; // integration point coordinates
    intcoord[0] = -SQRT3BY5; intcoord[1] = 0.0; intcoord[2] = SQRT3BY5;
    
    double intwts[3]; 
    intwts[0]= W15BY9; intwts[1] = W28BY9; intwts[2]=W35BY9;
    
    int intecounter = 0;
    
    
    for(int zinte = 0; zinte<shpstrct->ninte[0];zinte++){
      for(int yinte = 0; yinte <shpstrct->ninte[0] ; yinte++ ){
	for(int xinte = 0; xinte <shpstrct->ninte[0] ; xinte++ ){
	  shpstrct->intpts[intecounter][0] = intcoord[xinte];
	  shpstrct->intpts[intecounter][1] = intcoord[yinte];
	  shpstrct->intpts[intecounter][2] = intcoord[zinte];
	  shpstrct->wts[intecounter]       = intwts[xinte]*intwts[yinte]*intwts[zinte];
	  intecounter++;
	}
      }
    }

    assert(intecounter == shpstrct->tinte);

  }
  
  

  return 0;
}

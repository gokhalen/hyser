#ifndef HAVE_SHAPE_H
#define HAVE_SHAPE_H
#endif

#ifndef HAVE_SHAPESTRUCT_HPP
#include "../datastruct/shapestruct.hpp"
#endif

#ifndef HAVE_UTILS_HPP
#include "../utillib/utils.hpp"
#endif

// function declarations of shape function routines
// Nachiket Gokhale gokhalen@bu.edu Spring 2005

// function shape is calls other functions depending on dimension 

int shape   (shapestruct *shpstrct);
int shape1d (shapestruct *shpstrct);
int shape2d (shapestruct *shpstrct);
int shape2d9(shapestruct *shpstrct);
int shape3d (shapestruct *shpstrct);

void checkjacobiandet(shapestruct*, int );

#ifndef HAVE_INTEG_H
#define HAVE_INTEG_H
#endif

#ifndef HAVE_SHAPESTRUCT_HPP
#include "../datastruct/shapestruct.hpp"
#endif

#include <cassert>
using namespace std;

// integ decides which function to call, depending on the element type
int integ(shapestruct *shpstrct);

// support for 1 to 3 point gaussian integration.
int gauss1d(shapestruct *shpstrct); 
int gauss2d(shapestruct *shpstrct);
int gauss3d(shapestruct *shpstrct);

/*
 * Nachiket Gokhale gokhalen@bu.edu
 * 3rd April 2005
 * 
 *
 */

#ifndef HAVE_GLOBALS_H
#define HAVE_GLOBALS_H
#endif

// for open bufferes
#define MAX_BUFFER 4096

// general purpose epsilon value
#define DOUBLE_EPSILON   1E-12
#define LINEAR_SYS_REL_TOL 1E-10

// general purpose precision
#define DEFAULT_DOUBLE_PRECISION 16

// flags used to call localizer to localize the current state data
// or the last state data
#define CURRENT_LOAD_ITERATION 0
#define LAST_LOAD_ITERATION 1
#define ENHANCED 2

// output file formats 
// OUTPUT_STANDARD this is the format used in my phd
// OUTPUT_GNUPLOT the first ndime fields are coordinates, then the values
//                this is useful for reading and plotting in octave using
//                "splot 'sol0-0-0.out' using 1:2:3 with image"
#define OUTPUT_STANDARD 0  
#define OUTPUT_GNUPLOT  1  

#define SOLVE_TYPE_BVP            1   // straight boundary value problem bvp solve
#define SOLVE_TYPE_EIGENK         2   // eigen-analysis  Kx = wx  
#define SOLVE_TYPE_EIGENM         3   // eigen-analysis  Mx = wx
#define SOLVE_TYPE_EIGEN          4   // eigen-analysis  Mx = wMx
#define SOLVE_TYPE_DISPERSIONK2W  5   // dispersion analyis - hussein formulation k->w
#define SOLVE_TYPE_DISPERSIONW2K  6   // dispersion analyis - formulation w->k
#define SOLVE_TYPE_HOMOGENIZATION 7   // homogenization analysis

#define BC_TYPE_DIRICHLET      1
#define BC_TYPE_TRACTION       2
#define BC_TYPE_POINTFORC      3
#define BC_TYPE_PERIODIC       4

// common to elem_base and globals - must be consistent
#define BUILD_REAL_STIFFNESS        1 
#define BUILD_COMPLEX_STIFFNESS_K2W 2
#define BUILD_COMPLEX_STIFFNESS_W2K 3

// common to elem_base and globals - must be consistent
#define BUILD_CONSISTENT_MASS 1
#define BUILD_LUMPED_MASS 2 

#define PHYSICS_ELAS   1 
#define PHYSICS_PIEZO  2 

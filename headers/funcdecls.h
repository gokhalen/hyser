/**********************************************************************************
 *  Function declarations for Hyser, a serial hyperelastic finite element program. 
 *  Nachiket Gokhale gokhalen@bu.edu
 *
 *
 **********************************************************************************/

#ifndef HAVE_FUNCDECLS_H
#define HAVE_FUNCDECLS_H
#endif

//#ifndef HAVE_PETSC_HEADERS_H
//#include "petsc_headers.h"
//#endif

//#ifndef HAVE_TAO_HEADERS_H
//#include "tao_headers.h"
//#endif

#ifndef HAVE_ELEMLIB_HPP
#include "../lib/elemlib/elemlib.hpp"
#endif

#ifndef HAVE_MESH_HPP
#include "../lib/meshlib/mesh.hpp"
#endif

// mesh reading functions
// forward declaration, break the circular include
class mesh; 

// building stiffness
int buildtan(mesh *inmesh);

// building residual vector
int buildres(mesh *inmesh);

// build the enhanced parameters
int buildenh(mesh *inmesh);

// assembling stiffness
int assembletan(mesh* inmesh,long int ielem,elem_base* felem);

// assembling residual
int assembleres(mesh* inmesh,long int ielem,elem_base* felem);

// assembling enhanced modes
int assembleenh(mesh* inmesh, int ielem, elem_base* felem);

//get elemref: returns a reference to the base class of elements
elem_base* getelemref(mesh *inmesh, int, elem_base*);

//localize: puts data for element number ielem into appropriate class
elem_base* localizer(mesh *inmesh, int locielem,elem_base *felem, int loaditflag); 

//csr utils: 
int csr_add(double *val,long int *row,long int *col,long int nn);


// the forward calculation
int calc_forward(mesh *inmesh);

// the dual field calculation
int calc_dual(mesh *inmesh);




// finite differences gradient calculation
int fdgrad(mesh *inmesh);

// calculate functional and gradient
int calcfuncgrad(mesh *inmesh);

// function for homogenization
int homogenize(mesh *inmesh);

// calculate rhs for the dual
int builddualrhs(mesh *inmesh);


// the TAO wrapper to calculate the functional and the gradient.
//PetscErrorCode gradientcallbackfortao(TaoSolver taosol, Vec X, double *f,Vec G,void *ptr);

// the monitor to print out the properties at every iteration.
//PetscErrorCode propMonitor(TaoSolver taoapp,void *ptr);



// the function to calculate the norm of the measured field
// this is called only once
int calc_normmeas(mesh *);

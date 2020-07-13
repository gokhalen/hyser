/*-------------------------------------------------------------------------
 * Nachiket Gokhale gokhalen@bu.edu
 * 
 * A constitutive equation library.
 * 
 * Descripton:
 * In non-linear elasticity, there are infinite number of laws
 * that can relate stress to strain. It is virtually impossible to 
 * code this "with in an element routine" on an ad hoc basis.            
 * Nor is it practical to create different "element routines"                  
 * for each element. Therefore, this constitutive law library provides 
 * according to law number: 
 *
 *  1) The material stiffness (C_IJKL) at a given point. 
 *  2) Elemental contribution to the tangent stiffness.
 *  3) The elemental contribution to the gradient.
 * 
 * Documentation: As of May 27, 2005 these routines are no longer documented in
 * the docs directory. They are documented in the appendix of the dissertation.  * In much more detail with all the steps for the derivation.
 *
 * Data structures:
 *                 This also defines a few data structures that help passinng
 * arguments a lot easier. Maybe I'll change it to a class, a nice inherited 
 * structure 
 *                
 *------------------------------------------------------------------------*/
#ifndef HAVE_CONSLIB_H
#define HAVE_CONSLIB_H
#endif

#ifndef HAVE_UTILS_HPP
#include "../utillib/utils.hpp"
#endif

typedef struct {

  // the input parameters
  int grdparm;       // which parameter to take the gradient wrt.
  double matdisp;    // material displacement U(X)
  double matstrain;  // strain U(X) 
  double *stfparam;  // pointer to an array containing "stiffness parameters"

  // the output parameters
  
  double stress;     // matstrain*stiffness
  double stiffness;  // output parameter
  double ctang;      // contribution to the consistent tangent
  double gradcontrib; //contribution to the gradient (doesn't include shape func)
} conselas1d;


typedef struct {
  // we play little type games with pointers here (in elas2d and elas3d). 
  // passing arrays with pointers has to be done so carefully. 
  
  
  int ndime;         // dimension of the problem, 2/3d
  int grdparm;       // which parameter to take the gradient wrt.
  int grdflag;       // check if we are required to calculate the gradient

  // the input parameters
  double *disp;      // pointer to the first element of the displacement vector
  double *stfparam;  // pointer to an array containing "stiffness parameters"



  double **rcgreen;  
  double **defgrad;

  // the output parameters

  double ctang;            
  double ****cijklptr;    // (must be initialized, used to modify external ptr)
  double ****aijklptr;    // (must be inititalized used to modify external ptr)
  double ****sptan;       // spatial tangent modulus

  // the stress tensors 
  double **cauchystrs;
  double **krchstrs;
  double **firstpk;
  double **secpk;         // (must be initialized, used to modify external ptr)


  double **pk2grad;       // (gradient of the second piola kirchoff)
  double gradcontrib;     // stores the gradient contrib from the gradparam


  // the follow
  // the tangent modulii
  
  double **inclambda;   // the incremental \lambda
  double **incmu;       // the incremental \mu
  
  //matrix representation of aijkl cijkl
  double **AIJ;
  double **CIJ;
  

} conselas23d;

typedef struct {
  
} conheat2d;


 
int elas1d(int lawnumber, conselas1d *);
int elas2d(int lawnumber, conselas23d *);
int elas3d(int lawnumber, conselas23d *);
int heat2d(int lawnumber, conheat2d  *);

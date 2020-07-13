#ifndef HAVE_PETSC_SYSTEM_HPP
#define HAVE_PETSC_SYSTEM_HPP
#endif

#ifndef HAVE_EQUATION_SYSTEM_HPP
#include "equation_system.hpp"
#endif

#ifndef HAVE_PETSC_HEADERS_H
#include "../../headers/petsc_headers.h"
#endif

class petsc_system:public equation_system{

public:  
  Mat stiff;  
  Vec rhsvec;  
  Vec solvec;  
  PetscViewer view_mat,view_rhs,view_sol;
  PetscErrorCode ierr;

  
  
  petsc_system();
  petsc_system(equation_initdata temp);
  petsc_system(int num_equations_arg, int which_field_arg, int inewton, int iload,int printflag);
  virtual ~petsc_system();
  virtual void createcsrk(int *row_ptr,int *col_ptr,double *value_arr);
  virtual void createrhs(double **value_arr);
  virtual int solve();
  virtual void getsol(double **);
  
};

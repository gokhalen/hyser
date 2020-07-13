/**
 * Nachiket Gokhale gokhalen@bu.edu
 * This class, equation_system, defines a consistent interface
 * to use different solvers. Currently it uses PETSc Only
 * This is a declaration, the definition is in the .cpp file
 */

#ifndef HAVE_EQUATION_SYSTEM_HPP
#define HAVE_EQUATION_SYSTEM_HPP
#endif

#ifndef HAVE_CSR_HPP
#include "../utillib/csr.hpp"
#endif

#ifndef HAVE_GLOBALS_H
#include "../../headers/globals.h"
#endif

class equation_initdata{
public:
  int num_equations;
  int num_equations_real;
  int num_equations_imag;
  int which_field;
  int inewton;
  int iload;
  int printflag;
  int max_nnz_per_row;
  int which_seg;
  int which_pt;
};

class equation_system{

 public:
  int num_equations; 
  int num_equations_real;
  int num_equations_imag; 
  int max_nnz_per_row;

  int which_field;   // number of field that is being solved.

  // use these to print out matrices
  int inewton;
  int iload;
  int printflag;
  equation_initdata eqninitdata;

  equation_system();
  equation_system(equation_initdata temp);
  equation_system(int num_equations_arg, int which_field_arg, int inewton,int iload,int printflag);

  virtual void createcsrk(int *row_ptr,int *col_ptr, double *value_arr);
  virtual void createcsrkreal(int *row_ptr,int *col_ptr,double *value_arr);
  virtual void createcsrkimag(int *row_ptr,int *col_ptr,double *value_arr);
  virtual void createcsrm(int *row_ptr,int *col_ptr, double *value_arr);
  virtual void createcsrmreal(int *row_ptr,int *col_ptr, double *value_arr);
  virtual void createcsrmimag(int *row_ptr,int *col_ptr, double *value_arr);

  virtual void createrhs(double **value);
  virtual int  solve();
  virtual int  solvecomplex();

  virtual void getsol(double **);
  virtual int  geteigenvalue(double*,double*,double*,int);
  virtual int  geteigenvector(double **, double **,int);
  virtual int  geteigenvectordispersion(double **,double **,int);
  virtual int  getconveigen();


  virtual ~equation_system();
};

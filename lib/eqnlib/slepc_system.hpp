#ifndef HAVE_SLEPC_SYSTEM_HPP
#define HAVE_SLEPC_SYSTEM_HPP
#endif

#ifndef HAVE_EQUATION_SYSTEM_HPP
#include "equation_system.hpp"
#endif

#ifndef HAVE_SLEPC_HEADERS_H
#include "../../headers/slepc_headers.h"
#endif

class slepc_system:public equation_system{

public:  
  Mat stiff,mass;  
  Mat stiff_real,stiff_imag,stiff_complex;
  Mat mass_real,mass_imag,mass_complex;
  //Vec rhsvec;  
  Vec solvec;  
  PetscViewer view_mat,view_rhs,view_sol;
  PetscErrorCode ierr;
  PetscInt nconv;
  int num_eigen;
  EPS eps;
  EPSType epstype;
  Vec   eigvalreal,eigvalimag,eigenerror;

  slepc_system();
  slepc_system(equation_initdata eqninitdata);
  slepc_system(int num_equations_arg, int which_field_arg, int inewton, int iload,int printflag);
  virtual ~slepc_system();
  virtual void createcsrk(int *row_ptr,int *col_ptr,double *value_arr);
  virtual void createcsrkreal(int *row_ptr,int *col_ptr,double *value_arr);
  virtual void createcsrkimag(int *row_ptr,int *col_ptr,double *value_arr);
  virtual void createcsrm(int *row_ptr,int *col_ptr,double *value_arr);
  virtual void createcsrmreal(int *row_ptr,int *col_ptr,double *value_arr);
  virtual void createcsrmimag(int *row_ptr,int *col_ptr,double *value_arr);
  virtual void createrhs(double **value_arr);
  virtual int solve();
  virtual int solvecomplex();
  virtual void getsol(double **);
  virtual int geteigenvalue(double *, double *, double *, int);
  virtual int geteigenvector(double **, double **,int);
  virtual int geteigenvectordispersion(double **,double **,int);
  virtual void setsolvetype(int);
  virtual int getconveigen();
  

private:
  int allocmass, allocstiff;
  int allocmassimag, allocstiffimag;
  int solvetype;
  virtual int assemblecomplexmat(Mat *MatReal, Mat *MatImag, Mat *MatComplex);

};

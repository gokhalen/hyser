#ifndef HAVE_SLEPC_QUAD_SYSTEM_HPP
#define HAVE_SLEPC_QUAD_SYSTEM_HPP
#endif

#ifndef HAVE_EQUATION_SYSTEM_HPP
#include "equation_system.hpp"
#endif

#ifndef HAVE_SLEPC_HEADERS_H
#include "../../headers/slepc_headers.h"
#endif

class slepc_quad_system:public equation_system{
public:
  Mat MM,CC,KK;
  Mat MMr,CCr,KKr1,KKr2;
  Mat MMi,CCi,KKi;
  Vec solvec;
  PetscViewer view_mat,view_rhs,view_sol;
  PetscErrorCode ierr;
  PetscInt nconv;
  int num_eigen;
  int which_pt,which_seg;
  QEP qep;
  QEPType type;
  Vec eigvalreal,eigvalimag,eigenerror;
  PetscScalar ww;
  slepc_quad_system();
  slepc_quad_system(equation_initdata eqninitdata);
  slepc_quad_system(int num_equations_arg, int which_field_arg, int inewton, int iload,int printflag);

  virtual ~slepc_quad_system();

  virtual void createcsrkr1(int *,int *,double *);
  virtual void createcsrkr2(int *,int *,double *);
  virtual void createcsrki(int *,int *,double *);
  virtual void createcsrcreal(int *,int *,double *);
  virtual void createcsrcimag(int *,int *,double *);
  virtual void createcsrmr(int *,int *,double *);
  virtual void createcsrmi(int *,int *,double *);
  virtual int  solve();
  virtual int  solvecomplex();
  virtual void getsol(double **);
  virtual int geteigenvalue(double *, double *, double *, int);
  virtual int geteigenvectordispersion(double **,double **,int);
  virtual int getconveigen();

private:
  virtual int assemblecomplexPP();
  virtual int assemblecomplexQQ();  
  virtual int assemblecomplexRR();
    
};

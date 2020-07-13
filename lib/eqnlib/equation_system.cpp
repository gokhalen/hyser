/*****************************************************************
*   Nachiket Gokhale gokhalen@bu.edu
*   March 3 2005. 
*   This class, equation_system, defines a consistent interface
*   to use different solvers. Currently it uses PETSc/SlepC (2012)Only
*****************************************************************/

#ifndef HAVE_EQUATION_SYSTEM_HPP
#include "equation_system.hpp"
#endif



equation_system::equation_system(){}
equation_system::equation_system(int num_equations_arg, int which_field_arg, int inewton,int iload,int printflag){}
equation_system::~equation_system(){}


void equation_system::createcsrk(int *row_ptr,int *col_ptr, double *value_arr){}
void equation_system::createcsrkreal(int *row_ptr,int *col_ptr, double *value_arr){}
void equation_system::createcsrkimag(int *row_ptr,int *col_ptr, double *value_arr){}

void equation_system::createcsrm(int *row_ptr,int *col_ptr, double *value_arr){}
void equation_system::createcsrmreal(int *row_ptr,int *col_ptr, double *value_arr){}
void equation_system::createcsrmimag(int *row_ptr,int *col_ptr, double *value_arr){}

void equation_system::createrhs(double **value){}
int  equation_system::solve(){}
int  equation_system::solvecomplex(){}
int  equation_system::getconveigen(){}
void equation_system::getsol(double **){}
int  equation_system::geteigenvalue(double   *kreal, double *kimag,  double *kerror, int ieigen){}
int  equation_system::geteigenvector(double **kreal, double **kimag, int ieigen){}
int  equation_system::geteigenvectordispersion(double **kreal, double **kimag, int ieigen){}

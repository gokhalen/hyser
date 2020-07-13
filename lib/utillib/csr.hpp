/*==========================================================================
  Header file for a matrix in CSR format
  Nachiket Gokhale   gokhalen@bu.edu
============================================================================*/
#ifndef HAVE_CSR_HPP
#define HAVE_CSR_HPP
#endif

#ifndef HAVE_UTILS_HPP
#include "utils.hpp"
#endif



class MATCSR{

 public:
  double *value_arr;
  int *row_ptr;
  int *col_ptr;
  int neqn;
  int maxperrow;

  // size = no of rows, maxperrow max entries per row
public:
  MATCSR();
  MATCSR(int nn, int maxperrow);
  MATCSR(int nn, int band, int *row_ptr, int *col_ptr);
public:
  // adds value vval to location column and row. if zero then creates location.
  void add(int row,int col, double vval);
  void addfill(int row,int col,double vval);

public:
  // checks if there is an entry at row, col
  int isthere(int row,int col);

public:
  void viewcsr(char *);

public:
  // zeroes *value_arr;
  void zeromat();
  
 
public:  
  ~MATCSR();
  
};

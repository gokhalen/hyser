 /*======================================================================
  Nachiket Gokhale gokhalen@bu.edu
  
  // Feb 22. Need to Document This Much Better. This is the first stab at it.

  // Initialization: MATCSR(nrows,maxperrow)
  //                 nrow == number of rows in SQUARE matrix. 
  //                 maxperrow == 
  //
  
=======================================================================*/
#ifndef HAVE_CSR_HPP
#include "csr.hpp"
#endif

#include <iostream>
#include <fstream>
#include <cassert>
using namespace std;

MATCSR::MATCSR(){
}

MATCSR::MATCSR(int nrows, int maxperrow_arg){
  // uses zero based indices 

  assert(nrows>0);assert(maxperrow_arg>0);
  

  // allocate and zero storage
  neqn = nrows;
  this->maxperrow = maxperrow_arg;

  value_arr   = allocvec(neqn*maxperrow,value_arr);zerovec(neqn*maxperrow,value_arr);
  row_ptr = allocvec(neqn+1,row_ptr); zerovec(neqn+1,row_ptr);
  col_ptr = allocvec(neqn*maxperrow+1,col_ptr);zerovec(neqn*maxperrow+1,col_ptr);

  for(int ipos = 0;ipos<neqn*maxperrow+1;ipos++){
    col_ptr[ipos]=-1;
  }

  // set neqn+1 th entry to point to the last entry in col_ptr
  row_ptr[neqn] = 0;
  
}
//======================================================================
MATCSR::MATCSR(int nrows, int maxperrow_arg,int *row_ptr1, int *col_ptr1){
  // preferred method of rows; fill in is given
  // much faster than other approach


  neqn = nrows;
  this->maxperrow = maxperrow_arg;

  value_arr   = allocvec(neqn*maxperrow,value_arr);zerovec(neqn*maxperrow,value_arr);
  // pointer arithmetic: row_ptr gets the address of row_ptr1
  // now need to allocate: externally allocated
  row_ptr=row_ptr1;
  col_ptr=col_ptr1;
}

//======================================================================
void MATCSR::addfill(int row,int col,double insval){
  // adds in an element to matrix whose csr format has been set with 
  // the second constructor
  
  // sanity checks
  assert(row>=0);assert(col>=0);
  assert(row<neqn);assert(col<neqn);
  
  // follow row pointer
  
  int rowposincolptr = row_ptr[row];
  int nentriesrow = row_ptr[row+1]-row_ptr[row];


  assert(nentriesrow>0);

  for(int icolpos = rowposincolptr; icolpos<(rowposincolptr+nentriesrow);icolpos++){
    if(col_ptr[icolpos] == col){
      value_arr[icolpos] += insval;
      break;
    }
  }
  
}


//=====================================================================
void MATCSR::add(int row, int col, double insval){
  
  // uses zero based indexing. user must also input in zero based indices
  
  // make sure the entry values are reasonable
  assert(row>-1);assert(col>-1); 
  assert(row<neqn);assert(col<neqn);
  // get the location (in zero based indices) of the start of rowth row
  // in col_ptr and val_arr
  int collocation = row_ptr[row]; assert(collocation>-1);
  
  // find number of entries in that row

  int nentriesrow = row_ptr[row+1]-row_ptr[row];


  // a flag to check if the row, col position 
  // already exists in the matrix
  int colisthereflag = -2;
  int replacecol;
  for(int icol = collocation; icol< (collocation + nentriesrow);icol++){

    if(col_ptr[icol] == col){
      replacecol = icol;
      colisthereflag = 1; break;
      
    }
    if(col_ptr[icol] == -1){
      replacecol = icol;
      colisthereflag = 0; break;
      
    }
    
  }
  
  if(colisthereflag == 1){
    // when col is there just add value
    value_arr[replacecol] += insval;
    
  }
  else if (colisthereflag == 0){
    // case when col_ptr == -1
    value_arr[replacecol] = insval;
    col_ptr[replacecol]   = col;
    
    // and increment row pointers
    for(int irow = row+1; irow<neqn+1; irow++){
      row_ptr[irow]++;
    }
  }
  else{
    int inscolpos = row_ptr[row+1]; //the new position minus one, not for zero index
    
  
    // first, shift pointers for row

    for(int irow = row+1; irow<neqn+1; irow++){
      row_ptr[irow]++;
    }
   
    // then, shift pointers for col and value_arr 
    for(int shiftcolptr = row_ptr[neqn]; shiftcolptr >=row_ptr[row+1];shiftcolptr--){
      value_arr[shiftcolptr] = value_arr[shiftcolptr-1];
      col_ptr[shiftcolptr] = col_ptr[shiftcolptr-1];
    }

  
    //shifting is done, insert value
 
    
    col_ptr[inscolpos] = col;
    value_arr[inscolpos]   = insval;
    
    for(int icol = inscolpos;icol>collocation;icol--){
      if(col_ptr[icol]<col_ptr[icol-1]){
	int  itemp = col_ptr[icol];
	double dtemp    = value_arr[icol];

	col_ptr[icol] = col_ptr[icol-1]; 
	col_ptr[icol-1] = itemp;

	value_arr[icol]   = value_arr[icol-1];
	value_arr[icol-1] = dtemp;
      }
      else{
	break;
      }
    }

  }

}

//======================================================================
int MATCSR::isthere(int row,int col){
  return -1;
}

//=====================================================================
void MATCSR::viewcsr(char *filename){
  // print the matrix to a file
  ofstream outfilename(filename);
  outfilename.precision(12);



  //for (int irow = 0; irow< neqn+1; irow++){
  //  outfilename<<irow <<" "<<row_ptr[irow]<<endl;
  // }

  // cout<<__FILE__<<" "<<__LINE__<<" printing col_ptr"<<endl;

  for(int icolptr=0;icolptr<=row_ptr[neqn];icolptr++){
    outfilename<<icolptr <<" colptr "<<col_ptr[icolptr]<<" value "<<value_arr[icolptr]<<endl;
  }
  outfilename.close();
}
//=======================================================================
void MATCSR::zeromat(){
  zerovec(neqn*maxperrow,value_arr);
}


//=======================================================================
MATCSR::~MATCSR(){
  // we must delete value_arr and possible others (row,ptr) depending on which constructor was called. For now we assume that the second constructor (int, int, int *, int *)was called.
  deallocvec(value_arr);
} 

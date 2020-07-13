/*=================================================================
  Nachiket Gokhale gokhalen@bu.edu
  memutils.cpp Feb 8 2005
  makes it very easy to handle type independent allocation

  NOTE: It seems like implementations of template functions must be accessible as source where it is used. 
 
===================================================================*/

#ifndef HAVE_MEMUTILS_HPP
#define HAVE_MEMUTILS_HPP
#endif

#include <iostream>
#include <new>

using namespace std;

//type independent memory allocation routine for a vector
template <class X> X* allocvec(long int length, X *vec);

//type independent memory allocation routine for a matrix
template <class X> X** allocmat(long int rows, long int cols, X **mat);

//type independent memory allocation routine for a ***
template <class X> X*** alloc3D(long int rows, long int cols, long int depth, X ***tens);

//type independent memory allocation routine for a ****
template <class X> X**** alloc4D(int first,int sec,int thr, int fourth, X**** tens);

// implementations ==========================================================
// allocvec
template<class X> X* allocvec(long int length, X *vec){
  
  try{
    vec = (X *)new X[length];
  }
  catch(bad_alloc xa){
    cout <<"allocation failure in template vector allocation function"; 
  };
  return vec;
}

//============================================================================
//allocmat
template <class X> X** allocmat(long int rows,long int cols,X **mat){
  
  try{
    mat = (X**)new X*[rows];
  }catch(bad_alloc xa){
    cout<<"first allocation error allocmat";
  };
  
  for(long int irow = 0;irow<rows;irow++){
    try{
      mat[irow] = (X *)new X[cols];
    }catch(bad_alloc xa){
      cout<<"second allocation error allocmat";
    };
  }
  return mat;
}


// third method: allocate a 3D tensor (three indexed object)
template <class X> X*** alloc3D(long int rows,long int cols,long int depth,X ***tens){
  // allocate the number of rows
  try{
    tens = (X ***)new X**[rows];
  }catch(bad_alloc xa){
    cout<<"first allocation error alloc3D";
  };
  
  // allocate columns
  for(long int irow = 0;irow<rows;irow++){
    try{
      tens[irow] = (X **)new X*[cols];
    }catch(bad_alloc xa){
      cout<<"second allocation error in template function alloc3D";
    };
  }

  //allocate depth
  for(long int irow  =0;irow<rows;irow++){
    for(long int icol=0;icol<cols;icol++){
      try{
	tens[irow][icol] = (X *) new X[depth];
      } catch(bad_alloc xa){
	cout<<"third (depth) allocation error in template class alloc3D";
      };
    }
  }

 return tens;
}
//============================================================================
template <class X> X**** alloc4D(int first,int second,int third, int fourth, X**** tens){
  // allocate first dimension
  // Then allocate the other three using alloc3D
  try{
    tens = (X ****)new X***[first];
  }catch(bad_alloc xa){
    cout<<"first allocation error alloc4D";
  };
  
  // then can we make use of alloc3D ?
  // yes, quite elegant...:-)
  for(long int ifirst = 0;ifirst<first;ifirst++){
    tens[ifirst] = alloc3D(second,third,fourth,tens[ifirst]);
  }

  return tens;
}
//=========================================================================
//deallocation functions
//=========================================================================
template <class X> void deallocvec(X *vec){

  delete [] vec; // release the array vec;
}
//=========================================================================
template <class X> void deallocmat(long int rows,X **mat){

  for(long int irow = 0;irow<rows;irow++){
    delete [] mat[irow];
  }
  
  delete [] mat;
}
//========================================================================
template <class X> void dealloc3D(long int rows, long int cols, X ***tens){

  for(long int irow = 0;irow < rows ;irow ++){
    for(long int icol = 0;icol < cols;icol++){
      delete [] tens[irow][icol];
    }
    delete [] tens[irow];
  }

  delete [] tens;
}
//=========================================================================
template <class X> void dealloc4D(int first,int second, int third,int fourth, X**** tens){
  
  for(int ifirst = 0; ifirst<first ; ifirst++){
    for(int isecond = 0;isecond<second;isecond++){
      for(int ithird = 0;ithird<third;ithird++){
	delete [] tens[ifirst][isecond][ithird];
      }
      delete [] tens[ifirst][isecond];
    }
    delete [] tens [ifirst];
  }
  delete [] tens;
}

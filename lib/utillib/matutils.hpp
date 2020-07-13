/***********************************************
*  Matrix utilities:
*  Nachiket Gokhale gokhalen@bu.edu
*
*
*
*
*************************************************/

#ifndef HAVE_MATUTILS_HPP
#define HAVE_MATUTILS_HPP
#endif

// provides delta(ii,jj)

#ifndef HAVE_MISC_UTILS_HPP
#include "misc_utils.hpp"
#endif



#include <math.h>
#include <cassert>
#include <cstdlib>  // Feb 23, 2009, added cstdlib to provide abort()

using namespace std;


// can zero doubles ints whatever

template <class XX> void zerovec(int length, XX *mat){
  for (int ilength = 0;ilength <length;ilength++){
    mat[ilength] = 0;
  }
}

template <class XX> void copyvec(int length, XX *mat, XX *matcopy){
  for (int ilength = 0;ilength <length;ilength++){
    matcopy[ilength] = mat[ilength];
  }
}


template <class XX> void zeromat(int rows,int cols, XX **mat){
  for (int irow = 0;irow < rows;irow++){
    for(int icol =0;icol <cols;icol++){
      mat[irow][icol] = 0;
    }
  }
}

template <class XX> inline void zero3D(int rows,int cols, int depth, XX ***mat){
  for (int irow = 0;irow < rows;irow++){
    for(int icol =0;icol <cols;icol++){
      for(int idepth=0;idepth< depth;idepth++){
	mat[irow][icol][idepth] = 0;
      }
    }
  }
}

template <class XXX> inline void zero4D(int aa, int bb, int cc, int dd, XXX ****mat){
  
  for(int irow = 0; irow < aa ; irow++){
    zero3D(bb,cc,dd,mat[irow]);
  }

}



//transpose a square matrix non-inplace
template <class XX> void tpose(int nn,XX **matin,XX **matout){
  assert(nn>0);

  for(int ii = 0 ; ii<nn; ii++){
    for(int jj =0; jj<nn; jj++){
      matout[ii][jj] =  matin[jj][ii];
    }
  }

}



// template class transpose to transpose a square matrix in place
template <class XX> void tposeinplace(int nn, XX **mat){
  assert(nn>0);
  XX xtemp;
  for(int ii =0;ii<nn;ii++){
    // notice the jj<ii, other wise double swapping, no transpose
    for(int jj =0;jj<ii;jj++){
      xtemp = mat[ii][jj];
      mat[ii][jj] = mat[jj][ii];
      mat[jj][ii] = xtemp;
    }
  }
}

// template class to calculate F^{T}F i.e. A_{ij} = F_{mi}F_{mj} 

template <class XX> void tposetimesself(int nn, XX **matin, XX **matout){
  assert(nn>0);
  
  zeromat(nn,nn,matout);
  
  for(int irow = 0; irow<nn;irow++){
    for(int icol = 0;icol<nn;icol++){
      
      for(int mm = 0;mm<nn;mm++){
	matout[irow][icol] += matin[mm][irow]*matin[mm][icol];
      }
    }
  }
 

}



// evaluate determinants

template <class XX> XX determinant(int nn, XX **mat){
  XX deter;
  
  if(nn == 2){
    deter = (mat[0][0]*mat[1][1] - mat[1][0]*mat[0][1]);
  }

  if(nn==3){

    // hope this stuff (defining new variables, creating temporaries) 
    // is optimized away by the compiler
    // in any case since the number of multiplications 
    // stays the same, it should not matter
   
    // anyway, in the timings it doesn't appear 

    XX a1 =  mat[0][0]*mat[1][1]*mat[2][2];
    XX a2 =  mat[0][1]*mat[1][2]*mat[2][0];
    XX a3 =  mat[0][2]*mat[1][0]*mat[2][1];
    XX a4 = -mat[0][2]*mat[1][1]*mat[2][0];
    XX a5 = -mat[0][0]*mat[1][2]*mat[2][1];
    XX a6 = -mat[0][1]*mat[1][0]*mat[2][2];
    deter = (a1 + a2 + a3 + a4 + a5 +a6 );
  }

  return deter;
}

// matrix inverse for small (3X3) matrices 
template <class XX> void invertsmall(int nn, XX **mat, XX **inv){
  assert((nn==2)||(nn==3));
  
  if(nn == 2){
    XX dett = mat[0][0]*mat[1][1] - mat[0][1]*mat[1][0];
    if(dett !=0){
      inv[0][0] =  mat[1][1]/dett;
      inv[0][1] = -mat[0][1]/dett;
      inv[1][0] = -mat[1][0]/dett;
      inv[1][1] =  mat[0][0]/dett;
     }
    else{
      cout<<__FILE__<<__LINE__<<"Determinant of matrix is zero!!!!"<<endl;
      abort();
    }
  }
  if(nn == 3){
    XX dett =  determinant(3,mat);

    // throw an exception here
    // from kreyszig pg 354
    
    assert(dett !=0);
    
    inv[0][0]  =  (mat[1][1]*mat[2][2] - mat[2][1]*mat[1][2])/dett;
    inv[0][1]  = -(mat[1][0]*mat[2][2] - mat[2][0]*mat[1][2])/dett;
    inv[0][2]  =  (mat[1][0]*mat[2][1] - mat[2][0]*mat[1][1])/dett;
    
    inv[1][0]  = -(mat[0][1]*mat[2][2] - mat[2][1]*mat[0][2])/dett;
    inv[1][1]  =  (mat[0][0]*mat[2][2] - mat[2][0]*mat[0][2])/dett;
    inv[1][2]  = -(mat[0][0]*mat[2][1] - mat[2][0]*mat[0][1])/dett;
    
  
    inv[2][0]  =  (mat[0][1]*mat[1][2] - mat[1][1]*mat[0][2])/dett;
    inv[2][1]  = -(mat[0][0]*mat[1][2] - mat[1][0]*mat[0][2])/dett;
    inv[2][2]  =  (mat[0][0]*mat[1][1] - mat[1][0]*mat[0][1])/dett;

    
    tposeinplace(3,inv);

  }
   
}




// returns trace of a square matrix
template <class XX> XX tracemat(int nn, XX **mat){
  assert(nn>0);
  XX value = 0;
  for (int ii = 0;ii<nn;ii++){
    value += mat[ii][ii];
  }
  return value;
}



// to view a matrix


template <class XXX> void vecview(int length, XXX *vec){
  assert(length > 0);

  for(int ii = 0; ii<length; ii++){
    cout<<vec[ii]<<endl;
  }
}

template <class XX> void matview(int rows, int cols, XX **mat){
  
  assert(rows>0); assert(cols>0);

  for(int irow = 0; irow < rows ; irow++){
    for(int icol = 0; icol <cols; icol++){
      cout<<mat[irow][icol]<<"  ";
    }
    cout<<endl;
  }


}

// to calculate and return the norm of a vector

template <class XXX> XXX vecnorm(int length, XXX *vec){
  
  assert(length>0);
  
  double norm = 0.0;

  for(int ii = 0; ii<length ; ii++){
    norm += vec[ii]*vec[ii];
  }
  
  return sqrt(norm);
 
}

// create a unit matrix
template <class XXX> void matunit(XXX **inmatrix,int nn){

  for(int ii = 0; ii < nn ; ii++){
    for(int jj = 0; jj< nn; jj++){
      inmatrix[ii][jj] = delta(ii,jj);
    }
  }
  
}


// calculate the inverse of a matrix (use double precision only)
// gauss-jordan elimination without pivoting from Numerical Recipies in C 
// use for small index values only, and for nice matrices only

template <class XXX> void matinverse(int nn, XXX **inmatrix, XXX **outmatrix){

  XXX **copymat;
  copymat = allocmat(nn,nn,copymat);

  for (int ii = 0; ii < nn; ii++){
    for(int jj =0; jj < nn; jj++){
      copymat[ii][jj] = inmatrix[ii][jj];
    }
  }
  
  matunit(outmatrix,nn);
  
  for(int irow=0; irow< nn; irow++){
    for(int icol=0; icol< nn; icol++){

      if(irow == icol){ //pick diagonal element

	double pivot = copymat[irow][icol];

	if(pivot == 0){
	  cout<<endl;
	  matview(nn,nn,inmatrix);
	  cout<<endl;
	}
	
	assert(pivot !=0);
	

	// divide the irowth row by the pivot for both the copymat and outmatrix

	for(int kcol=0; kcol< nn; kcol++){
	  copymat[irow][kcol] = copymat[irow][kcol]/pivot;
	  outmatrix[irow][kcol] = outmatrix[irow][kcol]/pivot;
	}

	// loop through all the rows (except current pivot one) and make them zero
	for(int krow=0; krow<nn; krow++){
	  double multiplier = copymat[krow][icol];
	  if(krow != irow){
	    for(int kcol = 0; kcol< nn; kcol++){
	      copymat[krow][kcol] = copymat[krow][kcol]  - multiplier*copymat[irow][kcol];
	      outmatrix[krow][kcol] = outmatrix[krow][kcol] -multiplier*outmatrix[irow][kcol];
	    }
	  }
	}

      }
    }
  }
  
  deallocmat(nn,copymat);
}

// matmultiplybyscalar -- the input matrix is modified
template <class XXX> void matmultbyscalar(int rows, int cols, XXX factor,XXX **matrix){

  for(int irow = 0; irow < rows; irow++){
    for(int icol = 0; icol < cols; icol++){
      matrix[irow][icol] = matrix[irow][icol]*factor;
    }
  }

}



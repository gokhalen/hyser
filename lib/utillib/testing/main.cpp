#include <iostream>
#include <fstream>

#include "../utils.hpp"

int main(){
  
  double **mat;
  double **inv;
  double **correctinverse;

  mat = allocmat(3,3,mat);
  inv = allocmat(3,3,mat);
  correctinverse = allocmat(3,3,correctinverse);

  correctinverse[0][0] = -0.7;
  correctinverse[0][1] =  0.2;
  correctinverse[0][2] =  0.3;
  correctinverse[1][0] = -1.3;
  correctinverse[1][1] = -0.2;
  correctinverse[1][2] =  0.7;
  correctinverse[2][0] =  0.8;
  correctinverse[2][1] =  0.2;
  correctinverse[2][2] = -0.2;
 
  double dett;

  // from kreyszig  pg 354
  
  mat[0][0] = -1;
  mat[0][1] =  1;
  mat[0][2] =  2;

  mat[1][0] =  3;
  mat[1][1] = -1;
  mat[1][2] =  1;

  mat[2][0] = -1;
  mat[2][1] =  3;
  mat[2][2] =  4;
  
  //mat[0][0] =  1;
  //mat[0][1] =  2;
  //mat[0][2] =  3;
  //mat[1][0] = 11;
  //mat[1][1] = 17;
  //mat[1][2] = 21;
  //mat[2][0] = 21;
  //mat[2][1] =  3;
  //mat[2][2] =  4;

  dett = determinant(3,mat);
  cout<<"Inverting using inversmall()"<<endl;
  cout<<"Determinant == "<<dett<<endl;
  invertsmall(3,mat,inv);
  matview(3,3,inv);
  
  cout<<"Inverting using matinverse()"<<endl;
  matinverse(3,mat,inv);
  matview(3,3,inv);

  cout<<"Original matrix from matinverse()"<<endl;
  matview(3,3,mat);


  
  cout<<endl<<"Correct inverse is "<<endl;
  matview(3,3,correctinverse);

  deallocmat(3,mat);
  deallocmat(3,inv);
  deallocmat(3,correctinverse);

  return 0;
}

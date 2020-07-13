#include <iostream>
#include "../utils.hpp"
using namespace std;

#define SIZE 10

int main(){
  
  double *ptr1;
  double **ptr2;
  double ***ptr3;
  double ****ptr4;
  double CIJKL[SIZE][SIZE][SIZE][SIZE];
  double *cptr;

  cptr = &CIJKL[0][0][0][0];
  while(1){
    
    ptr1 = allocvec(SIZE,ptr1);
    ptr2 = allocmat(SIZE,SIZE,ptr2);
    ptr3 = alloc3D(SIZE,SIZE,SIZE,ptr3);
    ptr4 = alloc4D(SIZE,SIZE,SIZE,SIZE,ptr4);
    

    deallocvec(ptr1);
    deallocmat(SIZE,ptr2);
    dealloc3D(SIZE,SIZE,ptr3);
    dealloc4D(SIZE,SIZE,SIZE,SIZE,ptr4);
  }
  return 0;
}

#include <iostream>
#include "../utils.hpp"
using namespace std;

int main(){
  double ****cijkl;

  cijkl = alloc4D(2,2,2,2,cijkl);

  int counter = 0;
  for(int ii = 0;ii<2;ii++){
     for(int jj = 0;jj<2;jj++){
       for(int kk = 0;kk<2;kk++){
	  for(int ll = 0;ll<2;ll++){
	    cijkl[ii][jj][kk][ll] = counter;
	    counter++;
	    cout<<"CIJKL == "<<cijkl[ii][jj][kk][ll]<<endl;
	  }
       }
     }
  }

  dealloc4D(2,2,2,2,cijkl);
  
  return 0;
}

#include <iostream>
#include <../shapestrct.cpp>
using namespace std;

class teststruct{
public:
  int *ptr3;
  teststruct(){

    ptr3 = (int *) new (int) [3];
    
  }


  ~teststruct(){
    delete [] ptr3;
  }

};


int main(){
  

  while(1){
  
   
    delete tstrct;
    
   }

  return 0;
}

#ifndef HAVE_LINELAS2D9_HPP
#define HAVE_LINELAS2D9_HPP
#endif


#ifndef HAVE_ELEMLIB_HPP
#include "elemlib.hpp"
#endif


class linelas2d9:public elem_base{

public:
  
  linelas2d9(int*,int,long int);     // constructor
  ~linelas2d9();    // destructor 

  void elemstiff(int);
  //void calcfunc();
  void elemrhs();
  void calcgradient();
  void getinteg();
  void getshape();
  void printeltype();
  // calculate initial strain body force
  // homogenization functions
  void eleminitstrn();
  void homogenizehassanihinton();
  //void memalloc();
  //void memdealloc();
  
  //void setcoord(double, int, int);
  //void setprop(double, int, int);
  //void setideqn(int,int,long int);
  //void setidnode(int,long int);
  //void setisbc(int,int,int);
  //void setdirdata(double, int, int);
  //void settracdata(double, int, int);
  
};



#ifndef HAVE_NONLINELASTRAC2D_HPP
#define HAVE_NONLINELASTRAC2D_HPP
#endif


#ifndef HAVE_ELEMLIB_HPP
#include "elemlib.hpp"
#endif


class nonlinelastrac2d:public elem_base{

public:
  
  nonlinelastrac2d(int*,int,int,long int);     // constructor
  ~nonlinelastrac2d();    // destructor 

  void elemstiff(int);
  void calcfunc();
  void calcgradient();
  void addreggradient();
  void calcdual();
  void elemrhs();
  void getinteg();
  void getshape();
  void printeltype();
  //void memalloc();
  void memdealloc();
  void calcnormmeas();
  
  //void setcoord(double, int, int);
  //void setprop(double, int, int);
  //void setideqn(int,int,long int);
  //void setidnode(int,long int);
  //void setisbc(int,int,int);
  //void setdirdata(double, int, int);
  //void settracdata(double, int, int);
  
};



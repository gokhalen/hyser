#ifndef HAVE_NONLINELAS1D_HPP
#define HAVE_NONLINELAS1D_HPP
#endif


#ifndef HAVE_ELEMLIB_HPP
#include "elemlib.hpp"
#endif


class nonlinelas1d:public elem_base{

  // constructor will allocate set data of shp strct

public:
  
  nonlinelas1d(int*,long int);     // constructor
  ~nonlinelas1d();    // destructor 

  void elemstiff(int);
  void calcfunc();
  void elemrhs();
  void calcgradient();
  void getinteg();
  void getshape();
  void printeltype();
  void memalloc();
  void memdealloc();
  void setcoord(double, int, int);
  void setprop(double, int, int);
  void setideqn(int,int,long int);
  void setidnode(int,long int);
  void setisbc(int,int,int);
  void setdirdata(double, int, int);
  void settracdata(double, int, int);
  void setmeasdata(double,int,int);
  void setdualfield(double,int,int);
  void calcdual();
};

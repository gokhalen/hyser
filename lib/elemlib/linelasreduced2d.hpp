#ifndef HAVE_LINELASREDUCED2D_HPP
#define HAVE_LINELASREDUCED2D_HPP
#endif


#ifndef HAVE_ELEMLIB_HPP
#include "elemlib.hpp"
#endif


class linelasreduced2d:public elem_base{

public:
  
  linelasreduced2d(int*,int,long int);     // constructor
  ~linelasreduced2d();    // destructor 

  void elemstiff(int);
  void elemcmass(int);
  //void calcfunc();
  void elemrhs();
  void calcgradient();
  void getinteg();
  void getshape();
  void printeltype();
  //void memalloc();
  void eleminitstrn();
  void homogenizehassanihinton();
  //void averagedensity();
  //void memdealloc();
  
  //void setcoord(double, int, int);
  //void setprop(double, int, int);
  //void setideqn(int,int,long int);
  //void setidnode(int,long int);
  //void setisbc(int,int,int);
  //void setdirdata(double, int, int);
  //void settracdata(double, int, int);
  
};



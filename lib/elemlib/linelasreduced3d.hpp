#ifndef HAVE_LINELASREDUCED3D_HPP
#define HAVE_LINELASREDUCED3D_HPP
#endif

#ifndef HAVE_ELEMLIB_HPP
#include "elemlib.hpp"
#endif

class linelasreduced3d:public elem_base{

public:
  
  linelasreduced3d(int*,int,long int);     // constructor
  ~linelasreduced3d();    // destructor 

  void elemstiff(int);
  //void calcfunc();
  void elemrhs();
  void calcgradient();
  void getinteg();
  void getshape();
  void printeltype();

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



#ifndef HAVE_NONLINELAS2D9_HPP
#define HAVE_NONLINELAS2D9_HPP
#endif


#ifndef HAVE_ELEMLIB_HPP
#include "elemlib.hpp"
#endif


class nonlinelas2d9:public elem_base{

  // constructor will allocate set data of shp strct

public:
  
  nonlinelas2d9(int*,int,int,long int);     // constructor
  ~nonlinelas2d9();    // destructor 

  void elemstiff(int);
  void elemrhs();
  void calcgradient();
  void getinteg();
  void getshape();
  void printeltype();
  //void memalloc();
  //void memdealloc();
  //void calcdual();

  
  // I'm inheiriting all this from elem_base 
  // every element uses the same form of the function
  
  //void setcoord(double, int, int);
  //void setprop(double, int, int);
  //void setideqn(int,int,long int);
  //void setidnode(int,long int);
  //void setisbc(int,int,int);
  //void setdirdata(double, int, int);
  //void settracdata(double, int, int);
  //void setmeasdata(double, int, int);
  //void calcfunc();
  //void setdualfield(double ddata, int inode, int idofn);
  //void calculate_Tmeas_Tprim();

  // calcnormmeas is not defined here. it is used from elem_base.cpp



};

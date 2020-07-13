#ifndef HAVE_PIEZOREDUCED2D_HPP
#define HAVE_PIEZOREDUCED2D_HPP
#endif


#ifndef HAVE_ELEMLIB_HPP
#include "elemlib.hpp"
#endif


class piezoreduced2d:public elem_base{

public:
  
  piezoreduced2d(int*,int,long int);     // constructor
  ~piezoreduced2d();    // destructor 

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
  void computeaveprops();

  //void averagedensity();
  //void memdealloc();
  
  //void setcoord(double, int, int);
  //void setprop(double, int, int);
  //void setideqn(int,int,long int);
  //void setidnode(int,long int);
  //void setisbc(int,int,int);
  //void setdirdata(double, int, int);
  //void settracdata(double, int, int);
  private:
  double uscal,vscal,uvscal,uscal2,vscal2;
  double getpz(int,int,int,double,double,double);
  double getelas(double,double, int,int,int,int);
  double getkk(int,int,double,double);
  void getkl(int *,int *);
};



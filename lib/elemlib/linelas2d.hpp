#ifndef HAVE_LINELAS2D_HPP
#define HAVE_LINELAS2D_HPP
#endif

#ifndef HAVE_ELEMLIB_HPP
#include "elemlib.hpp"
#endif

class linelas2d:public elem_base{
public:
  
  linelas2d(int*,int,long int);     // constructor
  ~linelas2d();    // destructor 

  void elemstiff(int);
  void elemrhs();
  void calcgradient();
  void getinteg();
  void getshape();
  void printeltype();
  void eleminitstrn(); // calculate initial strain body force
  // homogenization functions
  void homogenizehassanihinton();
  
};



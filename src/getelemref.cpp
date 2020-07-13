/*==================================================================
 *  Nachiket Gokhale gokhalen@bu.edu  getelemref.cpp Feb 10 2005
 *  
 *  Gets a reference to a finite element of a particular type
 *  Allocates memory
 *
 *================================================================== */

#ifndef HAVE_MESH_HPP
#include "../lib/meshlib/mesh.hpp"
#endif  

#ifndef HAVE_ELEMLIB_HPP
#include "../lib/elemlib/elemlib.hpp"
#endif

#include <iostream>
#include <typeinfo>
using namespace std;

elem_base* getelemref(mesh *inmesh,int locelem, elem_base *felement){
  
    // set  the integration points
  int ninte[3];
  ninte[0] = ninte[1] = ninte[2] = inmesh->ninte;
  
  switch(inmesh->lmref[locelem][0]){

  case LINELAS1D:
    felement = (elem_base *)(new linelas1d(&ninte[0],locelem));
    break;
    
  case LINELAS2D:
    felement = (elem_base *)(new linelas2d(&ninte[0],inmesh->mregparm,locelem));
    break;

  case LINELAS2D9:
    felement = (elem_base *)(new linelas2d9(&ninte[0],inmesh->mregparm,locelem));
    break; 

  case PIEZOREDUCED2D:
    felement = (elem_base *)(new piezoreduced2d(&ninte[0],inmesh->mregparm,locelem));
    break;
    
  case LINELASREDUCED2D:
    felement = (elem_base *)(new linelasreduced2d(&ninte[0],inmesh->mregparm,locelem));
    break;
    
  case LINELAS3D:
    //cout<<endl<<"Calling Linelas3d"<<endl;
    felement = (elem_base *)(new linelas3d(&ninte[0],inmesh->mregparm,locelem));
    break;	

  case LINELASREDUCED3D:
    //cout<<endl<<"Calling Linelas3d"<<endl;
    felement = (elem_base *)(new linelasreduced3d(&ninte[0],inmesh->mregparm,locelem));
    break;			     

  case LINELASTRAC2D:
    felement = (elem_base *)(new linelastrac2d(&ninte[0],inmesh->nprop,inmesh->mregparm,locelem));
    break;

  case NONLINELAS1D:
    felement = (elem_base *)(new nonlinelas1d(&ninte[0],locelem));
    //return (elem_base* )(new nonlinelas1d(&ninte[0],locelem));
    break;

  case NONLINELAS2D:
    felement = (elem_base *)(new nonlinelas2d(&ninte[0],inmesh->nprop,inmesh->mregparm,locelem));
    break;

  case NONLINELAS2D9:
    felement = (elem_base *)(new nonlinelas2d9(&ninte[0],inmesh->nprop,inmesh->mregparm,locelem));
    break;
 
      
  case NONLINELASTRAC2D:
    felement = (elem_base *)(new nonlinelastrac2d(&ninte[0],inmesh->nprop,inmesh->mregparm,locelem));
    break;
    
  case NONLINHEAT2D:
    felement = (elem_base *)(new nonlinheat2d(&ninte[0],locelem));
    break;

  case NONLINELAS3D:
    felement = (elem_base *)(new nonlinelas3d(&ninte[0],inmesh->nprop,inmesh->mregparm,locelem));
    break;
    
  default:
    cout<<"getelemref.cpp"<<endl<<" Element type "<<inmesh->lmref[locelem][0]<<" not found"<<endl;
    abort();
  } 

    
    
  
  return felement; 
}

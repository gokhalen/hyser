/********************************************************************
 *  Nachiket Gokhale gokhalen@bu.edu:
 *  Called in main.c. Loops over elements and puts elemental stiffness
 *  and right-hand side contributions in the global matrix.
 *  It makes use of the element library and integration library
 **********************************************************************/
#ifndef HAVE_MESH_HPP
#include "../lib/meshlib/mesh.hpp"
#endif

#ifndef HAVE_FUNCDECLS_H
#include "../headers/funcdecls.h"
#endif

#ifndef HAVE_ELEMLIB_HPP
#include "../lib/elemlib/elemlib.hpp"
#endif

#ifndef HAVE_GLOBALS_H
#include "../headers/globals.h"
#endif


int buildtan(mesh *inmesh){
  
  int which_stiffness = BUILD_REAL_STIFFNESS;

  if (inmesh->issolvedispersion()){
    if ( inmesh->solvetype == SOLVE_TYPE_DISPERSIONW2K ){
      which_stiffness = BUILD_COMPLEX_STIFFNESS_W2K;
    }
    if ( inmesh->solvetype == SOLVE_TYPE_DISPERSIONK2W ){
      which_stiffness = BUILD_COMPLEX_STIFFNESS_K2W;
    }
  }
  
  for(int ielem = 0;ielem < inmesh->num_elem;ielem++){
    elem_base *felement; // a pointer to a base class element
    felement = getelemref(inmesh,ielem,felement); // allocate memory 
    // localize data
    felement = localizer(inmesh,ielem,felement,CURRENT_LOAD_ITERATION); 
    felement->getinteg();                         // get integration points
    felement->getshape();                         // get shape functions
    felement->elemstiff(which_stiffness);         // calc. elem stiffness
    felement->elemcmass(inmesh->which_mass);      // calc. elem mass
    assembletan(inmesh,ielem,felement);          
    delete felement;
  }


 
  return 0;
}

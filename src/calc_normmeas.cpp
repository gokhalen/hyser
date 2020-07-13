/********************************************************************
 * Nachiket Gokhale
 * September 7 2005
 *
 * Calculates the norm of the measured data
 * This quantity can be used to estimate the regularization parameter
 *
 * |T(U^{m} - T(U^{h})| = C\Delta|T(U^{m})|
 *******************************************************************/


#ifndef HAVE_FUNCDECLS_H
#include "../headers/funcdecls.h"
#endif


#ifndef HAVE_ELEMLIB_HPP
#include "../lib/elemlib/elemlib.hpp"
#endif

#ifndef HAVE_MESH_HPP
#include "../lib/meshlib/mesh.hpp"
#endif

#ifndef HAVE_GLOBALS_H
#include "../headers/globals.h"
#endif

int calc_normmeas(mesh *inmesh){

  for(int ifield = 0; ifield<inmesh->nfields; ifield++){
    
    inmesh->which_field = ifield;
    
    for(int ielem = 0; ielem<inmesh->num_elem; ielem++){

      elem_base *felement; // a element base class pointer 
      felement = getelemref(inmesh,ielem,felement);
      
      // localize
      
      felement = localizer(inmesh,ielem,felement,LAST_LOAD_ITERATION);
      felement->getinteg();                         // get integration points
      felement->getshape();                         // get shape functions
      
      felement->calcnormmeas();
      
      inmesh->accumulatenormmeas(felement->normmeasdata);
      delete felement;
    
    }//ielem

    inmesh->writemeasnorm();

  }//ifield

  return 0;
}

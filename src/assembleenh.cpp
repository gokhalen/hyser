/*=========================================================================
  Nachiket Gokhale gokhalen@bu.edu 
  Assembly of a finite element matrix
  ============================================================================*/
#ifndef HAVE_FUNCDECLS_H
#include "../headers/funcdecls.h"
#endif

#ifndef HAVE_PETSC_HEADERS_H
#include "../headers/petsc_headers.h"
#endif

#ifndef HAVE_MESH_HPP
#include "../lib/meshlib/mesh.hpp"
#endif

int assembleenh(mesh* inmesh, int ielem,elem_base* felem){
  
  for(int ii = 0; ii < 4; ii++){
    inmesh->enhvec[ielem][ii][inmesh->which_field] += felem->enhinc[ii];
  }
  
  return 0;
}


 


/*****************************************************************************
 * Nachiket Gokhale gokhalen@bu.edu
 * April 04 2004
 * 
 * This file calculates the forcing for the dual and stores it 
 * in a mesh object.
 ****************************************************************************/

#ifndef HAVE_FUNCDECLS_H
#include "../headers/funcdecls.h"
#endif

#ifndef HAVE_MESH_HPP
#include "../lib/meshlib/mesh.hpp"
#endif

#ifndef HAVE_ELEM_BASE_HPP
#include "../lib/elemlib/elem_base.hpp"
#endif

#ifndef HAVE_GLOBALS_H
#include "../headers/globals.h"
#endif

int builddualrhs(mesh *inmesh){
  // zero once outside loop.

  int which_stiffness = BUILD_REAL_STIFFNESS;
  if ( inmesh->solvetype == SOLVE_TYPE_DISPERSIONK2W ){
    which_stiffness = BUILD_COMPLEX_STIFFNESS_K2W;
  }

  inmesh->zerorhs();
  for(int ielem = 0;ielem<inmesh->num_elem;ielem++){
    
    elem_base *felement;                             // a base class pointer, deleted at the end
    felement = getelemref(inmesh,ielem,felement);    // allocate memory 
    
    // this is a hack. we localize the boundary data corresp. final soln.
    felement = localizer(inmesh,ielem,felement,LAST_LOAD_ITERATION);
    
    felement->getinteg();                         // get integration points
    felement->getshape();                         // get shape functions
    felement->elemstiff(which_stiffness);   // need matrices for static condensation
    felement->calcdual();                         // calculate the dual rhs
    
    // we have the rhs for the dual, assemble it.
    // use the same array to store the rhs vector for dual and primal
    
  
    int ieqn =0;
    for(int inode =0;inode<felement->shpstrct->numnodes;inode++){
      for(int idofn=0;idofn<felement->dofnatnode[inode];idofn++){
	long int ipos = felement->ideqn[inode][idofn] -1;
	
	if(ipos>=0){
	  inmesh->rhsvec[ipos][inmesh->which_field] += felement->erhs[ieqn];
	}
	
	ieqn++;
      }
    }
    delete felement;
  }

  inmesh->writedualrhs();
  return 0;
}

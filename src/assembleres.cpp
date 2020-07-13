/*=========================================================================
  Nachiket Gokhale gokhalen@bu.edu 
  Assembly of a finite element matrix
  ============================================================================*/
#ifndef HAVE_FUNCDECLS_H
#include "../headers/funcdecls.h"
#endif

#ifndef HAVE_MESH_HPP
#include "../lib/meshlib/mesh.hpp"
#endif


#ifndef HAVE_GLOBALS_H
#include "../headers/globals.h"
#endif

int assembleres(mesh* inmesh,long int ielem,elem_base* felem){
  
  int ieqn =0;
  for(int inode = 0;inode<felem->shpstrct->numnodes;inode++){
    for(int idofn=0;idofn<felem->dofnatnode[inode];idofn++){
	long int ipos = felem->ideqn[inode][idofn] -1;

	if(ipos>=0){
	  inmesh->rhsvec[ipos][inmesh->which_field] += felem->erhs[ieqn];
	}
	// nodal forces are added directly into the global vector
	ieqn++;
    }
  }
  
 
  
  
  return 0;
}


 


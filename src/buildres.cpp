/*Nachiket Gokhale gokhalen@bu.edu:
Called in main.c. Loops over elements and puts elemental stiffness
and right-hand side contributions in the global matrix.
It makes use of the element library and integration library
*/
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


int buildres(mesh *inmesh){

  int which_stiffness = BUILD_REAL_STIFFNESS;
  if ( inmesh->solvetype == SOLVE_TYPE_DISPERSIONW2K ){
    which_stiffness = BUILD_COMPLEX_STIFFNESS_W2K;
  }
  
  if ( inmesh->solvetype == SOLVE_TYPE_DISPERSIONK2W ){
    which_stiffness = BUILD_COMPLEX_STIFFNESS_K2W;
  }



  int which_field = inmesh->which_field;
  for(int ielem = 0;ielem < inmesh->num_elem;ielem++){
    //cout<<__FILE__<<__LINE__<<"ielem=="<<ielem<<endl;
    elem_base *felement; // a pointer to a base class element
   
    bool integrateflag;
    // by default the integration flag is set to true
    integrateflag  = true;

    // this allocation is outside
    felement = getelemref(inmesh,ielem,felement); // allocate memory 
    // localize data
    felement = localizer(inmesh,ielem,felement,CURRENT_LOAD_ITERATION); 
    // for nonlinear we calculate the rhs integral over the whole domain
    // for linear we calculate the rhs integral only over dirichlet bound.
    
    // skip this section for linear problems
    // the integrateflag flag is set only for the boundary elements

    if(inmesh->linearflag == 1){
      integrateflag = false;
      for(int inode = 0; inode < felement->shpstrct->numnodes; inode++){
	for(int idofn = 0; idofn <felement->dofnatnode[inode];idofn++){           
	  int globalnode=inmesh->ien[ielem][inode];

	  if (( inmesh->isdirichlet(globalnode-1,idofn,which_field) ) || ( inmesh->istraction(globalnode-1,idofn,which_field)) ){
	    integrateflag = true;
	  }

	}
      }

      // bodyforce forcing
      {
	felement->getinteg();
	felement->getshape();
	felement->eleminitstrn();
	assembleres(inmesh,ielem,felement);
      }
    }
    
    //cout<<__FILE__<<__LINE__<<"integrateflag== "<<integrateflag<<endl;
    
    if(integrateflag){
      //cout<<__FILE__<<__LINE__<<"RHS Called for element number "<<ielem<<endl;
      felement->getinteg();  // get integration points
      felement->getshape();  // get shape functions
      felement->elemstiff(which_stiffness); // build element stiffness --not necess for nonlin
      felement->elemrhs();   // calc. elem rhs
      assembleres(inmesh,ielem,felement);
    }  
    
   
    delete felement;
  } // closes global element loop

  // assemble nodal forces

  for(int inode = 0; inode < inmesh->num_nodes; inode++){
    for(int idofn = 0; idofn < inmesh->mdofn; idofn++){
      int pointlocation = inmesh->ispointbc[inode][idofn][inmesh->which_field];
      int ipos  = inmesh->id[inode][idofn][inmesh->which_field]-1;
      if(pointlocation > 0){
	inmesh->rhsvec[ipos][inmesh->which_field] += inmesh->pforc[pointlocation-1][inmesh->iloadits][inmesh->which_field];
      }
	
    }
  }
  return 0;
}

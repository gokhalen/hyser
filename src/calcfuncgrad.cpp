/***************************************************************************
 * Nachiket Gokhale gokhalen@bu.edu
 * April 3 2005
 * 
 * Calculates functional and gradient vector.
 * This assumes that the dual, primal are calculated.
 ***************************************************************************/

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

int calcfuncgrad(mesh *inmesh){
  // this procedure is very similar to the buildtanres
  // the idea is to get a element pointer (felement) 
  // localize data, initialize shape functions etc 
  // and finally call the functional routine in each element

  // functional and gradient are  zeroed externally

  int which_stiffness = BUILD_REAL_STIFFNESS;

  if ( inmesh->solvetype == SOLVE_TYPE_DISPERSIONW2K ){
    which_stiffness = BUILD_COMPLEX_STIFFNESS_W2K;
  }
  
  if ( inmesh->solvetype == SOLVE_TYPE_DISPERSIONK2W ){
    which_stiffness = BUILD_COMPLEX_STIFFNESS_K2W;
  }

  for(int ielem = 0;ielem<inmesh->num_elem;ielem++){
    //cout<<__FILE__<<" "<<__LINE__<<" ielem== "<<ielem<<endl;
    elem_base *felement;                             // a base class pointer
    felement = getelemref(inmesh,ielem,felement);    // allocate memory 
    
    // this is a hack. we localize the boundary data corresp. final soln.
    felement = localizer(inmesh,ielem,felement,LAST_LOAD_ITERATION);    //cout<<endl<<__FILE__<<__LINE__<<endl;
    felement->getinteg();   // get integration points     //cout<<endl<<__FILE__<<__LINE__<<" get integ"<<endl;
    felement->getshape();   // get shape functions        //cout<<endl<<__FILE__<<__LINE__<<" get shape"<<endl;
    felement->elemstiff(which_stiffness);                        // get the local stiffnesses
    felement->calcfunc();                         // calculate the functional

    // add the functionals

    inmesh->functarray[inmesh->which_field]  += felement->functional; 
    inmesh->functional   += felement->functional;

    inmesh->regfunctarray[inmesh->which_field] += felement->regfunctional;  //cout<<endl<<__FILE__<<__LINE__<<"calc func"<<endl;
  
    felement->calcgradient();      // calculate the gradient  
    felement->addreggradient();    // add gradient of reg func     //cout<<endl<<__FILE__<<__LINE__<<"add reggradient"<<endl;

    // store the gradient into the mesh
    int inodemax=felement->shpstrct->numnodes;
    if(felement->propdiscontinflag == 1){
      inodemax=1; 
    }
    for(int inode=0;inode<inodemax;inode++){
      int globalnode = inmesh->ien[ielem][inode];
      if(felement->propdiscontinflag==1){
	// +1 beause the index is globalnode -1;
	globalnode=ielem+1;
      }
      if ( inmesh->opttype == 0 ){
	for(int iprop=0;iprop<felement->nprop;iprop++){
	  inmesh->gradient[globalnode-1][iprop][inmesh->which_field] += felement->gradient[inode][iprop];
	}
      }
      else if (inmesh->opttype == 1 ) {
	inmesh->gradient[globalnode-1][0][inmesh->which_field] += felement->gradient[inode][0] + (inmesh->lamtomu)*felement->gradient[inode][1];
      }
      else if (inmesh->opttype == 2 ){ 
	inmesh->gradient[globalnode-1][1][inmesh->which_field] += felement->gradient[inode][1] + (inmesh->mutolam)*felement->gradient[inode][0];
      }
    }
    
    // destroy the element
    delete felement;
  }

  // we write functional and gradient values in gradientcallbackfortao.cpp

  cout<<__FILE__<<__LINE__<<"done with calcfuncgrad"<<endl;

  return 0;
}

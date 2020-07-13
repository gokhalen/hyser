// this is  homogenization based on hassani-and-hinton
// the structure of this file is very similar to calcfuncgradient.cpp
// Nachiket Gokhale gokhale@wai.com gokhalen@gmail.com
// March 24 2009

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


int homogenize (mesh* inmesh){
  
  cout << endl<<"Start Homognization Field "<< inmesh->which_field<<endl;
  for (int ielem = 0; ielem < inmesh->num_elem; ielem++){
    elem_base *felement;

    // get an element reference 
    felement = getelemref(inmesh,ielem,felement);
    felement = localizer(inmesh,ielem,felement,LAST_LOAD_ITERATION);

    // get integration points
    felement->getinteg();

    // initialize shape functions 
    felement->getshape();
    
    

    felement->homogenizehassanihinton();

    for (int ii = 0; ii < inmesh->ndime; ii++){
      for (int  jj = 0; jj < inmesh->ndime; jj++){
	for (int kk = 0; kk < inmesh->ndime; kk++){
	  for (int ll = 0; ll < inmesh->ndime; ll++){
	    inmesh->EHIJKL[ii][jj][kk][ll] += felement->BASE_EHIJKL[ii][jj][kk][ll];
	    // inmesh->EHIJKL[ii][jj][kk][ll] +=1;
	  }//ll
	  

	  if ( inmesh->physics == PHYSICS_PIEZO){
	    inmesh->PZHKIJ[kk][ii][jj] += felement->BASE_PZHKIJ[kk][ii][jj];
	  } 
	}//kk

	if (inmesh->physics ==PHYSICS_PIEZO ){
	  inmesh->KHIJ[ii][jj] += felement->BASE_KHIJ[ii][jj];
	}

      }//jj
    }//ii
    
 
    delete felement;
  }

  
  cout << endl<<"End Homognization Field "<< inmesh->which_field<<endl;
  return 0;
}

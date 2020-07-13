 /*Nachiket Gokhale gokhalen@bu.edu:

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


int buildenh(mesh *inmesh){

  int which_stiffness = BUILD_REAL_STIFFNESS;

  if ( inmesh->solvetype == SOLVE_TYPE_DISPERSIONK2W ){
    which_stiffness = BUILD_COMPLEX_STIFFNESS_K2W;
  }
  if ( inmesh->solvetype == SOLVE_TYPE_DISPERSIONW2K ){
    which_stiffness = BUILD_COMPLEX_STIFFNESS_W2K;
  }


  for(int ielem = 0;ielem < inmesh->num_elem;ielem++){
    //cout<<__FILE__<<__LINE__<<"ielem=="<<ielem<<endl;
    
    
    //for(int ii = 0; ii < 1; ii++){
      elem_base *felement; // a pointer to a base class element
      
      bool integrateflag;
      // by default the integration flag is set to true
      integrateflag  = true;

    
      
      if(inmesh->isenhanced(ielem)){
	// this allocation is outside
	felement = getelemref(inmesh,ielem,felement); // allocate memory 
	// localize data
	felement = localizer(inmesh,ielem,felement,ENHANCED); 
		
	felement->getinteg();                         // get integration points
	felement->getshape();  
      

	felement->elemstiff(which_stiffness);
	felement->calcenh();
	assembleenh(inmesh, ielem, felement);

	  
	delete felement;
	
      }
    
      //}
  }//ielem


 
  return 0;
}

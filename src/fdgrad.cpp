/********************************************************************
 * Nachiket Gokhale, gokhalen@bu.edu
 * 3rd April 2005
 * Finite difference gradient calculation.
 * This calculates the gradient with finite differences
 * And checks this against the adjoint calculation.
 *********************************************************************/

#ifndef HAVE_FUNCDECLS_H
#include "../headers/funcdecls.h"
#endif

#ifndef HAVE_MESH_HPP
#include "../lib/meshlib/mesh.hpp"
#endif

int fdgrad(mesh* inmesh){
  
  // loop over the properties
  //for(int inode = 0;inode<inmesh->num_nodes;inode++){
  
  double *functional1, *functional2;
  functional1 = allocvec(inmesh->nfields,functional1);
  functional2 = allocvec(inmesh->nfields,functional2);
  
  // independent properties
  int indprop = 1;
  if ( inmesh->opttype == 0 ){indprop = 2;}
    
  for(int inode = 0;inode<1;inode++){
    for(int iprop = 0;iprop<indprop;iprop++){

      if ( inmesh->opttype==0){
	inmesh->prop[inode][iprop] += inmesh->fdstep;
      }
      else if (inmesh->opttype==1){
	inmesh->prop[inode][0] += inmesh->fdstep;
	inmesh->prop[inode][1]  = inmesh->lamtomu*inmesh->prop[inode][0];
      }
      else if (inmesh->opttype==2){
	inmesh->prop[inode][1] += inmesh->fdstep;
	inmesh->prop[inode][0]  = inmesh->mutolam*inmesh->prop[inode][1];
      }

      for (int ifield = 0; ifield <inmesh->nfields; ifield++){
	inmesh->which_field = ifield;
	inmesh->make_fillin();
	inmesh->init_primal();
	inmesh->writeproperties();
	// need to solve and calculate functional
	calc_forward(inmesh);
	inmesh->delete_stiff();
      }
    
  
      // can move the homogenizezation routine to inmesh
      if ( inmesh->linearflag  == 1) {
	// zero the homogenized tensor -- this has to be done outside the fields loop
	zero4D(inmesh->ndime,inmesh->ndime,inmesh->ndime,inmesh->ndime,inmesh->EHIJKL);
	for (int ifield = 0; ifield < inmesh->nfields ; ifield++){
	  // this sets the field
	  inmesh->which_field = ifield;
	  // then call the homogenization code
	  homogenize(inmesh);
	}
	inmesh->writeehijkl();
      }
  
  
      for (int ifield = 0; ifield <inmesh->nfields; ifield++){
	inmesh->which_field = ifield;
	calcfuncgrad(inmesh);
	functional1[ifield] = inmesh->functional;
	//reset the value of the functional in inmesh - we're getting the regularized functional
	inmesh->functional = 0.0;
	
      }
      
      // step 2

      
      if ( inmesh->opttype==0){
	inmesh->prop[inode][iprop] -= 2*inmesh->fdstep;
      }
      else if (inmesh->opttype==1){
	inmesh->prop[inode][0] -= 2*inmesh->fdstep;
	inmesh->prop[inode][1]  = inmesh->lamtomu*inmesh->prop[inode][0];
      }
      else if (inmesh->opttype==2){
	inmesh->prop[inode][1] -= 2*inmesh->fdstep;
	inmesh->prop[inode][0]  = inmesh->mutolam*inmesh->prop[inode][1];
      }
      
      //inmesh->prop[inode][iprop] -= 2*inmesh->fdstep;
      
      for (int ifield = 0; ifield <inmesh->nfields; ifield++){
	inmesh->which_field = ifield;
	inmesh->make_fillin();
	inmesh->init_primal();
	inmesh->writeproperties();
	calc_forward(inmesh);
	inmesh->delete_stiff();
		
      }
      
      // can move the homogenizezation routine to inmesh
      if ( inmesh->linearflag  == 1) {
	// zero the homogenized tensor -- this has to be done outside the fields loop
	zero4D(inmesh->ndime,inmesh->ndime,inmesh->ndime,inmesh->ndime,inmesh->EHIJKL);
	for (int ifield = 0; ifield < inmesh->nfields ; ifield++){
	  // this sets the field
	  inmesh->which_field = ifield;
	  // then call the homogenization code
	  homogenize(inmesh);
	}
	inmesh->writeehijkl();
      }

       for (int ifield = 0; ifield <inmesh->nfields; ifield++){
	inmesh->which_field = ifield;
	calcfuncgrad(inmesh);
	functional2[ifield] = inmesh->functional;
	//reset the value of the functional in inmesh - we're getting the regularized functional
	inmesh->functional = 0.0;
	
      }
	

       
       
       for (int ifield = 0; ifield <inmesh->nfields; ifield++){
	 inmesh->which_field = ifield;
	 if(inmesh->opttype == 0){
	   inmesh->fdgradgrad[inode][iprop][inmesh->which_field]  = (functional1[ifield]-functional2[ifield])/(2.0*inmesh->fdstep);
	 }
	 else if (inmesh->opttype == 1){
	   inmesh->fdgradgrad[inode][0][inmesh->which_field]  = (functional1[ifield]-functional2[ifield])/(2.0*inmesh->fdstep);
	 }
	 else if (inmesh->opttype == 2){
	   inmesh->fdgradgrad[inode][1][inmesh->which_field]  = (functional1[ifield]-functional2[ifield])/(2.0*inmesh->fdstep);
	   cout<<"in fdgrad.cpp fdgrad=="<<inmesh->fdgradgrad[inode][0][inmesh->which_field]<<endl;
	 }
	 // write fdgradgrad for that field to a file
	 inmesh->writefdgradgrad();
       }

       if ( inmesh->opttype==0){
	inmesh->prop[inode][iprop] += inmesh->fdstep;
       }
       else if (inmesh->opttype==1){
	 inmesh->prop[inode][0] += inmesh->fdstep;
	 inmesh->prop[inode][1]  = inmesh->lamtomu*inmesh->prop[inode][0];
       }
       else if (inmesh->opttype==2){
	 inmesh->prop[inode][1] += inmesh->fdstep;
	 inmesh->prop[inode][0]  = inmesh->mutolam*inmesh->prop[inode][1];
       }

      // reset propreties
      //inmesh->prop[inode][iprop] += inmesh->fdstep;

    } // prop
  }  // inode
  
  return 0;
}

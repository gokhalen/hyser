/************************************************************************************************
 * Nachiket Gokhale gokhalen@bu.edu
 * Implementation of the callback function for TAO.
 * The tao routine(s) will call this function 
 * when they want/need the objective function(al) \pi and the 
 * gradient vector of \pi with respect to material properties.
 *
 * This function will, therefore, call:
 * a) calc_forward
 * b) calc_dual
 * c) calc_gradient
 *
 * NOTE: the void pointer argument field is used to pass a mesh* pointer
 *       the mesh pointer is passed by TAO to gradientcallbackfortao
 ***********************************************************************************************/

#ifndef HAVE_FUNCDECLS_H
#include "../headers/funcdecls.h"
#endif

#ifndef HAVE_TAO_HEADERS_H
#include "../headers/tao_headers.h"
#endif

// heavily based on FormFunctionGradient in rosenbrock1c

PetscErrorCode gradientcallbackfortao(TaoSolver taosol, Vec X, double *f,Vec G,void *ptr){
  
  cout<<__FILE__<<__LINE__<<"entering ...gradientcallbackfortao.cpp"<<endl;
  
  // pointer to state
  PetscScalar *stateptr, *gradptr;
  
  // get a pointer to the data in X & G
  VecGetArray(X,&stateptr); 
  VecGetArray(G,&gradptr); 
  
  // cast void pointer ptr to (mesh *)
  mesh *meshptr;
  meshptr = (mesh *)ptr;

  

  // set the properties in meshptr to those of stateptr

  {
   int ieqn = 0;
   for(int inode = 0; inode<meshptr->num_prop_nodes; inode++){

     if ( meshptr->opttype == 0 ){
       for(int iprop = 0;iprop<meshptr->nprop;iprop++){
	 meshptr->prop[inode][iprop] = stateptr[ieqn];
	 ieqn++;
       }
     }
     else if ( meshptr->opttype == 1 ) {
       meshptr->prop[inode][0] = stateptr[ieqn];
       meshptr->prop[inode][1] = meshptr->lamtomu*meshptr->prop[inode][0];
       ieqn++;
     }
     else if ( meshptr->opttype == 2 ) { 
       meshptr->prop[inode][1] = stateptr[ieqn];
       meshptr->prop[inode][0] = meshptr->mutolam*meshptr->prop[inode][1];
       ieqn++;
     }
   }
  }
  
  // write properties -- not here any more 
  // we call this from the monitor, so that this gets called every iteration
  //meshptr->writeproperties();
  

  // zero functional functarray and gradient
  // the combined gradient arrays are zeroed
  // before the combination is done
  
  meshptr->zerofunctional();
  meshptr->zerogradient();

  // loop over deformation fields
  for (int ifield = 0; ifield < meshptr->nfields; ifield++){
    
    meshptr->which_field = ifield;
    
    // make fillin allocates a stiffness matrix
    meshptr->make_fillin();
    meshptr->init_primal();

    //meshptr->init_primal(1.0/meshptr->loadits);
    

    //the forward problem for this deformation field.
    calc_forward(meshptr); 

    // delete the stiffness matrix
    meshptr->delete_stiff();
    
    cout<<__FILE__<<__LINE__<<" forward problem done "<<endl;
  }


  
  if ( meshptr->linearflag  == 1) {
    // zero the homogenized tensor -- this has to be done outside the fields loop
    zero4D(meshptr->ndime,meshptr->ndime,meshptr->ndime,meshptr->ndime,meshptr->EHIJKL);
    for (int ifield = 0; ifield < meshptr->nfields ; ifield++){

      // this sets the field
      meshptr->which_field = ifield;
      
      // then call the homogenization code
      homogenize(meshptr);
    }

    meshptr->writeehijkl();
  }

  cout<<endl<<"In gradientcallbackfortao.cpp "<<endl;
    cout<<"D11 = "<<meshptr->EHIJKL[0][0][0][0]<<endl;
    cout<<"D12 = "<<meshptr->EHIJKL[0][0][1][1]<<endl;
    cout<<"D22 = "<<meshptr->EHIJKL[1][1][1][1]<<endl;
    cout<<"D66 = "<<meshptr->EHIJKL[1][0][1][0]<<endl;
    cout<<endl;
    //exit(1);
  


  for ( int ifield = 0; ifield < meshptr->nfields; ifield++){
    meshptr->which_field = ifield;
    
    // make fillin allocates a stiffness matrix
    meshptr->make_fillin();
    //meshptr->init_primal();
    
    // IMPORANT: Read [1] at the end of file.
    // the dual problem for this deformation field.
    calc_dual(meshptr);

    // delete the stiffness matrix
    meshptr->delete_stiff();
    
    cout<<__FILE__<<__LINE__<<" dual problem done "<<endl;
  }

  
  for ( int ifield = 0; ifield < meshptr->nfields; ifield++){
    meshptr->which_field = ifield;
    // gradient calculation for this deformation field.
    calcfuncgrad(meshptr);

    // write functional  for this field
    meshptr->writeonefieldfunctional();

    // write gradient for this field.
    meshptr->writeonefieldgradient(); 
  }


  // combine all the gradients into combined gradient vector
  meshptr->combine_gradients();

  // wrie the combined gradient vector to a file.
  meshptr->writegradient();

  // write the combined functional to a file
  meshptr->writefunctional();
  
  // feed this combined gradient vector to TAO using gradptr
  
  {
    int ieqn = 0;
    for (int inode = 0;inode < meshptr->num_prop_nodes; inode ++){
      if ( meshptr->opttype == 0 ) {
	for (int iprop = 0; iprop < meshptr->nprop; iprop++){
	  gradptr[ieqn] = meshptr->combgrad[inode][iprop];
	  ieqn++;
	}
      }
      else if (meshptr->opttype == 1) {
	gradptr[ieqn] = meshptr->combgrad[inode][0];
	ieqn++;
      }
      else if (meshptr->opttype == 2) {
	gradptr[ieqn] = meshptr->combgrad[inode][1];
	  ieqn++;
      }
	
    }
    
  }
  

  // put functional in *f and give it to TAO
  *f = meshptr->functional;
  
  // update inverse iteration flag;
  
  meshptr->inverseiteration++;

  // restore pointers ---- otherwise all hell breaks loose.

  VecRestoreArray(X,&stateptr); 
  VecRestoreArray(G,&gradptr);
  
  cout<<__FILE__<<__LINE__<<"done with gradientcallbackfortao.cpp"<<endl;

  return 0;
}
//====================================================================
// END OF FILE: Documentation follows.
/*[1]: If the forward problem is solved in one load increment and one newton iteration (as in a newton iteration), then the tangent in calc_forward and calc_dual are not the same. This is because the tangent in calc_forward is calcualted  about the zero deformation state and the tangent in calc_dual is calculated about the deformation state at the end of calc_forward. Thus when non-linearities are present and then two different tangents are calculated.
 */

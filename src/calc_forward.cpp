/**************************************************************************
* Nachiket Gokhale gokhalen@bu.edu
* 
* The non-linear forward problem 
*
**************************************************************************/

#ifndef HAVE_FUNCDECLS_H
#include "../headers/funcdecls.h"
#endif

#ifndef HAVE_MESH_HPP
#include "../lib/meshlib/mesh.hpp"
#endif

#ifndef HAVE_PETSC_HEADERS_H
#include "../headers/petsc_headers.h"
#endif

#ifndef HAVE_PETSC_SYSTEM_HPP
#include "../lib/eqnlib/petsc_system.hpp"
#endif

#ifndef HAVE_SLEPC_SYSTEM_HPP
#include "../lib/eqnlib/slepc_system.hpp"
#endif

#ifndef HAVE_SLEPC_QUAD_SYSTEM_HPP
#include "../lib/eqnlib/slepc_quad_system.hpp"
#endif 

#ifndef HAVE_GLOBALS_H
#include "../headers/globals.h"
#endif

int calc_forward(mesh* inmesh){
  
  // tells calc_forward howmany steps the load is actually applied in
  // depending on primalguessflag and inverseiteration number it gets 
  // set to either 1 (the entire load applied at once)
  // or inmesh->loadits (the load is applied gradually)

  int loadits;
  
  loadits = inmesh->loadits;
  
  // print properties if we're solving a forward problem
  // otherwise they're printed as they come out of TAO by propMonitor in inverse.cpp 
  if ( inmesh->inverseflag == false){
    inmesh->writeproperties();
  }

  // For the first inverse iteration, we need to solve 
  // a complete forward problem
  // After that we can solve it with no load steps.

  if(inmesh->primalguessflg == 1){ 
    //  the primal guess flag is set
    if (inmesh->inverseflag == true){
      // if we are solving an inverse problem
      if((inmesh->inverseiteration == 1)&&(inmesh->firstguessflag==0)){
	// solving an inverse problem at the first iteration
	// apply the load in steps 
	loadits = inmesh->loadits;
	inmesh->init_primal_goodguess_flag = false;
	inmesh->localizelast     = false;
      }
      else{
	// we are solving an inverse problem and 
	// are not at the first iteration or have forced first iteration guess
	// localize last loadstate
	loadits = 1;
	inmesh->localizelast = true;
	inmesh->init_primal_goodguess_flag = true;
      }
    }
    else{
      // we are solving a forward problem with the good guess flag
      // don't localize last iteration...helps if you want to increase deformation to beyond 
      // which a solution is available
      //loadits = 1;
      //inmesh->localizelast  = false;
      inmesh->localizelast = true;
      inmesh->init_primal_goodguess_flag = true;
    }
  }
  
  if(inmesh->linearflag == 1){
    inmesh->localizelast = true;
    loadits  = 1;
  }

  // before we start the solve: initialize the primal
  // so that it is reasonably close to expected soln
  // init_primal is responsible for setting the good guess
  // init_primal takes care of fdgradflag too 

  inmesh->init_primal();
  
  // the enhanced parameters are set by the constructor 
  // or by explicitly zeroing them
  // so when the inversion program calls this (calc_forward.cpp)
  // this starts with the values of the enhanced modes that are already present
  // currently they are being explicitly zeroed.

  if(inmesh->primalguessflg == 0){
    inmesh->zeroenhanced();
  }

  // forward part of the inverse problem
  for(int iloadits = 0;iloadits<loadits;iloadits++){
    
    // set the load iteration number in inmesh
    if(inmesh->localizelast == false){
      inmesh->iloadits = iloadits;
    }
    else{
      inmesh->iloadits = inmesh->loadits -1;
    }
    
    // need a better condition to terminate loop (when residual
    // is small then don't keep newton iterating.

    for(int inewton = 0;inewton<inmesh->newtonits;inewton++){
      // set the newton iteration number in inmesh;
      inmesh->inewtonits = inewton;
      equation_system *linsystem;

      slepc_system       *slepcsystem;
      petsc_system       *petscsystem;
      slepc_quad_system  *slepcquadsystem;

      equation_initdata eqninitdata;
      
      eqninitdata.num_equations = inmesh->num_equations[inmesh->which_field];

      if ( inmesh->issolvedispersion() ){
	eqninitdata.num_equations_real = eqninitdata.num_equations ;
	eqninitdata.num_equations_imag = eqninitdata.num_equations ;
	eqninitdata.which_pt           = inmesh->which_disperpt; 
	eqninitdata.which_seg          = inmesh->which_disperseg;
      }

      eqninitdata.which_field   = inmesh->which_field;
      eqninitdata.inewton       = inmesh->inewtonits;
      eqninitdata.printflag     = inmesh->printflag;
      eqninitdata.iload         = inmesh->iloadits;
      eqninitdata.max_nnz_per_row   = inmesh->max_nnz_per_row;

      cout<<endl<<"calc_forward.cpp: Number of equations for field "<<inmesh->which_field<<" == "<<inmesh->num_equations[inmesh->which_field]<<endl;
      
      if ( inmesh->solvetype == SOLVE_TYPE_BVP ){
	petscsystem =  (new petsc_system(eqninitdata));
	linsystem   =  petscsystem;
	cout<<endl<<"Allocated Petsc Linear System for BVP Solve"<<endl;
      }
      else if (inmesh->solvetype == SOLVE_TYPE_HOMOGENIZATION){
	petscsystem =  (new petsc_system(eqninitdata));
	linsystem   =  petscsystem;
	cout<<endl<<"Allocated Petsc Linear System for HOMOGENIZATION Solve"<<endl;	
      }
      else if (inmesh->issolvedispersion()) {
	// eigensolve
	if ( inmesh->solvetype != SOLVE_TYPE_DISPERSIONW2K ){
	  slepcsystem = (new slepc_system(eqninitdata));
	  slepcsystem->num_eigen = inmesh->num_eigen;
	  slepcsystem->setsolvetype(inmesh->solvetype);
	  linsystem = slepcsystem;
	  cout<<endl<<"Allocated SlepC Linear System"<<endl;
	}
	else{
	  int which_pt  = inmesh->which_disperpt;
	  int which_seg = inmesh->which_disperseg;
	  slepcquadsystem            = (new slepc_quad_system(eqninitdata));
	  slepcquadsystem->num_eigen = inmesh->num_eigen;
	  slepcquadsystem->ww        = inmesh->disperwvecs[which_pt][which_seg];
	  linsystem                  = slepcquadsystem;
	  cout<<endl<<"Allocated SlepC QUAD System"<<endl;	  
	}
      }
      
      // if we are doing newton's iterations consistently and 
      // build the consistent tangent matrix, calculate residual vector
      // if we are at the first newton iteration we need to build it.
      // why is this < 2 why not inmesh->ctangflag==1
	
      if((inmesh->ctangflag < 2)||(inewton==0)){
	buildtan(inmesh);
      }
      // give a pointer to the matrix to the PETSc routines
      
      linsystem->createcsrk(inmesh->stiff->row_ptr,inmesh->stiff->col_ptr,inmesh->stiff->value_arr);

      if (( inmesh->solvetype == SOLVE_TYPE_EIGENK) || (inmesh->solvetype== SOLVE_TYPE_EIGENM) || (inmesh->solvetype== SOLVE_TYPE_EIGEN)) {
	linsystem->createcsrm(inmesh->mass->row_ptr,inmesh->mass->col_ptr,inmesh->mass->value_arr);
	cout <<"Done with build createcsrm"<<endl;
      }

      if ( inmesh->solvetype == SOLVE_TYPE_DISPERSIONK2W ){
	// stiffness
	linsystem->createcsrkreal(inmesh->stiff->row_ptr,inmesh->stiff->col_ptr,inmesh->stiff->value_arr);
	linsystem->createcsrkimag(inmesh->stiffimag->row_ptr,inmesh->stiffimag->col_ptr,inmesh->stiffimag->value_arr);
	// mass 
	linsystem->createcsrmreal(inmesh->mass->row_ptr,inmesh->mass->col_ptr,inmesh->mass->value_arr);
 	linsystem->createcsrmimag(inmesh->massimag->row_ptr,inmesh->massimag->col_ptr,inmesh->massimag->value_arr);
      }

      if ( inmesh->solvetype == SOLVE_TYPE_DISPERSIONW2K){
	// mass 
	slepcquadsystem->createcsrmr(inmesh->quadmr->row_ptr,inmesh->quadmr->col_ptr,inmesh->quadmr->value_arr);
	// stiffness 
	slepcquadsystem->createcsrkr1(inmesh->quadkr1->row_ptr,inmesh->quadkr1->col_ptr,inmesh->quadkr1->value_arr);
	slepcquadsystem->createcsrkr2(inmesh->quadkr2->row_ptr,inmesh->quadkr2->col_ptr,inmesh->quadkr2->value_arr);
	slepcquadsystem->createcsrki(inmesh->quadki->row_ptr,inmesh->quadki->col_ptr,inmesh->quadki->value_arr);
	// damping 
	//slepcquadsystem->createcsrcreal(inmesh->quadcreal->row_ptr,inmesh->quadcreal->col_ptr,inmesh->quadcreal->value_arr);
 	//slepcquadsystem->createcsrcimag(inmesh->quadcimag->row_ptr,inmesh->quadcimag->col_ptr,inmesh->quadcimag->value_arr);
      }

      // COPY the rhs into PETSC
      if (( inmesh->solvetype == SOLVE_TYPE_BVP )||(inmesh->solvetype == SOLVE_TYPE_HOMOGENIZATION) ) {
	buildres(inmesh);
	inmesh->writeprimalrhs();
	linsystem->createrhs(inmesh->rhsvec);
      }
	
      // solve the matrix equation
      if ( inmesh->issolvedispersion() ){
	linsystem->solvecomplex();
	// get number of converged eigenpairs 
	inmesh->num_conv_eigen=min(linsystem->getconveigen(),(int)inmesh->num_eigen);
      }
      else{
	linsystem->solve();
	// get number of converged eigenpairs 
	inmesh->num_conv_eigen=min(linsystem->getconveigen(),(int)inmesh->num_eigen);
      }
      
      if ( inmesh->issolveeigen()){
	if ( inmesh->num_conv_eigen < inmesh->num_eigen ){
	  cout<<"======================================================================"<<endl;
	  cout<<"WARNING: Number of converged eigenpairs less than the number requested"<<endl;
	  cout<<"======================================================================"<<endl;
	}
      }
      
      // get solution and transfer it into inmesh
      if (( inmesh->solvetype == SOLVE_TYPE_BVP)||(inmesh->solvetype == SOLVE_TYPE_HOMOGENIZATION)){
	linsystem->getsol(inmesh->solrhsreal);
	inmesh->view_and_store_sol();
        inmesh->updategoodguess();
        inmesh->printrhsnorm();
      }
      
      if (( inmesh->solvetype == SOLVE_TYPE_EIGENK) || (inmesh->solvetype== SOLVE_TYPE_EIGENM) || (inmesh->solvetype== SOLVE_TYPE_EIGEN)) {
	for (int ieigen = 0; ieigen<inmesh->num_conv_eigen; ieigen++){
	  double kreal,kimag,kerror;
	  linsystem->geteigenvalue(&kreal,&kimag,&kerror,ieigen);
	  inmesh->eigenvalreal[ieigen][inmesh->which_field]= kreal;
	  inmesh->eigenvalimag[ieigen][inmesh->which_field]= kimag;
	  inmesh->eigenerror[ieigen][inmesh->which_field]  = kerror;
	}
	inmesh->writeeigenvalues();
      } // solvetype

      if ( inmesh->issolvedispersion()){
	for ( int ieigen = 0; ieigen<inmesh->num_conv_eigen; ieigen++){
	  double kreal,kimag,kerror;
	  linsystem->geteigenvalue(&kreal,&kimag,&kerror,ieigen);
	  linsystem->geteigenvectordispersion(inmesh->solrhsreal,inmesh->solrhsimag,ieigen);
	  inmesh->seteigendisper(ieigen,kreal,kimag,kerror);
	  if ( inmesh->printflag == 1 ) {
	    inmesh->writeeigenvectorsvtk(ieigen);
	  }
	}
      }
      
      delete linsystem;
      buildenh(inmesh);
      // print the enhanced strain parameters after we have calculated them
      inmesh->writeenhancedparm();

      // zero the tangent if we are doing newton consistently
      //if we are at the last newton iteration zero it!
      if((inmesh->ctangflag==1)||(inewton==inmesh->newtonits)){
	inmesh->stiff->zeromat();
      }

      inmesh->zerorhs();
      
      //skip the remaining iterations if rhsnorm< 1E-15
      if(inmesh->rhsnorm < pow(10.0,-12.00)){
      	inewton  = inmesh->newtonits;
      }
      
      // if the residual is not sufficiently small, then continue iterating
      if(inmesh->intelligentnewton == 1){
	if(inmesh->linearflag == 0){
	  if(inmesh->rhsnorm > pow(10.0,-7.00)){
	    if(inewton == (inmesh->newtonits -1)){
	      inewton = inewton -1;
	      cout<<"Warning --- Residual in primal not less than 1E-12 at end of newton iterations for field" <<inmesh->which_field<<endl;
	    }
	  }
	}
      }
      
    
    }
  }

  return 0;
}

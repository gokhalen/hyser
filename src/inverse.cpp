/***************************************************************************
 *
 * Nachiket Gokhale gokhalen@bu.edu
 * inverse.cpp: main driver file for the inverse problem.
 *
 * In this file, we create the TAO application and choose the method 
 * to be used for the optimization and the callback which provides TAO 
 * with the functional and the gradient vector. 
 **************************************************************************/


#include <stdio.h>
#include <stdlib.h>

#include <fstream>

#include "tao.h"

#ifndef HAVE_FUNCDECLS_H
#include "../headers/funcdecls.h"
#endif

#ifndef HAVE_PETSC_HEADERS_H
#include "../headers/petsc_headers.h"
#endif

#ifndef HAVE_TAO_HEADERS_H
#include "../headers/tao_headers.h"
#endif

#ifndef HAVE_GLOBALS_H
#include "../headers/globals.h"
#endif

#ifndef HAVE_MESH_HPP
#include "../lib/meshlib/mesh.hpp"
#endif

#ifndef HAVE_EQUATION_HPP
//#include "../lib/eqnlib/equation.hpp"
#endif

using namespace std;

static char help[] = "Inverse Hyperelastic Finite Element Code";

static PetscErrorCode set_inverse_problem_bounds(TaoSolver,Vec,Vec,void*);

int main(int argc,char **argv){
  // We first declare and intitalize variables that are required 
  // Adapted from the rosenbrock example in tao
   
  int info;
  // PetscScalar zero = 0.0;
  // PetscScalar one  = 1.0;
  
  //PetscScalar lowerbound  = 0.1;
  //PetscScalar upperbound  = 100.0;

  //PetscTruth flg;

  Vec statevars;
  // create vector of lower bounds (XL), upper bounds UB
  Vec XL, XU;

  TaoSolver taosol;
  //TAO_APPLICATION taoapp;
  TaoSolverTerminationReason reason;
  //TaoMethod  method = "tao_blmvm"; /* default minimization method */

  PetscErrorCode ierr;

  int size,rank;

  PetscViewer pview_taosol;

  PetscInitialize(&argc,&argv,(char *)0,help);
  TaoInitialize(&argc,&argv,(char *)0,help);

  // clean out files: system dependency;
  //system("find . -name 'armero*' |xargs rm -rf {}");
  //system("find . -name '*.m' |xargs rm -rf {}");
  //system("find . -name '*.out' |xargs rm -rf {}");

  system("rm -rf *.out armero* *.m");
  
  // Allocate memory for the mesh
  mesh *inmesh;
  
  // Read in the mesh
  // The material properties are set here, from data.in 
  // And used as the first guess

  inmesh = new mesh();
  inmesh->inverseflag = true;
  inmesh->init_primal();
  inmesh->readmeasdata();      // read in the measured data  
  //inmesh->make_fillin();     // make the fillin matrix (CSR) representatn
  inmesh->calculate_optvars(); // calulate number of optimization vars
  inmesh->resetprops();        // reset the properties based on optimization type

  
  if(inmesh->fdgradflag == 1){

    
      // calculate the gradient vector by finite differences
      fdgrad(inmesh);
      // write fdgradgrad for that field to a file
      //inmesh->writefdgradgrad();
      // combine the finite difference gradients
      inmesh->combine_fdgradients();
      // write the combined finite difference gradient
      inmesh->writecombfdgrad();
    
  }

  calc_normmeas(inmesh);

  info = MPI_Comm_size(PETSC_COMM_WORLD,&size); CHKERRQ(info);
  info = MPI_Comm_rank(PETSC_COMM_WORLD,&rank); CHKERRQ(info);
  
  if (size >1) {
    if (rank == 0)
      PetscPrintf(PETSC_COMM_SELF,"This example is intended for single processor use!\n");
    SETERRQ(PETSC_COMM_WORLD,1,"Incorrect number of processors");
  } 

  
  // create PETSC viewer
  PetscViewerASCIIOpen(PETSC_COMM_SELF,"taosol.out",&pview_taosol);
  
  
  // allocate memory for the solution vector, bound vector;
  VecCreateSeq(PETSC_COMM_SELF,inmesh->numoptvars,&statevars);
  VecCreateSeq(PETSC_COMM_SELF,inmesh->numoptvars,&XL);
  VecCreateSeq(PETSC_COMM_SELF,inmesh->numoptvars,&XU);

  //choose tao solver and create tao application
  TaoCreate(PETSC_COMM_WORLD,&taosol);
  ierr = TaoSetType(taosol,"tao_blmvm"); CHKERRQ(ierr);
  //TaoApplicationCreate(PETSC_COMM_SELF,&taoapp);



  // set guess from the inmesh properties
  // this can also be used to restart the inversion.
  // no need to use another file for restaring
  int counter = 0;
  for (int inode  = 0;inode<inmesh->num_prop_nodes;inode++){

      if ( inmesh->opttype == 0 ){
	for (int iprop = 0;iprop<inmesh->nprop; iprop++){
	  VecSetValue(statevars,counter,inmesh->prop[inode][iprop],INSERT_VALUES);
	  counter++;
	}
      }
      else if ( inmesh->opttype == 1 ){
	VecSetValue(statevars,counter,inmesh->prop[inode][0],INSERT_VALUES);
	counter++;
      }
      else if ( inmesh->opttype  == 2 ){
	VecSetValue(statevars,counter,inmesh->prop[inode][1],INSERT_VALUES);
	counter++;
      }
  }
  
  TaoSetInitialVector(taosol,statevars);
  TaoSetObjectiveAndGradientRoutine(taosol,gradientcallbackfortao,(void *)inmesh);
  //TaoAppSetInitialSolutionVec(taoapp,statevars);
  //TaoAppSetObjectiveAndGradientRoutine(taoapp,gradientcallbackfortao,(void *)inmesh); 
  //TaoAppSetMonitor(taoapp, propMonitor, (void *)inmesh);

  // set tolerances
  TaoSetTolerances(taosol,0,0.0,0.0,0.0,0.0);
  
  // need to set bounds on the statevars at some point of time
  // set bounds on the variables
  // set the first lower bound to the upper bound 
  // to enforce the constraint

  TaoSetVariableBoundsRoutine(taosol,set_inverse_problem_bounds,(void*)inmesh);
  TaoSetMaximumIterations(taosol,inmesh->taoiterates);

  //TaoAppSetVariableBoundsRoutine(taoapp,set_inverse_problem_bounds,(void*)inmesh);
  //TaoSetMaximumIterates(taosolver,inmesh->taoiterates);
  
  // solve the problem

  TaoSetFromOptions(taosol);
  TaoSolve(taosol);
  TaoGetTerminationReason(taosol,&reason);CHKERRQ(ierr);
  
  if (reason <= 0)
    PetscPrintf(MPI_COMM_WORLD,"Try a different TAO method, adjust some parameters, or check the function evaluation routines\n");
  ierr=TaoView(taosol,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);

  /* Free TAO data structures */
  info = TaoDestroy(&taosol); CHKERRQ(info);
  //info = TaoAppDestroy(taoapp); CHKERRQ(info);

  // view the solution (material props) as seen by TAO.
  VecView(statevars,pview_taosol);			  

  info = VecDestroy(&statevars); CHKERRQ(info);
  // Sept 26 - 2011 - untested.
  info = VecDestroy(&XL); CHKERRQ(info);
  info = VecDestroy(&XU); CHKERRQ(info);
  
  // destroy the solution vector
  PetscViewerDestroy(&pview_taosol);
  TaoFinalize();
  PetscFinalize();
  
  // destroy the inmesh object
  /*cout<<"deleting stiff in inverse.cpp"<<endl;
   inmesh->delete_stiff();*/
  delete inmesh;
  return 0;
}

// ===============================================================================
PetscErrorCode set_inverse_problem_bounds(TaoSolver taosol, Vec XL,Vec XU,void *ptr){

  // set bounds on the optimization variable
  // TaoAppSetVariableBoundsRoutine passes inmesh into this function
  // through the void pointer

  // cast void pointer to inmesh;
  mesh *inmesh;
  inmesh = (mesh *)(ptr);

  int proplocation = 0;
    for(int inode = 0; inode < inmesh->num_prop_nodes; inode++){
      if ( inmesh->opttype == 0 ) {
	for(int iprop = 0; iprop < inmesh->nprop; iprop++){
	  VecSetValue(XL,proplocation,inmesh->propcons[inode][iprop][0],INSERT_VALUES);
	  VecSetValue(XU,proplocation,inmesh->propcons[inode][iprop][1],INSERT_VALUES);
	  proplocation++;
	
	}
      }
      else if ( inmesh->opttype == 1){
	VecSetValue(XL,proplocation,inmesh->propcons[inode][0][0],INSERT_VALUES);
	VecSetValue(XU,proplocation,inmesh->propcons[inode][0][1],INSERT_VALUES);
	proplocation++;
      }
      else if ( inmesh->opttype == 2){
	VecSetValue(XL,proplocation,inmesh->propcons[inode][1][0],INSERT_VALUES);
	VecSetValue(XU,proplocation,inmesh->propcons[inode][1][1],INSERT_VALUES);
	proplocation++;
      }
  }
  return 0;
}

//======================================================================
PetscErrorCode propMonitor(TaoSolver taosol, void *ptr){
  // void *ptr is a pointer of type mesh
  // cast the void pointer to inmesh

  mesh *inmesh;
  inmesh = (mesh *)(ptr);
  inmesh->writeproperties();

  return 0;
}

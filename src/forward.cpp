#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#ifndef HAVE_FUNCDECLS_H
#include "../headers/funcdecls.h"
#endif

#ifndef HAVE_PETSC_HEADERS_H
#include "../headers/petsc_headers.h"
#endif


#ifndef HAVE_SLEPC_HEADERS_H
#include "../headers/slepc_headers.h"
#endif

#ifndef HAVE_GLOBALS_H
#include "../headers/globals.h"
#endif

#ifndef HAVE_MESH_HPP
#include "../lib/meshlib/mesh.hpp"
#endif

// After the SlepC/Petsc/Tao compiling issue has been resolved we can have a unified lib
#ifndef HAVE_EQUATION_HPP
//#include "../lib/eqnlib/equation.hpp"
#endif




static char help[] = "HyserVII Finite Element Code";

int main(int argc,char **argv){
  //mtrace();
  clock_t clock_start,clock_end;
  clock_start = clock();
  //PetscInitialize(&argc,&argv,(char *)0,help);
  SlepcInitialize(&argc,&argv,(char *)0,help);

  // remove all outfiles
  //system("find . -name 'armero*' |xargs rm -rf {}");
  //system("find . -name '*.m' |xargs rm -rf {}");
  //system("find . -name '*.out' |xargs rm -rf {}");
  //system("find . -name '*.vtk' |xargs rm -rf {}"); 

  system("rm -rf *.out *.vtk  armero*"); 
  
  
  //Allocate memory for the mesh
  mesh *inmesh;
  
  //Read in the mesh

  cout<<"forward.cpp: Started reading mesh"<<endl;

  inmesh = new mesh(); 
  inmesh->inverseflag = false;
  inmesh->initvtk("mesh.vtk");
  
  cout<<"forward.cpp: Finished reading Mesh"<<endl;


  int nseg_disper = 1; 
  if ( inmesh->solvetype == SOLVE_TYPE_DISPERSIONK2W ){
    nseg_disper = inmesh->ndispersionseg;
  }
  for ( int iseg = 0; iseg < nseg_disper; iseg++){
    inmesh->which_disperseg = iseg;
    
    int npt_disper = 1;

    if ( inmesh->issolvedispersion()){
      npt_disper = inmesh->ndispersionpts[iseg];
    }

    for ( int ipt = 0; ipt < npt_disper ; ipt++) {
      inmesh->which_disperpt = ipt;
      for(int ifield = 0; ifield < inmesh->nfields;ifield++){
	clock_t clock_eqn_start = clock();
	
	inmesh->which_field = ifield;
	// make fillin allocates a stiffness matrix
	inmesh->make_fillin();
	
	inmesh->init_primal();
      
	cout<<"forward.cpp: Calling Calc_forward"<<endl;
	calc_forward(inmesh);
	
	// delete the stiffness mass and damping matrices
	inmesh->delete_stiff();
	inmesh->delete_mass();
	inmesh->delete_damping();

	clock_t clock_eqn_end = clock();
	double eqn_time    = (clock_eqn_end - clock_eqn_start)/CLOCKS_PER_SEC;
	double eqn_timemin = eqn_time / 60.0;
	cout<< endl<<"Solve took "<<eqn_time<<" seconds "<< " or "<< eqn_timemin <<" minutes" <<endl;
      } // ifield
    } // ipt
  } // iseg
  
  // write vtk file 
  inmesh->writevtk();

  // this calculates and prints out the homogenized tensor 
  // based on hassani-and hinton formulation 
  // output to be taken seriously only if the correct boundary conditions and
  // displacement fields are imposed. 
  
  if ( inmesh->linearflag  == 1) {
    // zero the homogenized tensor -- this has to be done outside the fields loop
    zero4D(inmesh->ndime,inmesh->ndime,inmesh->ndime,inmesh->ndime,inmesh->EHIJKL);
  

    // homogenize
    if ( inmesh->solvetype == SOLVE_TYPE_HOMOGENIZATION ){
      for (int ifield = 0; ifield < inmesh->nfields ; ifield++){
	// this sets the field
	inmesh->which_field = ifield;
	// then call the homogenization code
	homogenize(inmesh);
      }
    }
    
    // compute volume and density only once; no looping over fields
    inmesh->computevolumemass();
    inmesh->computeaveprops();
    
    if ( inmesh->solvetype == SOLVE_TYPE_HOMOGENIZATION ){
      if (inmesh->ndime == 2 ){ inmesh->writeehijkl();   }
      if (inmesh->ndime == 3 ){ inmesh->writeehijkl3d(); }
    }

  }

  // write the dispersion curve 
  if ( inmesh->issolvedispersion()){
    inmesh->writedispersion();
  }

  // all arrays are lost
  delete(inmesh);
  //PetscFinalize();
  SlepcFinalize();
  
  clock_end = clock();
  double programtime    = (double)(clock_end-clock_start)/CLOCKS_PER_SEC;
  double programtimemin = programtime / 60.0;
  cout<< endl<<"Forward took "<<programtime<<" seconds "<< " or "<< programtimemin << " minutes "<<endl;
  return 0;
}


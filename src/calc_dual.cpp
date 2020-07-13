/*********************************************************************
 * Nachiket Gokhale gokhalen@bu.edu
 * April 04 2004
 *
 * This calculates the dual field 
 * 
 ********************************************************************/
#ifndef HAVE_FUNCDECLS_H
#include "../headers/funcdecls.h"
#endif

#ifndef HAVE_MESH_HPP
#include "../lib/meshlib/mesh.hpp"
#endif

#ifndef HAVE_PETSC_HEADERS_H
#include "../headers/petsc_headers.h"
#endif

#ifndef HAVE_EQUATION_HPP
//#include "../lib/eqnlib/equation.hpp"
#endif

#ifndef HAVE_PETSC_SYSTEM_HPP
#include "../lib/eqnlib/petsc_system.hpp"
#endif



int calc_dual(mesh *inmesh){
  

  petsc_system *petscdualsystem;
  petscdualsystem = new petsc_system(inmesh->num_equations[inmesh->which_field],inmesh->which_field, inmesh->inewtonits,inmesh->iloadits,inmesh->printflag);
  
  // print the solution vector.
  { 
    char filename[100];
    sprintf(filename,"solindcalcdual%d.out",inmesh->which_field);
    ofstream solindualfile(filename);
    solindualfile.precision(12);
    for (int inode = 0;inode<inmesh->num_nodes;inode++){
      for (int idofn =0;idofn<inmesh->mdofn; idofn++){
	solindualfile<<inmesh->solvec[inode][idofn][inmesh->which_field]<<"  ";
      }
      
      solindualfile<<endl;
    }
    solindualfile.close();
  }
  
  // hack: build the tangent for the final iteration (of forward prob) only
  // correct coice of load increment number and non-linear iteration number 
  buildtan(inmesh);
  builddualrhs(inmesh);


  // check the solution
  inmesh->stiff->viewcsr("csrincalcdual.out");
  petscdualsystem->createcsrk(inmesh->stiff->row_ptr,inmesh->stiff->col_ptr,inmesh->stiff->value_arr);
  petscdualsystem->createrhs(inmesh->rhsvec);
  petscdualsystem->solve();
  petscdualsystem->getsol(inmesh->solrhsreal);

  inmesh->storedual(); // stores the dual field 
  inmesh->writedual(); // writes the dual field


  // clean up
  delete petscdualsystem;
  inmesh->stiff->zeromat();
  inmesh->zerorhs();
  
  return 0;
}

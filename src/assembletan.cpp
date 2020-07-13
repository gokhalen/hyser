/*=========================================================================
  Nachiket Gokhale gokhalen@bu.edu 
  Assembly of a finite element matrix
  ============================================================================*/
#ifndef HAVE_FUNCDECLS_H
#include "../headers/funcdecls.h"
#endif

#ifndef HAVE_PETSC_HEADERS_H
#include "../headers/petsc_headers.h"
#endif

#ifndef HAVE_MESH_HPP
#include "../lib/meshlib/mesh.hpp"
#endif

#ifndef HAVE_GLOBALS_H
#include "../headers/globals.h"
#endif

int assembletan(mesh* inmesh,long int ielem,elem_base* felem){
  
  {
    int ieqn = 0;
    int jeqn = 0;

    // assumes mdofn as the number of degrees of freedom at each node
    PetscInt mrow=felem->shpstrct->numnodes*felem->mdofn;
    
    PetscInt *idxmrow;
    idxmrow =allocvec(mrow,idxmrow);

    PetscInt ncol=felem->shpstrct->numnodes*felem->mdofn;

    PetscInt *idxncol;
    idxncol=allocvec(ncol,idxncol);

    PetscScalar *vvs=allocvec(mrow*ncol,vvs);
    
    
    int ipos=0;
    for(int inode = 0;inode<felem->shpstrct->numnodes;inode++){
      for(int idofn=0;idofn<felem->dofnatnode[inode];idofn++){

	int jpos =0;
	for(int jnode = 0;jnode<felem->shpstrct->numnodes;jnode++){
	  for(int jdofn=0;jdofn<felem->dofnatnode[jnode];jdofn++){
	    ieqn = felem->ideqn[inode][idofn];
	    jeqn = felem->ideqn[jnode][jdofn];
	    
	    if((ieqn>0)&&(jeqn>0)){
	      
	      inmesh->stiff->addfill(ieqn-1,jeqn-1,felem->estiff[ipos][jpos]);

	      if ( inmesh->solvetype == SOLVE_TYPE_DISPERSIONK2W ) {
		inmesh->stiffimag->addfill(ieqn-1,jeqn-1,felem->estiffimag[ipos][jpos]);
		inmesh->mass->addfill(ieqn-1,jeqn-1,felem->ecmass[ipos][jpos]);
	      }

	      if ( inmesh->solvetype == SOLVE_TYPE_DISPERSIONW2K ){
		inmesh->quadkr1->addfill(ieqn-1,jeqn-1,felem->stiffkr1[ipos][jpos]);
		inmesh->quadkr2->addfill(ieqn-1,jeqn-1,felem->stiffkr2[ipos][jpos]);
		inmesh->quadki->addfill(ieqn-1,jeqn-1,felem->stiffki[ipos][jpos]);
		inmesh->quadmr->addfill(ieqn-1,jeqn-1,felem->ecmass[ipos][jpos]);
	      }

	    }
	    jpos++;
	    
	  }//jdofn
	}//jnode
	ipos++;
      } //idofn
      
    }//inode

    deallocvec(idxmrow);
    deallocvec(idxncol);
    deallocvec(vvs);

    
    //MatSetValues(*stiff,mrow,idxmrow,ncol,idxncol,vvs,ADD_VALUES);
  }

  
  
  return 0;
}


 


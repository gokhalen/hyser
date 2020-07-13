#ifndef HAVE_ELEMLIB_HPP
#include "../lib/elemlib/elemlib.hpp"
#endif

#ifndef HAVE_MESH_HPP
#include "../lib/meshlib/mesh.hpp"
#endif

#ifndef HAVE_GLOBALS_H
#include "../headers/globals.h"
#endif

elem_base* localizer(mesh *inmesh,int locielem,elem_base *felem, int loaditflag){
  // localizes element data 
  // localize constitutive equation number 
  felem->conseqnno = inmesh->conseqnno;
  felem->opttype   = inmesh->opttype;
  
  // localize propdiscontinflag gets used only in the calculation of the grad
  felem->propdiscontinflag = inmesh->propdiscontinflag;

  // localize element options
  felem->elemopt   = inmesh->lmref[locielem][2];

  // localize the volume
  felem->unitcellvolume = inmesh->unitcellvolume;
  felem->nfields    = inmesh->nfields;

  int which_loadit = inmesh->iloadits;
  
  if(inmesh->iloadits == 0){
    felem->condense_matrix = true;
  }
  else{
    felem->condense_matrix = true;
  }

  if ( inmesh->issolvedispersion() ) {
    int which_pt  = inmesh->which_disperpt;
    int which_seg = inmesh->which_disperseg;
  
    if ( inmesh->solvetype == SOLVE_TYPE_DISPERSIONK2W ){
      for ( int idime = 0; idime < felem->shpstrct->ndime ; idime++){
	double kvalue = inmesh->disperkvecs[idime][which_pt][which_seg];
	felem->setkvec(kvalue,idime); // sets the idimeth component of kvec in felem
      }
    } // if  k2w

    if ( inmesh->solvetype == SOLVE_TYPE_DISPERSIONW2K){
      double kx = inmesh->w2karray[which_seg][4];  felem->setkvec(kx,0); 
      double ky = inmesh->w2karray[which_seg][5];  felem->setkvec(ky,1); 
      double ww = inmesh->disperwvecs[which_pt][which_seg]; felem->setwdisp(ww);
    } // if w2k

  } // if dispersion 
  
  // localize the minimization norm 
  // and the regularization norm
  
  felem->which_functional     = inmesh->functionaltype;
  felem->which_regularization = inmesh->regfunctionaltype;
  
  if(inmesh->localizelast == true) {
    // calc_forward has explicitly asked to localize last load iteration
    // loaditflag gets set in calc_forward;
    loaditflag    = LAST_LOAD_ITERATION;
  }

  if(loaditflag == LAST_LOAD_ITERATION){
    which_loadit = inmesh->loadits-1;
  }

  // set coordinate data and material properties
  for(int inode = 0;inode < felem->shpstrct->numnodes;inode++){
    int globalnode = inmesh->ien[locielem][inode];
    
    // localize primal 

    double primval;
    
    for(int idofn = 0;idofn<felem->dofnatnode[inode];idofn++){
      //if dirichlet
      if(inmesh->isdirichlet(globalnode-1,idofn,inmesh->which_field)){
	primval = inmesh->getdirichlet(globalnode-1,idofn,which_loadit,inmesh->which_field);
      }
      else{
	primval = inmesh->solvec[globalnode-1][idofn][inmesh->which_field];
	if(loaditflag == ENHANCED){
	  primval = primval - inmesh->lastupdate[globalnode-1][idofn][inmesh->which_field];
	}
      }
      felem->setprimal(primval,inode,idofn);      
      double lastup = inmesh->lastupdate[globalnode-1][idofn][inmesh->which_field];
      felem->setlastupdate(lastup,inode,idofn);
    }//idofn
    

    // localize primalwhole
    for ( int ifield=0; ifield < inmesh->nfields; ifield++){
      for(int idofn = 0;idofn<felem->dofnatnode[inode];idofn++){

	if( inmesh->isdirichlet(globalnode-1,idofn,ifield)){
	  felem->primalwhole[inode][idofn][ifield] = inmesh->getdirichlet(globalnode-1,idofn,which_loadit,ifield);
	}
	else{
	  felem->primalwhole[inode][idofn][ifield] = inmesh->solvec[globalnode-1][idofn][ifield];
	}

      }//idofn    
    }//ifield

    
    
    // assume material property nodes are same as the element nodes
    // and therefore, isoparametrically interpolated. can be changed..
    
    
    for(int iprop =0;iprop < felem->nprop;iprop++){
      if(inmesh->propdiscontinflag == 0){
	felem->setprop(inmesh->prop[globalnode-1][iprop],iprop,inode);
      }
      if(inmesh->propdiscontinflag == 1){
	felem->setprop(inmesh->prop[locielem][iprop],iprop,inode);
      }
    }

    // set density 
    if ( inmesh->propdiscontinflag == 0 ) {
      // continuous 
      felem->setdensity(inmesh->density[globalnode-1][0],inode);
    }
    if ( inmesh->propdiscontinflag == 1 ) {
      // discontinuous - elementwise constant
      felem->setdensity(inmesh->density[locielem][0],inode);
    }

    // localize coordinates
    for(int idime = 0;idime < felem->shpstrct->ndime;idime++){
      felem->setcoord(inmesh->coord[globalnode-1][idime],idime,inode);

      if((inmesh->lmref[locielem][0] == 261)&&(inmesh->lmref[locielem][2]==2)){
	double coordspat = inmesh->coord[globalnode-1][idime] + primval;
	felem->setcoord(coordspat,idime,inode);
      }
      
    }


  }//inode

  // localize enhancements for piltner-taylor and other elements which need them
  {
    if ( inmesh->isenhanced(locielem)){
      for(int ii = 0; ii < 4; ii++){
	felem->enh[ii] = inmesh->enhvec[locielem][ii][inmesh->which_field];
      }
    }
  } 
  

  //  abort();
  
  // localize fem arrays
  //  first the id array (local node number to global node number)
  for(int inode =0;inode< felem->shpstrct->numnodes;inode++){
    felem->setidnode(inode,inmesh->conn[locielem][inode]);
  }
    
  felem->viewid();
  //localize boundary data


  // localize the type of initial strain forcing
  felem->initstrntype = inmesh->initstrn[inmesh->which_field];
 
  // boundary data localization. 
  // notice that this accounts for variable degrees of freedom

  for(int inode=0;inode<felem->shpstrct->numnodes;inode++){
    for(int idofn=0;idofn<felem->dofnatnode[inode];idofn++){
      
      int globalnode = felem->idnode[inode];
      // arrays zeroed by constrctor
      //check if BC is specified at all
      //if dirichlet
      if(inmesh->isdirichlet(globalnode-1,idofn,inmesh->which_field)){
	
	felem->setisbc(inode,idofn,1);
	double dirval = inmesh->getdirichlet(globalnode-1,idofn,which_loadit,inmesh->which_field);
	felem->setdirdata(dirval,inode,idofn);
      }
	
      // traction
      if(inmesh->istraction(globalnode-1,idofn,inmesh->which_field)){
	felem->setisbc(inode,idofn,-1);
	double tracval = inmesh->gettraction(globalnode-1,idofn,which_loadit,inmesh->which_field);
	felem->settracdata(tracval,inode,idofn);
      }
    } // idofn
  } //inode

  // localize equation number 
  for(int inode = 0;inode<felem->shpstrct->numnodes;inode++){
    int globalnode = inmesh->ien[locielem][inode];
    
    for(int idofn = 0;idofn<felem->dofnatnode[inode];idofn++){
      
      felem->setideqn(inode,idofn,inmesh->id[globalnode-1][idofn][inmesh->which_field]);
    }
  }

  
  // localize the measured data

  for(int inode = 0;inode<felem->shpstrct->numnodes;inode++){
    int globalnode = inmesh->ien[locielem][inode];
    for(int idofn = 0;idofn<felem->dofnatnode[inode];idofn++){
      if(inmesh->inverseflag == true){
	double measval = inmesh->measdata[globalnode-1][idofn][inmesh->which_field];
	felem->setmeasdata(measval,inode,idofn);
      }
    }
  }

  // localize the functional modifying tensor 
  for (int idime = 0; idime < inmesh->ndime; idime++){
    for(int jdime =0; jdime < inmesh->ndime; jdime++){
      felem->functens[idime][jdime] = inmesh->functens[idime][jdime];
    }
  }


  // This function should be called after both measured data
  // and the primal data  and functens are localized
  // Calculate Tmeas and Tprim
  felem->calculate_Tmeas_Tprim();
    

  // localize the dual field
  // dual must be stored in mesh which it is not being stored.

  for(int inode = 0; inode<felem->shpstrct->numnodes;inode++){
    int globalnode = inmesh->ien[locielem][inode];
    for(int idofn  = 0;idofn<felem->dofnatnode[inode];idofn++){
      if(inmesh->inverseflag == true){
	double dualval = inmesh->dualfield[globalnode-1][idofn][inmesh->which_field];
	felem->setdualfield(dualval,inode,idofn);

	for (int ifield = 0; ifield < inmesh->nfields; ifield++){
	  felem->dualwhole[inode][idofn][ifield]   =  inmesh->dualfield[globalnode-1][idofn][ifield];
	}
      }
    }
  }

 
  
  // localize the regularization parameter
  for(int iprop = 0; iprop < inmesh->nregprop; iprop++){
    for(int iparam = 0; iparam < inmesh->mregparm; iparam++){
      felem->regparm[iprop][iparam] = inmesh->regparm[iprop][iparam];
    }
  }

  // localize the functionalweight
  felem->funcweight = inmesh->fieldfuncweight[inmesh->which_field];


  // locallize the measured material property tensor
   for (int ii = 0; ii < inmesh->ndime; ii++){
      for (int  jj = 0; jj < inmesh->ndime; jj++){
	for (int kk = 0; kk < inmesh->ndime; kk++){
	  for (int ll = 0; ll < inmesh->ndime; ll++){
	    felem->BASE_PMIJKL[ii][jj][kk][ll] = inmesh->PMIJKL[ii][jj][kk][ll];
	    felem->BASE_EHIJKL[ii][jj][kk][ll] = inmesh->EHIJKL[ii][jj][kk][ll];
	  }
	}
      }
   }
  
  
  
  return felem;
}

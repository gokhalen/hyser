/*======================================================================
 *  Nachiket Gokhale gokhalen@bu.edu
 *
 *  Two dimensional Linear Elasticity Element
 *  Builds stiffness matrix and right hand side
 *
 **========================================================================*/
#ifndef HAVE_NONLINELASTRAC2D_HPP
#include "nonlinelastrac2d.hpp"
#endif

#include <cassert>
#include <cstdlib>

using namespace std;


// implement constrctor 

nonlinelastrac2d::nonlinelastrac2d(int *cinte,int cnprop, int cmregparm, long int celemno){
  this->nprop = cnprop;
  this->mregparm = cmregparm;
  mdofn = 2;
  nonlinflag = true; // why is the purpose of this flag?
  shpstrct = new shapestruct(2,1,cinte,2,NONLINELASTRAC2D,celemno);
  
  // allocate the rest after number of nodes has been set
  prop        = allocmat(shpstrct->numnodes,nprop,prop); 
  estiff      = allocmat(shpstrct->numnodes*mdofn,shpstrct->numnodes*mdofn,estiff);
  dofnatnode  = allocvec(shpstrct->numnodes,dofnatnode);

  // set the dofnatnode array, one degree of freedom at each node
  for(int inode = 0;inode<shpstrct->numnodes;inode++){
    dofnatnode[inode] = mdofn;
  }

  // allocate primal
  primal  = allocmat(shpstrct->numnodes,mdofn,primal);
  zeromat(shpstrct->numnodes,mdofn,primal);

  // allocate lastupdate
  lastupdate = allocmat(shpstrct->numnodes,mdofn,lastupdate);
  zeromat(shpstrct->numnodes,mdofn,lastupdate);

  // allocate dirichlet data and traction data vectors
  dirdata  = allocmat(shpstrct->numnodes,mdofn,dirdata);
  zeromat(shpstrct->numnodes,mdofn,dirdata);

  tracdata = allocmat(shpstrct->numnodes,mdofn,tracdata);
  zeromat(shpstrct->numnodes,mdofn,tracdata);

  // allocate right hand side vector
  erhs = allocvec(shpstrct->numnodes*mdofn,erhs);
  

  // allocate finite element arrays
  ideqn  = allocmat(shpstrct->numnodes,mdofn,ideqn); // (node num, dofn) -> glb eqn no
  idnode = allocvec(shpstrct->numnodes,idnode); // local -> global node num.
  isbc   = allocmat(shpstrct->numnodes,mdofn,isbc); //
  zeromat(shpstrct->numnodes,mdofn,isbc);


  // allocate space for inversion data
  measdata  = allocmat(shpstrct->numnodes,mdofn,measdata);
  zeromat(shpstrct->numnodes,mdofn,measdata);

  // allocate dual field data
  dualfield = allocmat(shpstrct->numnodes,mdofn,dualfield);
  zeromat(shpstrct->numnodes,mdofn,dualfield);

  // allocate gradient matrix
  
  gradient = allocmat(shpstrct->numnodes,nprop,gradient);
  zeromat(shpstrct->numnodes,nprop,gradient);

  // allocate the regparm array
  // can be done in the base class once the virtual function stuff
  // is figured out
  //regparm  = allocvec(this->nprop,regparm);

  
  functens = allocmat(shpstrct->ndime+1,shpstrct->ndime+1,functens);

  Tprim    = allocmat(shpstrct->numnodes,shpstrct->ndime,Tprim);
  Tmeas    = allocmat(shpstrct->numnodes,shpstrct->ndime,Tmeas);

  zeromat(shpstrct->numnodes,shpstrct->ndime,Tprim);
  zeromat(shpstrct->numnodes,shpstrct->ndime,Tmeas);
  
  memalloc();
}
//==================================================================
nonlinelastrac2d::~nonlinelastrac2d(){
  deallocmat(shpstrct->numnodes,primal);
  deallocmat(shpstrct->numnodes,lastupdate);
  deallocmat(shpstrct->numnodes,prop);
  deallocmat(shpstrct->numnodes*mdofn,estiff);
  deallocvec(dofnatnode);
  deallocmat(shpstrct->numnodes,dirdata);
  deallocmat(shpstrct->numnodes,tracdata);
  deallocvec(erhs);

  deallocmat(shpstrct->numnodes,ideqn);
  deallocvec(idnode);
  deallocmat(shpstrct->numnodes,isbc);
  //shpstrct->~shapestruct();
  
   //deallocate measdata
  deallocmat(shpstrct->numnodes,measdata);

  // deallocate dualfield
  deallocmat(shpstrct->numnodes,dualfield);

  // deallocate gradient 
  deallocmat(shpstrct->numnodes,gradient);

  //deallocvec(regparm);
  deallocmat(shpstrct->ndime+1,functens);
  deallocmat(shpstrct->numnodes,Tprim);
  deallocmat(shpstrct->numnodes,Tmeas);
  
  
  delete shpstrct;
}
//================================================================
void nonlinelastrac2d::elemstiff(int which_stiffness){

  zeromat(shpstrct->numnodes*mdofn,shpstrct->numnodes*mdofn,estiff);
  
    


}//closes elemstiff
//====================================================================
void nonlinelastrac2d::calcfunc(){
  this->functional = 0.0;
  this->regfunctional = 0.0;
}
//====================================================================
void nonlinelastrac2d::elemrhs(){
  // zero right hand side vector
  // The equation that we solve is 
  // K_{n}\delta{U}_{n+1} = -N_{n}
  // where n is the residual and is defined by
  
  zerovec(shpstrct->numnodes*mdofn,erhs);

  double tractrans[2][2];
  
  // Cook's membrane hack
  // need to figure out inclination with Y direction
  
  double x1 = shpstrct->coord[0][0] + primal[0][0];
  double y1 = shpstrct->coord[0][1] + primal[0][1];
  double x2 = shpstrct->coord[1][0] + primal[1][0];
  double y2 = shpstrct->coord[1][1] + primal[1][1];

  double length = (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2);
  length = sqrt(length);

  double lengthratio = fabs(y2 - y1)/length;

  double theta = atan2(y2-y1,x2-x1);
  double tracmag  = tracdata[0][1];
  double newtrac[2];

  newtrac[0] = tracmag*cos(theta);
  newtrac[1] = tracmag*sin(theta);
  
  
  for(int iinte = 0; iinte < shpstrct->tinte; iinte++){
    double jacodet = shpstrct->jdet[iinte];
    double weight  = shpstrct->wts[iinte];

    for(int idofn = 0; idofn<mdofn; idofn++){
      for(int inode =0; inode<shpstrct->numnodes; inode++){    
	//selected the the node and dofn of the the weight
	  int ieqn = (inode -1 + 1)*mdofn + idofn;

	  

	  double temptrac = 0;

	  if((elemopt == 0)||(elemopt == 2)){
	    for (int jnode =0; jnode<shpstrct->numnodes;jnode++){
	      temptrac += shpstrct->shape[iinte][jnode]*tracdata[jnode][idofn];
	    }
	  }
	  else if(elemopt == 1){
	    for (int jnode =0; jnode<shpstrct->numnodes;jnode++){
	      temptrac += shpstrct->shape[iinte][jnode]*newtrac[idofn];	   
	    }
	  }
	  
	  
	  erhs[ieqn] += shpstrct->shape[iinte][inode]*temptrac*jacodet*weight;
	//cout<<"tracdata == "<<tracdata[inode][idofn]<<" "<<__FILE__<<__LINE__<<endl;
      }//idofn
    }//inode
    
  }//iinte

  //for(int ii = 0; ii<4;ii++){
  //  cout<<__FILE__<<" "<<__LINE__<<" "<<erhs[ii]<<endl;
  //}



}//elemrhs()
//====================================================================
void nonlinelastrac2d::addreggradient(){
  for(int inode = 0; inode < shpstrct->numnodes; inode++){
      for(int iprop = 0; iprop < nprop; iprop++){
	this->gradient[inode][iprop] = 0.0;
      }
  }
}
//===================================================================
void nonlinelastrac2d::calcgradient(){
  for(int inode = 0; inode < shpstrct->numnodes; inode++){
    for(int iprop = 0; iprop < nprop; iprop++){
      this->gradient[inode][iprop] = 0.0;
    }
  }
}
//====================================================================
void nonlinelastrac2d::calcdual(){
   for (int pdime = 0; pdime < shpstrct->ndime; pdime++){
     for (int pnode =0; pnode < shpstrct->numnodes; pnode++){
       int peqn = mdofn*(pnode - 1 + 1) + pdime; 
       this->erhs[peqn] = 0.0;
      }
   }
}
//====================================================================
void nonlinelastrac2d::calcnormmeas(){
  this->normmeasdata = 0;
}
//====================================================================
void nonlinelastrac2d::getinteg(){
  gauss1d(shpstrct);
}
//====================================================================
void nonlinelastrac2d::getshape(){
  shape1d(shpstrct);
}
//====================================================================
void nonlinelastrac2d::printeltype(){
  //cout<<endl<<"This is NONLINELASTRAC2D"<<endl;
}
//====================================================================
//void nonlinelastrac2d::memalloc(){

  // inherited
//}
//====================================================================
void nonlinelastrac2d::memdealloc(){
  // dealloc memory
}
//===================================================================

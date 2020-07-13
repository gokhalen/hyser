/*======================================================================
 *  Nachiket Gokhale gokhalen@bu.edu
 *
 *  Two dimensional Linear Elasticity Element
 *  Builds stiffness matrix and right hand side
 *
 **========================================================================*/
#ifndef HAVE_LINELASTRAC2D_HPP
#include "linelastrac2d.hpp"
#endif

#include <cassert>
#include <cstdlib>

using namespace std;


// implement constrctor 

linelastrac2d::linelastrac2d(int *cinte,int cnprop, int cmregparm,long int celemno){
  this->nprop = cnprop;
  this->mregparm = cmregparm;

  mdofn = 2;
  nonlinflag = false;
  shpstrct = new shapestruct(2,1,cinte,2,LINELASTRAC2D,celemno);
  
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
linelastrac2d::~linelastrac2d(){
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
void linelastrac2d::elemstiff(int which_stiffness){

  zeromat(shpstrct->numnodes*mdofn,shpstrct->numnodes*mdofn,estiff);
  
    


}//closes elemstiff
//====================================================================
void linelastrac2d::calcfunc(){
  this->functional = 0.0;
  this->regfunctional = 0.0;
}
//====================================================================
void linelastrac2d::elemrhs(){
  // zero right hand side vector
  zerovec(shpstrct->numnodes*mdofn,erhs);
  
  for(int iinte = 0; iinte < shpstrct->tinte; iinte++){
    double jacodet = shpstrct->jdet[iinte];
    double weight  = shpstrct->wts[iinte];

    for(int idofn = 0; idofn<mdofn; idofn++){
      for(int inode =0; inode<shpstrct->numnodes; inode++){

	int ieqn = (inode -1 + 1)*mdofn + idofn;
	double temptrac = 0;
	for (int jnode =0; jnode<shpstrct->numnodes;jnode++){
	  temptrac += shpstrct->shape[iinte][jnode]*tracdata[jnode][idofn];
	}
	  
	
	erhs[ieqn] += shpstrct->shape[iinte][inode]*temptrac*jacodet*weight;
	
	//cout<<"tracdata == "<<tracdata[inode][idofn]<<" "<<__FILE__<<__LINE__<<endl;
      }//idofn
    }//inode
    
  }//iinte
}//elemrhs()
//====================================================================
//====================================================================
void linelastrac2d::getinteg(){
  gauss1d(shpstrct);
}
//====================================================================
void linelastrac2d::getshape(){
  shape1d(shpstrct);
}
//====================================================================
void linelastrac2d::printeltype(){
  //cout<<endl<<"This is LINELASTRAC2D"<<endl;
}
//====================================================================
//void linelastrac2d::memalloc(){

  // inherited
//}
//====================================================================
void linelastrac2d::memdealloc(){
  // dealloc memory
}
//===================================================================
void linelastrac2d::addreggradient(){
  for(int inode = 0; inode < shpstrct->numnodes; inode++){
      for(int iprop = 0; iprop < nprop; iprop++){
	this->gradient[inode][iprop] = 0.0;
      }
  }
}
//===================================================================
void linelastrac2d::calcgradient(){
  for(int inode = 0; inode < shpstrct->numnodes; inode++){
    for(int iprop = 0; iprop < nprop; iprop++){
      this->gradient[inode][iprop] = 0.0;
    }
  }
}
//====================================================================
void linelastrac2d::calcdual(){
   for (int pdime = 0; pdime < shpstrct->ndime; pdime++){
     for (int pnode =0; pnode < shpstrct->numnodes; pnode++){
       int peqn = mdofn*(pnode - 1 + 1) + pdime; 
       this->erhs[peqn] = 0.0;
      }
   }
}
//====================================================================
void linelastrac2d::calcnormmeas(){
  this->normmeasdata = 0;
}

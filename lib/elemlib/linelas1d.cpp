/*======================================================================
  Nachiket Gokhale gokhalen@bu.edu
  One Dimensional Linear Elasticity Element
  Builds stiffness matrix and right hand side
=======================================================================*/

#ifndef HAVE_LINELAS1D_HPP
#include "linelas1d.hpp"
#endif

#include <cassert>
#include <cstdlib>

//=========================================================================
//implement constructor
linelas1d::linelas1d(int* cinte,long int celemno){
  // have all information necessary. call constrctor of shapestrct
  nprop = 1;
  mdofn = 1;
  nonlinflag = false;
  
  shpstrct = new shapestruct(1,1,cinte,2,LINELAS1D,celemno);

  // allocate the rest after number of nodes has been set
  prop       = allocmat(shpstrct->numnodes,nprop,prop); 
  estiff      = allocmat(shpstrct->numnodes*mdofn,shpstrct->numnodes*mdofn,estiff);
  dofnatnode = allocvec(shpstrct->numnodes,dofnatnode);

  // set the dofnatnode array, one degree of freedom at each node
  for(int inode = 0;inode<shpstrct->numnodes;inode++){
    dofnatnode[inode] = 1;
  }

  // allocate primal
  primal  = allocmat(shpstrct->numnodes,mdofn,primal);
  zeromat(shpstrct->numnodes,mdofn,primal);

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
  }
//==========================================================================
linelas1d::~linelas1d(){

  deallocmat(shpstrct->numnodes,prop);
  deallocmat(shpstrct->numnodes*mdofn,estiff);
  deallocvec(dofnatnode);
  deallocmat(shpstrct->numnodes,dirdata);
  deallocmat(shpstrct->numnodes,tracdata);
  deallocmat(shpstrct->numnodes,ideqn);
  deallocvec(idnode);
  deallocmat(shpstrct->numnodes,isbc);
  deallocmat(shpstrct->numnodes,primal);
  //shpstrct->~shapestruct();
  delete shpstrct;
}
//========================================================================
// implement element stiffness routine
void linelas1d::elemstiff(int which_stiff){
  
  // zero element stiffness K_{A,B} = \int_{\Omega}(N_{A,x}k(x)N_{B,X})
  zeromat(shpstrct->numnodes*mdofn,shpstrct->numnodes*mdofn,estiff);

  // loop over integration points
  for(int iinte = 0;iinte<shpstrct->tinte;iinte++){

    // get jacobian determinant and weight
    
    double jacodet = shpstrct->jdet[iinte];
    double weight  = shpstrct->wts[iinte];

    // interpolate stiffness at this integration point
    double stiff = 0;
    for(int knode = 0;knode<shpstrct->numnodes;knode++){
      stiff +=  shpstrct->shape[iinte][knode]*prop[knode][0];      
    }


    for(int inode =0;inode<shpstrct->numnodes;inode++){     // pick the node corresp. A
      double der1 = shpstrct->glbder[iinte][inode][0];

      for(int jnode=0;jnode<shpstrct->numnodes;jnode++){     // pick the node corresp. B
	double der2 = shpstrct->glbder[iinte][jnode][0];
	
	estiff[inode][jnode] += der1*stiff*der2*jacodet*weight;
      }
    }
  }

}

//=======================================================================
void linelas1d::calcfunc(){
  cout<<"calcfunc in linelas1d()"<<endl;
}

//=======================================================================
void linelas1d:: elemrhs(){

  // zero righthand side vector
  zerovec(shpstrct->numnodes*mdofn,erhs);
  // elemstiff should be called before calling elemrhs
  // dirichlet data: take product of matrix with right hand side

  for(int inode = 0; inode < shpstrct->numnodes;inode++){
    for(int jnode = 0; jnode< shpstrct->numnodes;jnode++){
      erhs[inode] -= estiff[inode][jnode]*dirdata[jnode][0];
    }
  }

  // neumann data

}

//=======================================================================
void linelas1d:: calcgradient(){
  cout<<"calcgradient in linelas1d()"<<endl;
}

//======================================================================
void linelas1d::getinteg(){
  // get integration points
  gauss1d(shpstrct);
}
//====================================================================
void linelas1d::getshape(){
  //cout<<"getshape in linelas1d()"<<endl;
  shape1d(shpstrct);
}
//======================================================================
void linelas1d::printeltype(){
  cout <<endl<<"This is LINELAS1D"<<endl;
}
//=====================================================================
void linelas1d::memalloc(){
  // material properties, connectivity, stiffnessmatrix, righthand side
  // ideqn,idnode,isbc
  // The best way is to inherit this..
  
}
//=====================================================================
void linelas1d::memdealloc(){
  // dealloc memory

}
//=====================================================================
void linelas1d::setcoord(double cval,int idime,int inode){
  // cval = coordinate value
  // inode = local nodenumber
  // idime = dimension
  // Again should, be in base class
  shpstrct->coord[inode][idime] = cval;
}
//====================================================================
void linelas1d::setprop(double propval,int iprop,int inode){
  // cval = coordinate value
  // inode = local nodenumber
  // idime = dimension
  // Again should, be in base class
  prop[inode][iprop] = propval;
}
//===================================================================
void linelas1d::setideqn(int inode,int idofn,long int eqnno){

  assert(inode<shpstrct->numnodes);  // 
 
 
  assert(idofn<dofnatnode[inode]);   // 
  


  ideqn[inode][idofn] = eqnno;
}
//==================================================================
void linelas1d::setidnode(int inode,long int globalnode){
  assert(inode<shpstrct->numnodes);
  idnode[inode]= globalnode;
}
//=================================================================
void linelas1d::setisbc(int inode,int idofn,int value){
  assert(inode<shpstrct->numnodes);
  assert(idofn<dofnatnode[inode]);
  
  isbc[inode][idofn] = value;

}
//==============================================================
void linelas1d::setdirdata(double ddata,int inode,int idofn){
   assert(inode<shpstrct->numnodes);
   assert(idofn<dofnatnode[inode]);
   
   dirdata[inode][idofn] = ddata;
}
//==============================================================
void linelas1d::settracdata(double tdata,int inode,int idofn){
  assert(inode<shpstrct->numnodes);
  assert(idofn<dofnatnode[inode]);

  tracdata[inode][idofn] = tdata;
}
//=============================================================

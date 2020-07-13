/*======================================================================
  Nachiket Gokhale gokhalen@bu.edu
  Two Dimensional Non-Linear Heat Conduction Element
  Builds stiffness matrix and right hand side
=======================================================================*/

#ifndef HAVE_NONLINHEAT2D_HPP
#include "nonlinheat2d.hpp"
#endif

#include <cassert>
#include <cstdlib>

//=========================================================================
//implement constructor
nonlinheat2d::nonlinheat2d(int* cinte,long int celemno){
  // have all information necessary. call constrctor of shapestrct
  nprop = 3;
  mdofn = 1;
  nonlinflag = true;

  shpstrct = new shapestruct(2,2,cinte,4,NONLINHEAT2D,celemno);

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
nonlinheat2d::~nonlinheat2d(){

  deallocmat(shpstrct->numnodes,prop);
  deallocmat(shpstrct->numnodes*mdofn,estiff);
  deallocvec(dofnatnode);
  deallocmat(shpstrct->numnodes,dirdata);
  deallocmat(shpstrct->numnodes,tracdata);
  deallocmat(shpstrct->numnodes,ideqn);
  deallocvec(idnode);
  deallocmat(shpstrct->numnodes,isbc);
  //shpstrct->~shapestruct();
  delete shpstrct;
}
//========================================================================
// implement element stiffness routine
void nonlinheat2d::elemstiff(int which_stiffness){
  
  // zero element stiffness K_{A,B} = \int_{\Omega}(N_{A,x}k(x)N_{B,X})
  zeromat(shpstrct->numnodes*mdofn,shpstrct->numnodes*mdofn,estiff);
  zerovec(shpstrct->numnodes*mdofn,erhs);

  // declarations specific to nonlinear problems;
  double *KK;
  double **defgrad;
  double **gradu;
  double **secondpiola;
  double **rightcauchygreen;
  double **rightcauchygreeninvtpose;


  // allocations, deallocated at the end of elemstiff()
  KK                  = allocvec(this->nprop,KK);
  defgrad             = allocmat(shpstrct->ndime,shpstrct->ndime,defgrad);
  gradu               = allocmat(shpstrct->ndime,shpstrct->ndime,gradu);
  secondpiola         = allocmat(shpstrct->ndime,shpstrct->ndime,secondpiola);
  rightcauchygreen    = allocmat(shpstrct->ndime,shpstrct->ndime,rightcauchygreen);
  rightcauchygreeninvtpose = allocmat(shpstrct->ndime,shpstrct->ndime,rightcauchygreeninvtpose);

  // loop over integration points
  for(int iinte = 0;iinte<shpstrct->tinte;iinte++){

    // get jacobian determinant and weight
    
    double jacodet = shpstrct->jdet[iinte];
    double weight  = shpstrct->wts[iinte];

    // interpolate conductivities at this integration point
    zerovec(this->nprop,KK);

    for(int iprop = 0;iprop<this->nprop;iprop++){
      for(int knode = 0;knode<shpstrct->numnodes;knode++){
	KK[iprop] +=  shpstrct->shape[iinte][knode]*prop[knode][iprop];      
      }
    }

    // calculate temperature at this point
    double temperature = 0;
    for(int inode = 0;inode<shpstrct->numnodes;inode++){
      temperature += shpstrct->shape[iinte][inode]*primal[inode][0];
    }

    // calculate D
    double DD[2][2];
    for(int idime = 0;idime<shpstrct->ndime;idime++){
      for(int jdime =0;jdime<shpstrct->ndime;jdime++){
	DD[idime][jdime] = delta(idime,jdime)*(KK[0] + KK[1]*temperature + KK[2]*temperature*temperature);
      }
    }

    // calcualte gradient of D wrt temp
    double DDGRD[2][2];
    for(int idime = 0;idime<shpstrct->ndime;idime++){
      for(int jdime =0;jdime<shpstrct->ndime;jdime++){
	DDGRD[idime][jdime] = delta(idime,jdime)*(KK[1] + 2*KK[2]*temperature);
      }
    }

    // build stiffness

    double sumvec[2]; sumvec[0]=0;sumvec[1] = 0;
    for(int knode = 0;knode<shpstrct->numnodes;knode++){

      double BC[2];
      BC[0] = shpstrct->glbder[iinte][knode][0];
      BC[1] = shpstrct->glbder[iinte][knode][1];
      
      for(int kdime =0;kdime<shpstrct->ndime;kdime++){
	sumvec[kdime] += BC[kdime]*primal[knode][0];
      }
    }
    
    for (int inode = 0;inode<shpstrct->numnodes;inode++){
      double BA[2];
      BA[0] = shpstrct->glbder[iinte][inode][0];
      BA[1] = shpstrct->glbder[iinte][inode][1];
      
      for(int jnode = 0;jnode<shpstrct->numnodes;jnode++){
	double BB[2];
	BB[0] = shpstrct->glbder[iinte][jnode][0];
	BB[1] = shpstrct->glbder[iinte][jnode][1];

	for(int idime =0;idime<shpstrct->ndime;idime++){
	  for(int jdime = 0;jdime<shpstrct->ndime;jdime++){
	    estiff[inode][jnode] += BA[idime]*DD[idime][jdime]*BB[jdime]*weight*jacodet;
	  }
	}
	
	double shap = shpstrct->shape[iinte][jnode];
	for(int idime =0;idime<shpstrct->ndime;idime++){
	  for(int jdime=0;jdime<shpstrct->ndime;jdime++){
	    
	    estiff[inode][jnode] += BA[idime]*DDGRD[idime][jdime]*shap*sumvec[jdime]*jacodet*weight;
	  }
	}

      }//jnode

      // calculate residual
      
      for(int idime =0;idime<shpstrct->ndime;idime++){
	for(int jdime=0;jdime<shpstrct->ndime;jdime++){
	  erhs[inode] = erhs[inode] - BA[idime]*DD[idime][jdime]*sumvec[jdime]*jacodet*weight;
	  
	}
      }

    }//inode


  }// iinte

  deallocvec(KK);
  deallocmat(shpstrct->ndime,gradu);
  deallocmat(shpstrct->ndime,defgrad);
  deallocmat(shpstrct->ndime,rightcauchygreen);
  deallocmat(shpstrct->ndime,secondpiola);
  deallocmat(shpstrct->ndime,rightcauchygreeninvtpose);

}

//=======================================================================
void nonlinheat2d::calcfunc(){
  cout<<"calcfunc in nonlinheat2d()"<<endl;
}

//=======================================================================
void nonlinheat2d:: elemrhs(){

}

//=======================================================================
void nonlinheat2d:: calcgradient(){
  cout<<"calcgradient in nonlinheat2d()"<<endl;
}

//======================================================================
void nonlinheat2d::getinteg(){
  // get integration points
  gauss2d(shpstrct);
}
//====================================================================
void nonlinheat2d::getshape(){
  //cout<<"getshape in nonlinheat2d()"<<endl;
  shape2d(shpstrct);
}
//======================================================================
void nonlinheat2d::printeltype(){
  cout <<endl<<"This is NONLINHEAT2D"<<endl;
}
//=====================================================================
void nonlinheat2d::memalloc(){
  // material properties, connectivity, stiffnessmatrix, righthand side
  // ideqn,idnode,isbc
  // The best way is to inherit this..
  
}
//=====================================================================
void nonlinheat2d::memdealloc(){
  // dealloc memory

}
//=====================================================================
void nonlinheat2d::setcoord(double cval,int idime,int inode){
  // cval = coordinate value
  // inode = local nodenumber
  // idime = dimension
  // Again should, be in base class
  shpstrct->coord[inode][idime] = cval;
}
//====================================================================
void nonlinheat2d::setprop(double propval,int iprop,int inode){
  // cval = coordinate value
  // inode = local nodenumber
  // idime = dimension
  // Again should, be in base class
  prop[inode][iprop] = propval;
}
//===================================================================
void nonlinheat2d::setideqn(int inode,int idofn,long int eqnno){

  assert(inode<shpstrct->numnodes);  // 
 
 
  assert(idofn<dofnatnode[inode]);   // 
  


  ideqn[inode][idofn] = eqnno;
}
//==================================================================
void nonlinheat2d::setidnode(int inode,long int globalnode){
  assert(inode<shpstrct->numnodes);
  idnode[inode]= globalnode;
}
//=================================================================
void nonlinheat2d::setisbc(int inode,int idofn,int value){
  assert(inode<shpstrct->numnodes);
  assert(idofn<dofnatnode[inode]);
  
  isbc[inode][idofn] = value;

}
//==============================================================
void nonlinheat2d::setdirdata(double ddata,int inode,int idofn){
   assert(inode<shpstrct->numnodes);
   assert(idofn<dofnatnode[inode]);
   
   dirdata[inode][idofn] = ddata;
}
//==============================================================
void nonlinheat2d::settracdata(double tdata,int inode,int idofn){
  assert(inode<shpstrct->numnodes);
  assert(idofn<dofnatnode[inode]);

  tracdata[inode][idofn] = tdata;
}
//=============================================================

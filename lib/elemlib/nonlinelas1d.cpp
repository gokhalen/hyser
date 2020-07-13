/*======================================================================
  Nachiket Gokhale gokhalen@bu.edu
  One Dimensional Non-Linear Elasticity Element
  Builds consistent tangent and right hand side
  Builds dual rhs, functional, calculates gradient
=======================================================================*/

#ifndef HAVE_NONLINELAS1D_HPP
#include "nonlinelas1d.hpp"
#endif

#include <cassert>
#include <cstdlib>

//=========================================================================
//implement constructor
nonlinelas1d::nonlinelas1d(int* cinte,long int celemno){
  // have all information necessary. call constrctor of shapestrct
  // set nprop from inmesh!
  nprop = 1;
  mdofn = 1;
  nonlinflag = true;
  shpstrct = (shapestruct *)(new shapestruct(1,1,cinte,2,NONLINELAS1D,celemno));

  // allocate the rest after number of nodes has been set
  prop       = allocmat(shpstrct->numnodes,nprop,prop); 
  propvec    = allocvec(nprop,propvec);
  estiff     = allocmat(shpstrct->numnodes*mdofn,shpstrct->numnodes*mdofn,estiff);
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

  // allocate space for inversion data
  measdata  = allocmat(shpstrct->numnodes,mdofn,measdata);
  zeromat(shpstrct->numnodes,mdofn,measdata);

  // allocate dual field data
  dualfield = allocmat(shpstrct->numnodes,mdofn,dualfield);
  zeromat(shpstrct->numnodes,mdofn,dualfield);

  // allocate gradient matrix
  
  gradient = allocmat(shpstrct->numnodes,nprop,gradient);
  zeromat(shpstrct->numnodes,nprop,gradient);

    
  }
//==========================================================================
nonlinelas1d::~nonlinelas1d(){

  // the destructor

  deallocmat(shpstrct->numnodes,prop);
  deallocvec(propvec);
  deallocmat(shpstrct->numnodes*mdofn,estiff);
  deallocvec(dofnatnode);
  deallocmat(shpstrct->numnodes,primal);
  deallocmat(shpstrct->numnodes,dirdata);
  deallocmat(shpstrct->numnodes,tracdata);
  deallocvec(erhs);
  deallocmat(shpstrct->numnodes,ideqn);
  deallocvec(idnode);
  deallocmat(shpstrct->numnodes,isbc);

  //deallocate measdata
  deallocmat(shpstrct->numnodes,measdata);

  // deallocate dualfield
  deallocmat(shpstrct->numnodes,dualfield);

  // deallocate gradient 
  deallocmat(shpstrct->numnodes,gradient);
  
  // deallocate shpstrct
  // use the delete operation, do not call destructor directly..
  //shpstrct->~shapestruct();
  delete shpstrct;
}
  //========================================================================
// implement element stiffness routine
void nonlinelas1d::elemstiff(int which_stiffness){
  
  // declare conslib structure & allocate mem. 
  conselas1d *elasptr; 
  elasptr = (conselas1d *)malloc(sizeof(conselas1d));
  
  

  
  // zero element stiffness K_{A,B} = \int_{\Omega}(N_{A,x}k(x)N_{B,X})
  zeromat(shpstrct->numnodes*mdofn,shpstrct->numnodes*mdofn,estiff);
  zerovec(shpstrct->numnodes*mdofn,erhs);
  // loop over integration points
  for(int iinte = 0;iinte<shpstrct->tinte;iinte++){

    // get jacobian determinant and weight
    
    double jacodet = shpstrct->jdet[iinte];
    double weight  = shpstrct->wts[iinte];

    // interpolate stiffness at this integration point
    double stiff = 0;
    double K1 = 0;

    for(int iprop = 0;iprop<nprop;iprop++){
      // zero the iprop-th propvec
      propvec[iprop] = 0;
      for(int knode = 0;knode<shpstrct->numnodes;knode++){
	// assume: stress (one d) = K*U,x;
	// where K = K1 + K1*U,X;
	// interpolate K1;
	K1 +=  shpstrct->shape[iinte][knode]*prop[knode][0]; 
	
	propvec[iprop] += shpstrct->shape[iinte][knode]*prop[knode][0];
      }
    }

    // compute deformation gradient
    // F = grad(\phi) = grad(X + U) = 1 + U,X (in one d)
    
    double defgrad = 0; // deformation gradient  
    double gradu   = 0; // gradient of u
    double disp    = 0; // displacement
    
    for (int inode = 0;inode<shpstrct->numnodes;inode++){
      gradu +=  shpstrct->glbder[iinte][inode][0]*this->primal[inode][0];
      
      if(isbc[inode][0] == 1){
	  // we have hit dirichlet node
	 disp += dirdata[inode][0]*shpstrct->shape[iinte][inode];
	}
	else{
	  // we have hit primal node
	  disp += primal[inode][0]*shpstrct->shape[iinte][inode];
	}
    }

    // add the one
    defgrad = gradu + 1;
    
    assert(defgrad !=0);
    
    // calculate the complete stiffness: K1 + K1*U,X 

    //stiff = K1 + K1*gradu;

    //double first   = (gradu*K1)/defgrad;
    //double second  = -stiff*gradu/(defgrad*defgrad);
    //double third   = stiff/defgrad;

    
    // we will now set the conselas1d data structure
    elasptr->grdparm     = 1;
    elasptr->matdisp      = disp;
    elasptr->matstrain    = gradu;
    elasptr->stfparam     = &this->propvec[0];
    
    // we call the function to get consistent tangent contrib
    elas1d(this->conseqnno,elasptr);

    
    
    // evaluate the consistent (to what?) tangent
    
    for(int inode =0;inode<shpstrct->numnodes;inode++){     
      // picked the node corresp. to index A and its shape function

      double der1 = shpstrct->glbder[iinte][inode][0];
      
      for(int jnode=0;jnode<shpstrct->numnodes;jnode++){     
	// pick the node corresp. to index B and its shape function

	double der2    = shpstrct->glbder[iinte][jnode][0];
	
	//estiff[inode][jnode] += der1*( first + second + third )*der2*jacodet*weight;
	estiff[inode][jnode] +=der1*(elasptr->ctang)*der2*jacodet*weight;
      }
      
      //erhs[inode] -= (der1*stiff*gradu/defgrad)*jacodet*weight;
      erhs[inode] -= (der1*elasptr->stiffness*gradu/defgrad)*jacodet*weight;
    }
  }
  // dealloc elasptr
  free(elasptr);
}

//=======================================================================
void nonlinelas1d::calcfunc(){
  // calculate functional 

  // begin by zeroing functional defined in elem_base
  this->functional = 0;

  // loop over each integration point

  for(int iinte = 0;iinte<shpstrct->tinte;iinte++){
    // get jacobian determinant and weight
    double jacodet = shpstrct->jdet[iinte];
    double weight  = shpstrct->wts[iinte];


    // only one degree of freedom
    // at this integ pt we need to interp primal and meas data
	
      double prim = 0;
      double meas = 0;
      for (int inode = 0;inode<shpstrct->numnodes;inode++){
	if(isbc[inode][0] == 1){
	  // we have hit dirichlet node
	  prim += dirdata[inode][0]*shpstrct->shape[iinte][inode];
	}
	else{
	  // we have hit primal node
	  prim += primal[inode][0]*shpstrct->shape[iinte][inode];
	}

	// interp meas data
	meas  += measdata[inode][0]*shpstrct->shape[iinte][inode];
      }
	
      // now calculate functional
      
      this->functional += 0.5*(meas-prim)*(meas-prim)*jacodet*weight;
  }
  
}

//=======================================================================
void nonlinelas1d:: elemrhs(){

  // zero righthand side vector
  
  // calculate residual vector
  
  


}

//=======================================================================
void nonlinelas1d:: calcgradient(){

  // declare a pointer to a conselas1d structure & allocate memory
  
  conselas1d *elasptr; elasptr = (conselas1d *)malloc(sizeof(conselas1d));
  
  
  // loop over integration points
  // this->gradient has been zeroed in the constructor

  for(int iinte = 0;iinte<shpstrct->tinte;iinte++){
    double jacodet = shpstrct->jdet[iinte];
    double weight  = shpstrct->wts[iinte];
    // now pick the property field to calculate the gradient with respect to


    // at this integration point: interpolate W,X and U,X

    double dualder = 0; double primder = 0;
    double disp = 0; 
    for (int knode = 0; knode<shpstrct->numnodes ;knode++){

      // interp dualder
      dualder += this->dualfield[knode][0]*shpstrct->glbder[iinte][knode][0];
      
      // the primal and bc's are stored separately
      if(this->isbc[knode][0] == 1){ 
	// we have hit dirichlet node
	primder += this->dirdata[knode][0]*shpstrct->glbder[iinte][knode][0]; // calculates displacement derivative
	disp    += this->dirdata[knode][0]*shpstrct->shape[iinte][knode];     // calculates displacement
      }
      else{
	primder += this->primal[knode][0]*shpstrct->glbder[iinte][knode][0]; // calculates displacement derivative
	disp    += this->primal[knode][0]*shpstrct->shape[iinte][knode];     // calculates displacement
      }
      
    }


    // interpolate properties
    for(int iprop = 0;iprop<nprop;iprop++){
      // zero the iprop-th propvec
      propvec[iprop] = 0;
      for(int knode = 0;knode<shpstrct->numnodes;knode++){
	// assume: stress (one d) = K*U,x;
	// where K = K1 + K1*U,X;
	// interpolate K1;
	//K1 +=  shpstrct->shape[iinte][knode]*prop[knode][0]; 

	propvec[iprop] += shpstrct->shape[iinte][knode]*prop[knode][0];
      }
    }

    // using primder initialize defgrad

    double defgrad = (1+primder);
    
    elasptr->matdisp    = disp;
    elasptr->matstrain  = primder;
    elasptr->stfparam   = &this->propvec[0]; 

    for(int iprop = 0;iprop<this->nprop;iprop++){
      // currently this law has only one property and so things are fine
      // however, in the near future this call must be replaced by 
      // a call to a "conslib" routine which will compute and return
      // the appropriate product
      
      // set which parameter we are taking the gradient with respectto
      elasptr->grdparm = iprop;

      // we call this function to get gradcontrib value 
      elas1d(this->conseqnno,elasptr);
      
      for(int inode = 0;inode<shpstrct->numnodes;inode++){

	// the stiffness function in this code is 
	// k = K1(1 + U,X) 
	// wonder if this is stable etc. 
	
	//double gradshp = shpstrct->shape[iinte][inode]*(1+primder);
	double gradshp   = shpstrct->shape[iinte][inode];
	
	//this->gradient[inode][iprop] += dualder*(gradshp)*primder*(1.0/defgrad)*jacodet*weight;
	this->gradient[inode][iprop] += dualder*(gradshp)*(elasptr->gradcontrib)*primder*(1.0/defgrad)*jacodet*weight;

      }// inode

    } //iprop


  } //iinte

  // deallocate a pointer to the conselas1d structure
  free(elasptr);
}

//======================================================================
void nonlinelas1d::getinteg(){
  // get integration points
  gauss1d(shpstrct);
}
//====================================================================
void nonlinelas1d::getshape(){
  //cout<<"getshape in nonlinelas1d()"<<endl;
  shape1d(shpstrct);
}
//======================================================================
void nonlinelas1d::printeltype(){
  cout <<endl<<"This is NONLINELAS1D"<<endl;
}
//=====================================================================
void nonlinelas1d::memalloc(){
  // material properties, connectivity, stiffnessmatrix, righthand side
  // ideqn,idnode,isbc
  // The best way is to inherit this..
  
}
//=====================================================================
void nonlinelas1d::memdealloc(){
  // dealloc memory

}
//=====================================================================
void nonlinelas1d::setcoord(double cval,int idime,int inode){
  // cval = coordinate value
  // inode = local nodenumber
  // idime = dimension
  // Again should, be in base class
  shpstrct->coord[inode][idime] = cval;
}
//====================================================================
void nonlinelas1d::setprop(double propval,int iprop,int inode){
  // cval = coordinate value
  // inode = local nodenumber
  // idime = dimension
  // Again should, be in base class
  prop[inode][iprop] = propval;
}
//===================================================================
void nonlinelas1d::setideqn(int inode,int idofn,long int eqnno){

  assert(inode<shpstrct->numnodes);  // 
 
 
  assert(idofn<dofnatnode[inode]);   // 
  


  ideqn[inode][idofn] = eqnno;
}
//==================================================================
void nonlinelas1d::setidnode(int inode,long int globalnode){
  assert(inode<shpstrct->numnodes);
  idnode[inode]= globalnode;
}
//=================================================================
void nonlinelas1d::setisbc(int inode,int idofn,int value){
  assert(inode<shpstrct->numnodes);
  assert(idofn<dofnatnode[inode]);
  
  isbc[inode][idofn] = value;

}
//==============================================================
void nonlinelas1d::setdirdata(double ddata,int inode,int idofn){
   assert(inode<shpstrct->numnodes);
   assert(idofn<dofnatnode[inode]);
   
   dirdata[inode][idofn] = ddata;
}
//==============================================================
void nonlinelas1d::settracdata(double tdata,int inode,int idofn){
  assert(inode<shpstrct->numnodes);
  assert(idofn<dofnatnode[inode]);

  tracdata[inode][idofn] = tdata;
}
//=============================================================
void nonlinelas1d::setmeasdata(double mdata,int inode,int idofn){
  assert(inode<shpstrct->numnodes);
  assert(idofn<dofnatnode[inode]);
  
  measdata[inode][idofn] = mdata;

}
//===============================================================
void nonlinelas1d::calcdual(){
  // begin by zeroing rhsvector
  zerovec(shpstrct->numnodes*mdofn,erhs);

  // loop over each integration point
  for(int iinte = 0;iinte<shpstrct->tinte;iinte++){
    // get jacobian determinant and weight
    double jacodet = shpstrct->jdet[iinte];
    double weight  = shpstrct->wts[iinte];

    // only one degree of freedom
    // at this integ pt we need to interp primal and meas data
	
    double prim = 0;
    double meas = 0;
    for (int inode = 0;inode<shpstrct->numnodes;inode++){
      if(this->isbc[inode][0] == 1){
	// we have hit dirichlet node
	prim += dirdata[inode][0]*shpstrct->shape[iinte][inode];
      }
      else{
	// we have hit primal node
	prim += primal[inode][0]*shpstrct->shape[iinte][inode];
      }
      
      // interp meas data
      meas  += measdata[inode][0]*shpstrct->shape[iinte][inode];
    }//inode

    // the shape function
    for(int inode = 0;inode<shpstrct->numnodes;inode++){
      // this equation has sign mistake in prospectus
      this->erhs[inode] += -(prim-meas)*shpstrct->shape[iinte][inode]*jacodet*weight;
    }
  }
}
//============================================================================
void nonlinelas1d::setdualfield(double ddata,int inode,int idofn){
  assert(inode<shpstrct->numnodes);
  assert(idofn<dofnatnode[inode]);
  
  dualfield[inode][idofn] = ddata;
}

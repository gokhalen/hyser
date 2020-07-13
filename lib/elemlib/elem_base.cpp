/*========================================================================
  Nachiket Gokhale gokhalen@bu.edu: Feb 2005 elem_base.cpp 
  Base class for finite elements. Keep the common data and functions here.
  This defines virtual functions that other elements must implement
==========================================================================*/

#ifndef HAVE_ELEM_BASE_CPP
#define HAVE_ELEM_BASE_CPP
#endif

#include <math.h>

#ifndef HAVE_SHAPESTRUCT_HPP
#include "../datastruct/shapestruct.hpp"
#endif

#ifndef HAVE_ELEMLIB_HPP
#include "elemlib.hpp"
#endif


// FINITE ELEMENT STUFF
  
void elem_base::setmatprops() {}           // set material properties
void elem_base::setfemarrays(){}           // set femarrays
   
//VIRTUAL FUNCTIONS MUST BE IMPLEMENTED
//FORWARD PROBLEM
void elem_base::elemstiff(int dummy){}     // calculate elemental stiffness 
void elem_base::elemcmass(int dummy){} // consistent mass matrix
void elem_base::elemrhs(){}       // calculate right hand side
void elem_base::calcenh(){}       // calculate the enhanced parameters
void elem_base::eleminitstrn(){}  // elem init strn
void elem_base::homogenizehassanihinton(){} // the homogenization integral
void elem_base::computeaveprops(){}

//==========================================================================
//INVERSE PROBLEM 
void elem_base::addreggradient(){

  // loop over integration points

  double *propgrad;
  
  propgrad = allocvec(shpstrct->ndime,propgrad);

  for(int iinte = 0; iinte < shpstrct->tinte; iinte++){
    double jacodet = shpstrct->jdet[iinte];
    double weight  = shpstrct->wts[iinte];
    
    // loop over the properties
    for (int iprop = 0; iprop < this->nregprop; iprop++){
      
      // calculate the gradient of the properties
      // these regularizations do not work correctly when the properties are linked to each other
      // these regularizations do not consider opttype - same for regfunctional

      zerovec(shpstrct->ndime,propgrad);
      for(int idime = 0; idime <shpstrct->ndime ;    idime++){
	for(int inode = 0; inode < shpstrct->numnodes; inode++){
	  propgrad[idime] += shpstrct->glbder[iinte][inode][idime]*this->prop[inode][iprop];
	}
      }

      // calculate the iprop^th property at this integration point
      double propval = 0.0;
      for(int inode = 0; inode < shpstrct->numnodes; inode++){
	propval += shpstrct->shape[iinte][inode]*this->prop[inode][iprop];
      }
      
      int knodemax = shpstrct->numnodes;
      double gradshp;
      if(this->propdiscontinflag == 1){
	knodemax = 1;
      }
      // pick the node wrt which we need to get the gradient
      for (int knode = 0; knode < knodemax ; knode++){
	gradshp = shpstrct->shape[iinte][knode];
	if(this->propdiscontinflag == 1){
	  gradshp =1.0;
	}
	switch(this->which_regularization){
	  
	case L2NORM:
	  {
	    double alpha = this->regparm[iprop][0];
	    this->gradient[knode][iprop] += alpha*propval*gradshp*jacodet*weight;
	  }
	  break;
	  
	  // H1 SEMINORM 
	case H1SEMINORM:
	  {
	    double alpha = this->regparm[iprop][0];
	    for(int idime = 0; idime < shpstrct->ndime ; idime++){
	      for(int jdime = 0; jdime <shpstrct->ndime; jdime++){
		if(idime == jdime){
		  
		  this->gradient[knode][iprop] += alpha*propgrad[idime]*shpstrct->glbder[iinte][knode][idime]*jacodet*weight;
		}
	      }
	    }
	  }
	  break;

	case TOTALVARIATION:
	  {
	    double alpha = this->regparm[iprop][0];
	    double cc    = this->regparm[iprop][1];
	    for(int idime = 0; idime < shpstrct->ndime ; idime++){
	      for(int jdime = 0; jdime <shpstrct->ndime; jdime++){
		if(idime == jdime){
		  double dotprdprop = pow(vecnorm(this->nregprop,propgrad),2.0);
		  double oosqrt = (1.0/sqrt(dotprdprop + (cc*cc)));
		  this->gradient[knode][iprop] += 0.5*alpha*oosqrt*propgrad[jdime]*shpstrct->glbder[iinte][knode][jdime]*jacodet*weight;
		}
	      }//jime
	    }//idime
	    break;
	  }

	case PAULNORM:
	  {
	    double alpha = this->regparm[iprop][0];
	    double cc    = this->regparm[iprop][1];
	    for(int idime = 0; idime < shpstrct->ndime; idime++){
	      for (int jdime = 0; jdime < shpstrct->ndime; jdime++){
		if(idime == jdime){
		  double last = propgrad[idime]*shpstrct->glbder[iinte][knode][jdime]*jacodet*weight;
		  double dotprdprop = pow(vecnorm(this->nregprop,propgrad),2.0);
		  double deno = (dotprdprop + (cc*cc));
		  this->gradient[knode][iprop] += alpha*(1.0/(deno) - (dotprdprop)/((deno)*(deno)))*last;
		}
	      }
	    }
	  }

	case PENTAREG1:
	  {

		double alpha = this->regparm[iprop][0];
		double mumax = this->regparm[iprop][1];
		double mumin = this->regparm[iprop][2];
		
		//cout<<"elem_base.cpp"<<endl;
		//cout<<"alpha == "<<alpha<<endl;
		//cout<<"mumax == "<<mumax<<endl;
		//cout<<"mumin == "<<mumin<<endl;
		
		this->gradient[knode][iprop] += alpha*(propval - mumax )*(propval - mumin)*( (2.0*propval) - mumax - mumin )*shpstrct->shape[iinte][knode]*jacodet*weight;

	    
	  }
	  
	break;
	  
	}

      }
      
      
    }//iprop
  } //iinte

  deallocvec(propgrad);
}

//==========================================================================
void elem_base::calcdual(){
  // calculate dorcing for dual
  // calculate the forcing for the dualfield (adjoint field)
  // again, closely modelled on nonlinelas1d::calcdual;
  // the old code for the dual forcing (without the tensor T) is 
  // at the end of the file
  // begin by zeroing rhsvector
  zerovec(shpstrct->numnodes*mdofn,erhs);
  
  // loop over each integration point
  
  for(int iinte = 0;iinte<shpstrct->tinte;iinte++){
    double jacodet = shpstrct->jdet[iinte];
    double weight  = shpstrct->wts[iinte];
    double combweight = jacodet*weight*funcweight;
    // pick the node and dime/dofn (pdime) to calculate the forcing for
    // this determines where things go in the forcing vector

    for (int pdime = 0; pdime < shpstrct->ndime; pdime++){
      for (int pnode =0; pnode < shpstrct->numnodes; pnode++){
	// this gets the eqn number in the local forcing vector
	int peqn = mdofn*(pnode - 1 + 1) + pdime; 
	
	//pick the dimension j 
	if(which_functional == 0){
	  for (int jdime =0;jdime < shpstrct->ndime; jdime++){
	    double tprim  = 0.0;
	    double tmeas  = 0.0;
	    
	    // eval the difference(or, tprim & tmeas) at this intg pt
	    for (int jnode=0;jnode<shpstrct->numnodes ;jnode++){
	      tprim += Tprim[jnode][jdime]*shpstrct->shape[iinte][jnode];
	      tmeas += Tmeas[jnode][jdime]*shpstrct->shape[iinte][jnode];
	    }//jnode
	    
	    
	    // multiply with T_{jdime,pdime} and N_{pnode} and you have your forcing
	    //double combweight = jacodet*weight*funcweight;
	    double forc = -(tprim-tmeas)*functens[jdime][pdime]*shpstrct->shape[iinte][pnode]*combweight;
	    
	    this->erhs[peqn]  += forc;
	    
	  }// jdime++
	}

	// need to account for H1 part of functional
	
	if(which_functional == 1){
	  for(int idime = 0; idime < shpstrct->ndime; idime++){
	    for(int jdime = 0; jdime < shpstrct->ndime; jdime++){
	      double tprim = 0.0;
	      double tmeas = 0.0;
	      double forc  = 0.0;
	      for(int inode = 0; inode < shpstrct->numnodes; inode++){
		tprim += Tprim[inode][idime]*shpstrct->glbder[iinte][inode][jdime];
		tmeas += Tmeas[inode][idime]*shpstrct->glbder[iinte][inode][jdime];
	      }
	      
	      forc = -(tprim - tmeas)*functens[idime][pdime]*shpstrct->glbder[iinte][pnode][jdime]*combweight;
	      this->erhs[peqn] +=forc;
	      
	      
	    }
	  }
	} // which_functional ==1

	if ( ( which_functional == PENTANORM ) || ( which_functional == PENTANORM2 ) || (which_functional == PENTANORM3 ) ) {
	  if (( which_functional == PENTANORM2 ) || (which_functional == PENTANORM3 )){
	    cout<<"Warning: elem_base.cpp: PENTANORMS 2 and 3 ARE NOT CORRECT!-INVESTIAGE"<<endl;
	  }
	  
	  double forc = 0.0;

         // calculate lambda and mu at this integration point
	  double lambda = 0; double mu = 0;
	  for(int inode = 0;inode<shpstrct->numnodes;inode++){
	    lambda += shpstrct->shape[iinte][inode]*prop[inode][0];
	    mu     += shpstrct->shape[iinte][inode]*prop[inode][1];
	  }
	  
	  //cout<<"GETS TO PENTANORM  == "<<initstrntype;


	  for ( int pp = 0; pp < shpstrct->ndime; pp++){
	    for ( int qq = 0; qq < shpstrct->ndime; qq++){
	      for ( int jj = 0; jj < shpstrct->ndime; jj++){
		//for ( int kk = 0; kk < shpstrct->ndime;  kk++){
		//for (int ll = 0; ll < shpstrct->ndime; ll++){
		    int kk, ll;
		    if ( this->initstrntype == 1 ) { kk = 0; ll =0;} 
		    if ( this->initstrntype == 2 ) { kk = 1; ll =1;} 
		    if ( this->initstrntype == 3 ) { kk = 0; ll =1;} 
		    
		//cout<<"elem_base.cpp initstrntype == "<<initstrntype;

		    int ii = pdime;
		    double gpqkl  = 0.0;

		    
		    if ( this->which_functional == PENTANORM ){
		      if ( this->BASE_PMIJKL[pp][qq][kk][ll] > -5.0 ){
			gpqkl  = (this->BASE_EHIJKL[pp][qq][kk][ll] - this->BASE_PMIJKL[pp][qq][kk][ll]);
		      }
		    }

		    if ( this->BASE_PMIJKL[pp][qq][kk][ll] > 0.0 ) {

		      if ( this->which_functional == PENTANORM2 ){
			gpqkl  = (this->BASE_EHIJKL[pp][qq][kk][ll] - this->BASE_PMIJKL[pp][qq][kk][ll]);
			gpqkl  = gpqkl / (this->BASE_PMIJKL[pp][qq][kk][ll]*this->BASE_PMIJKL[pp][qq][kk][ll]);
		      }

		      if ( this->which_functional == PENTANORM3 ){
			gpqkl  = (this->BASE_EHIJKL[pp][qq][kk][ll] - this->BASE_PMIJKL[pp][qq][kk][ll]);
			gpqkl  = gpqkl / this->BASE_PMIJKL[pp][qq][kk][ll];
		      }
		    }

		    double microeijpq  = lambda*(delta(ii,jj)*delta(pp,qq)) + mu*( delta(ii,pp)*delta(jj,qq) + delta(jj,pp)*delta(ii,qq));
		    forc +=  -(gpqkl/(this->unitcellvolume))*microeijpq*shpstrct->glbder[iinte][pnode][jj]*jacodet*weight;

	      }
	    }
	  }
	  //cout <<"force in calcdual.cpp == "<<forc<<endl;
          this->erhs[peqn] +=forc;
	}
 
                   
 

      }//pdofn
    } // pdime

  }//iinte
}

//============================================================================
void elem_base::calcfunc(){

  // calculate the functional 
  // done with nonlinelas1d as a model
  // begin by zeroing the functional

  this->functional    = 0;
  this->regfunctional = 0;

  for(int iinte =0; iinte<shpstrct->tinte;iinte++){
    double jacodet = shpstrct->jdet[iinte];
    double weight  = shpstrct->wts[iinte];
    double combweight = funcweight*jacodet*weight;

  
    
    // compute the tmeas and tprim at this integration point
    if(this->which_functional == 0) {
           
      for (int idime = 0; idime < shpstrct->ndime; idime++){
	double tmeas = 0.0;
	double tprim = 0.0;
	
	for (int inode = 0; inode < shpstrct->numnodes; inode++){
	  tmeas += Tmeas[inode][idime]*shpstrct->shape[iinte][inode];
	  tprim += Tprim[inode][idime]*shpstrct->shape[iinte][inode];
	}
	
	
	this->functional += 0.5*(tmeas-tprim)*(tmeas-tprim)*combweight;

      }
      
    }

    if(this->which_functional==1){

	for(int idime = 0; idime < shpstrct->ndime; idime++){
	  for(int jdime = 0; jdime < shpstrct->ndime; jdime++){
	    
	    double tmeas = 0.0;
	    double tprim = 0.0;
	    // compute derivatives of quantities
	    for (int inode = 0; inode < shpstrct->numnodes; inode++){
	      tmeas += Tmeas[inode][idime]*shpstrct->glbder[iinte][inode][jdime];
	      tprim += Tprim[inode][idime]*shpstrct->glbder[iinte][inode][jdime];
	      
	    }
	    
	    this->functional += 0.5*(tmeas-tprim)*(tmeas-tprim)*combweight;
	    
	    //this->functional +=combweight;

	    }//idime 
	  }//jdime
    
    }
 

    if ( (this->which_functional == PENTANORM) || (this->which_functional == PENTANORM2) || (this->which_functional == PENTANORM3 ) ){

      for (int ii = 0; ii < shpstrct->ndime; ii++){
	for (int  jj = 0; jj < shpstrct->ndime; jj++){
	  for (int kk = 0; kk < shpstrct->ndime; kk++){
	    for (int ll = 0; ll < shpstrct->ndime; ll++){	   

	      this->BASE_GIJKL[ii][jj][kk][ll]    = 0.0;

	      if ( this->which_functional == PENTANORM ){
		if ( this->BASE_PMIJKL[ii][jj][kk][ll] > -5.0 ){
		  this->BASE_GIJKL[ii][jj][kk][ll] = (this->BASE_EHIJKL[ii][jj][kk][ll] - this->BASE_PMIJKL[ii][jj][kk][ll]);
		}
	      }
	      
	      if ( this->BASE_PMIJKL[ii][jj][kk][ll] > 0.0 ){
		if (this->which_functional == PENTANORM2){
		  this->BASE_GIJKL[ii][jj][kk][ll]    = (this->BASE_EHIJKL[ii][jj][kk][ll] - this->BASE_PMIJKL[ii][jj][kk][ll])/(this->BASE_PMIJKL[ii][jj][kk][ll]);;
		}
		else if (this->which_functional == PENTANORM3){
		  this->BASE_GIJKL[ii][jj][kk][ll] = (this->BASE_EHIJKL[ii][jj][kk][ll] - this->BASE_PMIJKL[ii][jj][kk][ll]);
		  this->BASE_GIJKL[ii][jj][kk][ll] = this->BASE_GIJKL[ii][jj][kk][ll] / sqrt(this->BASE_PMIJKL[ii][jj][kk][ll]);
		}
	      }

	      
	      this->functional += ( 0.5 /  ( this->unitcellvolume * this->nfields)) * BASE_GIJKL[ii][jj][kk][ll] * BASE_GIJKL[ii][jj][kk][ll] * jacodet * weight;
	      


	    }
	  }
	}
      }
	    
    } // pentanorm




    
    // the regularization to the functional
    // H1-semi norm regularization is used. 

    for (int iprop = 0; iprop < this->nregprop ; iprop++){
      double regh1semifunclocal = 0.0;
      double regl2funclocal         = 0.0;

      // calculate the iprop^th property at this integration point
      double propval = 0.0;
      for(int inode = 0; inode < shpstrct->numnodes; inode++){
	propval += shpstrct->shape[iinte][inode]*this->prop[inode][iprop];
      }
      
	for(int idime = 0; idime < shpstrct->ndime; idime++ ){
	  double gradprop      = 0.0;
	  
	  // interpolate the derivative of prop 
	  for (int inode = 0; inode < shpstrct->numnodes ; inode++){
	    // opttype does not need to be considered here  - taken care of in calcfuncgrad.cpp
	    gradprop += shpstrct->glbder[iinte][inode][idime]*this->prop[inode][iprop];
	  }
	  regh1semifunclocal += gradprop*gradprop;
	}//idime
	
	{
	  double proplocal = 0.0;
	  // now interpolate the material property at this point
	  for (int inode = 0; inode < shpstrct->numnodes ; inode++){
	    proplocal += shpstrct->shape[iinte][inode]*this->prop[inode][iprop];
	  }
	  regl2funclocal     += proplocal*proplocal;
	}
	
	
	switch(this->which_regularization){

	case L2NORM:
	  {
	    double alpha = this->regparm[iprop][0];
	    functional += 0.5*alpha*regl2funclocal*jacodet*weight;
	    this->regfunctional += 0.5*alpha*regl2funclocal*jacodet*weight;
	    break;
	  }

	case H1SEMINORM:
	  {
	    double alpha = this->regparm[iprop][0];
	    functional += 0.5*alpha*regh1semifunclocal*jacodet*weight;
	    this->regfunctional += 0.5*alpha*regh1semifunclocal*jacodet*weight;
	    break;
	  }
	  
	case TOTALVARIATION:
	  {
	    double alpha = this->regparm[iprop][0];
	    double cc    = this->regparm[iprop][1];

	    assert(cc > 0); 

	    double oosqrt = sqrt(regh1semifunclocal + (cc*cc));
	    double common = 0.5*alpha*oosqrt*jacodet*weight;
	    
	    functional += common;
	    this->regfunctional += common;
	    break;
	  }

	case PAULNORM:  
	  {
	    double alpha = this->regparm[iprop][0];
	    double cc    = this->regparm[iprop][1];
	    double denominator = (regh1semifunclocal + (cc*cc));
	    double common      = 0.5*alpha*(regh1semifunclocal/denominator)*jacodet*weight;
	    functional += common;
	    this->regfunctional  += common;
	  }
          
        case PENTAREG1:
	  {
		double alpha = this->regparm[iprop][0];
		double mumax = this->regparm[iprop][1];
		double mumin = this->regparm[iprop][2];

		double common        = 0.5*alpha*pow(( mumax - propval),2.0)*pow(mumin - propval,2.0)*jacodet*weight;
		functional          += common;
		this->regfunctional += common;

					   
	  }
	  
	  break;
	  
	}
	
    }// iprop
    
  } // iinte

  //cout<<__FILE__<<" "<<__LINE__<<" functional == "<<this->functional<<endl;
}      // calculate functional
//=============================================================================
void elem_base::calcgradient(){}  // calculate gradient

//INTEGRATION POINTS          // gets integration points 
void elem_base::getinteg(){}
   
//SHAPE FUNCTIONS
void elem_base::getshape(){}
//=====================================================================
void elem_base::setcoord(double cval,int idime,int inode){
  // cval = coordinate value
  // inode = local nodenumber
  // idime = dimension
  // Again should, be in base class
  shpstrct->coord[inode][idime] = cval;
}
//====================================================================
void elem_base::setkvec(double kval,int idime){
  kvec[idime] = kval;
}
//====================================================================
void elem_base::setwdisp(double ww){
  this->wdisp = ww;
}
//====================================================================
void elem_base::setprop(double propval,int iprop,int inode){
  // cval = coordinate value
  // inode = local nodenumber
  // idime = dimension
  // Again should, be in base class
  prop[inode][iprop] = propval;
}
//===================================================================
void elem_base::setdensity(double propval, int inode){
  density[inode] = propval;
}
//===================================================================
void elem_base::setideqn(int inode,int idofn,long int eqnno){

  assert(inode<shpstrct->numnodes);  // 
  assert(idofn<dofnatnode[inode]);   // 
  ideqn[inode][idofn] = eqnno;       //

}
//==================================================================
void elem_base::setidnode(int inode,long int globalnode){
  assert(inode<shpstrct->numnodes);
  idnode[inode]= globalnode;
}
//=================================================================
void elem_base::setisbc(int inode,int idofn,int value){
  assert(inode<shpstrct->numnodes);
  assert(idofn<dofnatnode[inode]);
  
  isbc[inode][idofn] = value;

}
//==============================================================
void elem_base::setdirdata(double ddata,int inode,int idofn){
   assert(inode<shpstrct->numnodes);
   assert(idofn<dofnatnode[inode]);
   
   dirdata[inode][idofn] = ddata;
}
//==============================================================
void elem_base::settracdata(double tdata,int inode,int idofn){
  assert(inode<shpstrct->numnodes);
  assert(idofn<dofnatnode[inode]);

  tracdata[inode][idofn] = tdata;
}
//=============================================================
void elem_base::setmeasdata(double mdata,int inode,int idofn){
  assert(inode<shpstrct->numnodes);
  assert(idofn<dofnatnode[inode]);
  
  measdata[inode][idofn] = mdata;

}
//=============================================================
void elem_base::setdualfield(double ddata,int inode,int idofn){  
  assert(inode<shpstrct->numnodes);
  assert(idofn<dofnatnode[inode]);
  
  dualfield[inode][idofn] = ddata;
}


elem_base::elem_base(){

  // constructors are called in this sequence 
  

  // initialize rhsnorm
  this->rhsnorm = 0;
  
  
}
   
elem_base::~elem_base(){   // the virtual destructor is extended not overridden
  //memdealloc();
}

   
//==========================
void elem_base::viewid(){
  for(int inode = 0;inode < shpstrct->numnodes;inode++){
    //  cout<<endl<<__FILE__<<" global node "<<idnode[inode]<<" " <<shpstrct->numnodes;
  }
}

//=========================
void elem_base::setprimal(double vval,int inode,int idofn){
  
  this->primal[inode][idofn] = vval;
}
//=============================
void elem_base::setlastupdate(double vval,int inode,int idofn){
  this->lastupdate[inode][idofn] = vval;
}

//=========================
void elem_base::memalloc(){
  // allocate the regparm array
  // can be done in the base class once the virtual function stuff
  // is figured out
  //cout<<"initializing regparm from elem_base.cpp"<<endl;
  kvec          = allocvec(shpstrct->ndime,kvec);
  regparm       = allocmat(this->nregprop,this->mregparm,this->regparm);
  BASE_EHIJKL   = alloc4D(shpstrct->ndime, shpstrct->ndime, shpstrct->ndime, shpstrct->ndime, BASE_EHIJKL);
  BASE_EHAVIJKL = alloc4D(shpstrct->ndime, shpstrct->ndime, shpstrct->ndime, shpstrct->ndime, BASE_EHAVIJKL);

  BASE_CIJKL    = alloc4D(shpstrct->ndime, shpstrct->ndime, shpstrct->ndime, shpstrct->ndime, BASE_CIJKL);
  BASE_PMIJKL   = alloc4D(shpstrct->ndime, shpstrct->ndime, shpstrct->ndime, shpstrct->ndime, BASE_PMIJKL);
  BASE_GIJKL    = alloc4D(shpstrct->ndime, shpstrct->ndime, shpstrct->ndime, shpstrct->ndime, BASE_GIJKL);
  BASE_GHATIJKL = alloc4D(shpstrct->ndime, shpstrct->ndime, shpstrct->ndime, shpstrct->ndime, BASE_GHATIJKL);
  BASE_PZHKIJ   = alloc3D(shpstrct->ndime, shpstrct->ndime, shpstrct->ndime,BASE_PZHKIJ);
  BASE_PZHAVKIJ = alloc3D(shpstrct->ndime, shpstrct->ndime, shpstrct->ndime,BASE_PZHAVKIJ);

  BASE_KHIJ     = allocmat(shpstrct->ndime, shpstrct->ndime, BASE_KHIJ);
  BASE_KHAVIJ   = allocmat(shpstrct->ndime, shpstrct->ndime, BASE_KHAVIJ);

}

//==========================
void elem_base::memdealloc(){
  deallocvec(this->kvec);
  deallocmat(this->nregprop,this->regparm);
  dealloc4D(shpstrct->ndime, shpstrct->ndime, shpstrct->ndime, shpstrct->ndime, BASE_EHIJKL);
  dealloc4D(shpstrct->ndime, shpstrct->ndime, shpstrct->ndime, shpstrct->ndime, BASE_EHAVIJKL);
  dealloc4D(shpstrct->ndime, shpstrct->ndime, shpstrct->ndime, shpstrct->ndime, BASE_CIJKL);
  dealloc4D(shpstrct->ndime, shpstrct->ndime, shpstrct->ndime, shpstrct->ndime, BASE_PMIJKL);
  dealloc4D(shpstrct->ndime, shpstrct->ndime, shpstrct->ndime, shpstrct->ndime, BASE_GIJKL);
  dealloc4D(shpstrct->ndime, shpstrct->ndime, shpstrct->ndime, shpstrct->ndime, BASE_GHATIJKL);
  dealloc3D(shpstrct->ndime, shpstrct->ndime, BASE_PZHKIJ);
  dealloc3D(shpstrct->ndime, shpstrct->ndime, BASE_PZHAVKIJ);
  deallocmat(shpstrct->ndime, BASE_KHIJ);
  deallocmat(shpstrct->ndime, BASE_KHAVIJ);

}
//===========================
void elem_base::calculate_Tmeas_Tprim(){
  // zero Tmeas and Tprim
    zeromat(shpstrct->numnodes,shpstrct->ndime,Tprim);
    zeromat(shpstrct->numnodes,shpstrct->ndime,Tmeas);
    
    // put the values of Tmeas and Tprim in them
    // by multiplying functens and prim and meas
    // this can/should be done by a cute matrix multiplying
    // routine but we'll leave it for now Jul 18 2005
 
      for(int inode  = 0;inode < shpstrct->numnodes;inode++){
	for(int idime = 0;idime < shpstrct->ndime ; idime++){

	  // multiply functens and prim and meas
	  for(int jdime = 0;jdime < shpstrct->ndime; jdime++){

	    double first = functens[idime][jdime];

	    Tmeas[inode][idime] += first*this->measdata[inode][jdime];
	    Tprim[inode][idime] += first*this->primal[inode][jdime];
	  }
	
	}
      }
}
//===========================================

void elem_base::calcnormmeas(){

  // calculate Tmeas Tprim must be called beforehand
  // this is called by the localizer

  this->normmeasdata = 0;

   for (int iinte = 0; iinte < shpstrct->tinte; iinte++){
    double jacodet = shpstrct->jdet[iinte];
    double weight  = shpstrct->wts[iinte];
    double combweight = funcweight*jacodet*weight;

    if(this->which_functional == 0){
      for (int idime = 0; idime <shpstrct->ndime; idime++){
	
	double tmeas  = 0.0;
	
	// interpolate idime component of tmeas at this point
	for (int jnode=0;jnode<shpstrct->numnodes ;jnode++){
	  tmeas += Tmeas[jnode][idime]*shpstrct->shape[iinte][jnode];
	}//jnode
	
	this->normmeasdata += tmeas*tmeas*combweight;
      }
    } // L2 norm
    if (this->which_functional == 1){
       for(int idime = 0; idime < shpstrct->ndime; idime++){
	 for(int jdime = 0; jdime < shpstrct->ndime; jdime++){
	   double tmeas = 0.0;
	   // calculate derivative at this integration point
	   for (int inode = 0; inode < shpstrct->numnodes; inode++){
	     tmeas += Tmeas[inode][idime]*shpstrct->glbder[iinte][inode][jdime];     	   }
	   this->normmeasdata += (tmeas)*(tmeas)*combweight;
      
	 } // jdime
       } //idime
    } // H1-norm

  } //iinte


}  
//=============================================================
void elem_base::computevolumemass(){
  // homogenize the density
  this->elemmass   = 0.0;
  this->elemvolume = 0.0;

  
  for(int iinte=0; iinte < shpstrct->tinte; iinte++){
      double jacodet = shpstrct->jdet[iinte];
      double weight  = shpstrct->wts[iinte];
      
      double rhotmp = 0.0;
      
      for(int inode=0; inode<shpstrct->numnodes; inode++){
	rhotmp +=  this->density[inode]*shpstrct->shape[iinte][inode];
      }

      

      this->elemvolume += jacodet*weight;
      this->elemmass   += rhotmp*jacodet*weight;
  }
  
  //this->elemmass = -1.0;
  //this->elemvolume = -1.0;
}

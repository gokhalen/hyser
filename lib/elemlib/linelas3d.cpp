/*======================================================================
  Nachiket Gokhale gokhalen@bu.edu
  Three dimensional Linear Elasticity Element
  Builds stiffness matrix and right hand side
  
  No assumptions on special structure of the stiffness matrices are made.
  This enables easy extension to the B_{bar} formulation
========================================================================*/
#ifndef HAVE_LINELAS3D_HPP
#include "linelas3d.hpp"
#endif

#include <cassert>
#include <cstdlib>

using namespace std;


// implement constrctor 

linelas3d::linelas3d(int *cinte,int cmregparm, long int celemno){
  
  //cout<<"Calling Constructor for Linelas3d"<<endl;
  
  nprop          = 2;
  mdofn          = 3;
  this->mregparm = cmregparm;
  nonlinflag     = false;

  shpstrct = new shapestruct(3,3,cinte,8,LINELAS3D,celemno);

  // allocate the rest after the number of nodes has been set
  prop       = allocmat(shpstrct->numnodes,nprop,prop); 
  estiff     = allocmat(shpstrct->numnodes*mdofn,shpstrct->numnodes*mdofn,estiff);
  dofnatnode = allocvec(shpstrct->numnodes,dofnatnode);

  // set the dofnatnode array, one degree of freedom at each node
  for(int inode = 0;inode<shpstrct->numnodes;inode++){
    dofnatnode[inode] = mdofn;
  }

  // primal
  primal  = allocmat(shpstrct->numnodes,mdofn,primal);
  zeromat(shpstrct->numnodes,mdofn,primal);
  
  lastupdate = allocmat(shpstrct->numnodes,mdofn,lastupdate);
  zeromat(shpstrct->numnodes,mdofn,lastupdate);
  
  initstrn = allocvec(3*(shpstrct->ndime-1) , initstrn);


  // allocate dirichlet data and traction data vectors
  dirdata  = allocmat(shpstrct->numnodes,mdofn,dirdata);
  zeromat(shpstrct->numnodes,mdofn,dirdata);

  tracdata = allocmat(shpstrct->numnodes,mdofn,tracdata);
  zeromat(shpstrct->numnodes,mdofn,tracdata);

  // allocate right hand side vector
  erhs = allocvec(shpstrct->numnodes*mdofn,erhs);
  
  // allocate finite element arrays
  // (node num, dofn) -> glb eqn no
  ideqn  = allocmat(shpstrct->numnodes,mdofn,ideqn); 

  // local -> global node num.
  idnode = allocvec(shpstrct->numnodes,idnode);

  // checks if bc. (1 if dirichlet, -1 if traction 0 otherwise
  isbc   = allocmat(shpstrct->numnodes,mdofn,isbc); 
  zeromat(shpstrct->numnodes,mdofn,isbc);

  // allocate space for inversion of data
  measdata  = allocmat(shpstrct->numnodes,mdofn,measdata);
  zeromat(shpstrct->numnodes,mdofn,measdata);

  // allocate dual field data
  dualfield = allocmat(shpstrct->numnodes,mdofn,dualfield);
  zeromat(shpstrct->numnodes,mdofn,dualfield);

 // allocate the whole data - used in gradient calculation only
 // HACK for now - change with a constructor argument 
  primalwhole = alloc3D(shpstrct->numnodes,mdofn,6,primalwhole);
  zero3D(shpstrct->numnodes,mdofn,6,primalwhole);
  
  dualwhole = alloc3D(shpstrct->numnodes,mdofn,6,dualwhole);
  zero3D(shpstrct->numnodes,mdofn,6,dualwhole);

  // allocate gradient matrix
  
  gradient = allocmat(shpstrct->numnodes,nprop,gradient);
  zeromat(shpstrct->numnodes,nprop,gradient);

  functens = allocmat(shpstrct->ndime,shpstrct->ndime,functens);

  Tprim    = allocmat(shpstrct->numnodes,shpstrct->ndime,Tprim);
  Tmeas    = allocmat(shpstrct->numnodes,shpstrct->ndime,Tmeas);

  zeromat(shpstrct->numnodes,shpstrct->ndime,Tprim);
  zeromat(shpstrct->numnodes,shpstrct->ndime,Tmeas);  

  //cout<<"Calling Memalloc from linelas3d.cpp"<<endl;
  
  memalloc();
  
}
//==================================================================
linelas3d::~linelas3d(){

  deallocmat(shpstrct->numnodes,prop);
  deallocmat(shpstrct->numnodes*mdofn,estiff);
  deallocvec(dofnatnode);
  deallocmat(shpstrct->numnodes,primal);
  deallocmat(shpstrct->numnodes,lastupdate);
  deallocvec(initstrn);

  deallocmat(shpstrct->numnodes,dirdata);
  deallocmat(shpstrct->numnodes,tracdata);
  deallocvec(erhs);
  deallocmat(shpstrct->numnodes,ideqn);
  deallocvec(idnode);
  deallocmat(shpstrct->numnodes,isbc);
  deallocmat(shpstrct->numnodes,measdata);
  deallocmat(shpstrct->numnodes,dualfield);
  dealloc3D(shpstrct->numnodes,mdofn,primalwhole);
  dealloc3D(shpstrct->numnodes,mdofn,dualwhole);
  deallocmat(shpstrct->numnodes,gradient);

  deallocmat(shpstrct->ndime,functens);
  deallocmat(shpstrct->numnodes,Tprim);
  deallocmat(shpstrct->numnodes,Tmeas);
  
  memdealloc();
  delete shpstrct;
}
//================================================================
void linelas3d::elemstiff(int which_stiffness){
  // zero element stiffness K_{A,B} = \int_{\Omega}(N_{A,x}k(x)N_{B,X})
  zeromat(shpstrct->numnodes*mdofn,shpstrct->numnodes*mdofn,estiff);

  //loop over integration points
  
  for(int iinte=0;iinte<shpstrct->tinte;iinte++){
    double jacodet = shpstrct->jdet[iinte];
    double weight  = shpstrct->wts[iinte];
    
    //cout<<__FILE__<<"WTJAC =="<<shpstrct->jdet[0]*shpstrct->wts[0];

    double lambda = 0.0; double mu = 0.0;

    // interpolate \lambda and \mu at this point
    for(int inode = 0;inode<shpstrct->numnodes;inode++){
      lambda += shpstrct->shape[iinte][inode]*prop[inode][0];
      mu     += shpstrct->shape[iinte][inode]*prop[inode][1];
    }// close interpolation of \lambda & \mu


    // cout<<__FILE__<<" lambda mu "<<lambda<<" "<<mu;
    double DD[6][6]; 
    
    // zero DD
    for(int ii= 0;ii<6;ii++){
      for(int jj=0;jj<6;jj++){
	DD[ii][jj] = 0.0;
      }
    }
     
    DD[0][0] = lambda  + 2*mu;
    DD[0][1] = lambda;
    DD[0][2] = lambda;
    DD[1][0] = lambda;
    DD[1][1] = lambda + 2*mu;
    DD[1][2] = lambda;
    DD[2][0] = lambda;
    DD[2][1] = lambda;
    DD[2][2] = lambda + 2*mu;
    DD[3][3] = mu;
    DD[4][4] = mu;
    DD[5][5] = mu;
    
    
    for(int inode = 0; inode<shpstrct->numnodes; inode++){
      // set up B_{a} and initialize it
      double BA[6][3];
      //zero it 
      for(int ii = 0;ii<6;ii++){
	for(int jj =0;jj<3;jj++){
	  BA[ii][jj] = 0;
	}
      }
      
      BA[0][0] = shpstrct->glbder[iinte][inode][0];
      BA[1][1] = shpstrct->glbder[iinte][inode][1];
      BA[2][2] = shpstrct->glbder[iinte][inode][2];

      BA[3][1] = shpstrct->glbder[iinte][inode][2];
      BA[3][2] = shpstrct->glbder[iinte][inode][1];

      BA[4][0] = shpstrct->glbder[iinte][inode][2];
      BA[4][2] = shpstrct->glbder[iinte][inode][0];

      BA[5][0] = shpstrct->glbder[iinte][inode][1];
      BA[5][1] = shpstrct->glbder[iinte][inode][0];

      for(int jnode =0; jnode<shpstrct->numnodes; jnode++){
	double BB[6][3];
	// quite unnecessary can be moved out of loop
	for(int ii = 0;ii<6;ii++){
	  for(int jj =0;jj<3;jj++){
	    BB[ii][jj] = 0;
	  }
	}
	
	BB[0][0] = shpstrct->glbder[iinte][jnode][0];
	BB[1][1] = shpstrct->glbder[iinte][jnode][1];
	BB[2][2] = shpstrct->glbder[iinte][jnode][2];
	
	BB[3][1] = shpstrct->glbder[iinte][jnode][2];
	BB[3][2] = shpstrct->glbder[iinte][jnode][1];

	BB[4][0] = shpstrct->glbder[iinte][jnode][2];
	BB[4][2] = shpstrct->glbder[iinte][jnode][0];

	BB[5][0] = shpstrct->glbder[iinte][jnode][1];
	BB[5][1] = shpstrct->glbder[iinte][jnode][0];
	

	
	
	for(int idofn = 0;idofn<dofnatnode[inode];idofn++){
	  // need the (inode  - 1 + 1) because of zero based indexing
	  int ieqn = (dofnatnode[inode]*(inode-1 + 1) + idofn);
	  for(int jdofn = 0;jdofn<dofnatnode[jnode];jdofn++){
	    int jeqn = (dofnatnode[jnode]*(jnode-1 + 1) + jdofn);

	    // approx 83% of time is spent the following lines of code
	    for(int ii=0;ii<6;ii++){
	      for(int jj=0;jj<6;jj++){
		
 		estiff[ieqn][jeqn] += BA[ii][idofn]*DD[ii][jj]*BB[jj][jdofn]*jacodet*weight;
		  //rmat_ba(kk,idofn)*rmat_d(kk,ll)*
		  //$                          rmat_bb(ll,jdofn)*wtjac
	      }
	    }
		 
	    
	  }// jdofn
	} // idofn
	

      } // closes jnode
      
    }//closes inode
    
  } // closes iinte

  //cout<<endl;

  //    for(int ii = 0;ii<8;ii++){
  //   for(int jj =0;jj<8;jj++){
  //cout<<__FILE__<<" stiff "<<estiff[0][0]<<" ";
  //   }
  //  cout<<endl;
  //  }
  //fflush(stdin);
  //exit(1);

}//closes elemstiff
//====================================================================
//void linelas3d::calcfunc(){
//  cout<<"calcfunc in linelas3d()"<<endl;
//}
//====================================================================
void linelas3d::elemrhs(){
  // zero right hand side vector
  zerovec(shpstrct->numnodes*mdofn,erhs);
  
  // elemstiff should be call from outside before calling elemrhs
  // dirichlet data: take product of matrix with right hand side
  
  for(int inode = 0;inode<shpstrct->numnodes;inode++){
    for(int idofn = 0;idofn<dofnatnode[inode];idofn++){
      
      // zero based indexing. therefore, -1 + 1 
      int ieqn = (inode-1 + 1)*(dofnatnode[inode]) + idofn;
      //take matrix-vector product

      for(int jnode=0;jnode<shpstrct->numnodes;jnode++){
	for(int jdofn=0;jdofn<dofnatnode[jnode];jdofn++){
	  // zero based indexing. therefore, -1 + 1
	  int jeqn = (jnode-1+1)*dofnatnode[jnode] + jdofn;
	  erhs[ieqn] -= estiff[ieqn][jeqn]*dirdata[jnode][jdofn];
	} // jnode
      }   // jdofn
    }     // idofn
  }       // inode

  //for(int jj =0;jj< 24 ;jj++){
  // cout<<"linelas3d.cpp erhs "<<erhs[jj]<<" "<<endl;
  //}
  
}//elemrhs()
//====================================================================
void linelas3d::calcgradient(){
  cout<<"calcgradient in "<<__FILE__;
}
//====================================================================
void linelas3d::getinteg(){
  gauss3d(shpstrct);
}
//====================================================================
void linelas3d::getshape(){
  shape3d(shpstrct);
}
//====================================================================
void linelas3d::printeltype(){
  cout<<endl<<"This is LINELAS3D"<<endl;
}
//====================================================================
//void linelas3d::memalloc(){

  // inherited
//}
//====================================================================
//void linelas3d::memdealloc(){
  // dealloc memory
//}
//===================================================================
//void linelas3d::setcoord(double cval,int idime,int inode){
  // cval = coordinate value
  // inode = local nodenumber
  // idime = dimension
  // Again should, be in base class
  //shpstrct->coord[inode][idime] = cval;
//}
//====================================================================
//void linelas3d::setprop(double propval,int iprop,int inode){
  // cval = coordinate value
  // inode = local nodenumber
  // idime = dimension
  // Again should, be in base class
  //prop[inode][iprop] = propval;
//}
//===================================================================
//void linelas3d::setideqn(int inode,int idofn,long int eqnno){

//  assert(inode<shpstrct->numnodes);  // 
//
//  assert(idofn<dofnatnode[inode]);   // 
//  
//
//
//  ideqn[inode][idofn] = eqnno;
//}
//==================================================================
//void linelas3d::setidnode(int inode,long int globalnode){
//  assert(inode<shpstrct->numnodes);
//  idnode[inode]= globalnode;
//}
//=================================================================
// void linelas3d::setisbc(int inode,int idofn,int value){
//  assert(inode<shpstrct->numnodes);
//  assert(idofn<dofnatnode[inode]);
//  
//  isbc[inode][idofn] = value;
//
// }
//==============================================================
// void linelas3d::setdirdata(double ddata,int inode,int idofn){
//   assert(inode<shpstrct->numnodes);
//   assert(idofn<dofnatnode[inode]);
//   
//   dirdata[inode][idofn] = ddata;
// }
//==============================================================
// void linelas3d::settracdata(double tdata,int inode,int idofn){
//  assert(inode<shpstrct->numnodes);
//  assert(idofn<dofnatnode[inode]);
//
//  tracdata[inode][idofn] = tdata;
// }
//=============================================================
void linelas3d::eleminitstrn(){

  //` Hassani and Hinton show that the initial strain forcing is given by 
  //  f = \int B^{T} D \epsilon_{0}
  // where \epsilon_{0} is the initial strain

  
  // set the initial strain forcing if the initstrn type is >0 
  // otherwise this reduces to zero forcing

  zerovec(shpstrct->numnodes*mdofn,erhs);
  zerovec(3*(shpstrct->ndime-1),initstrn);
  //for(int ii=0;ii<6;ii++){cout<<"in linelas3d.cpp inistrn[ "<<ii<<"]=="<<initstrn[ii]<<endl;  }

  if ( ( initstrntype <= 3*(shpstrct->ndime-1) ) && (initstrntype > 0)) {
      // note: the components of the strainvector are (e11,e22,e33,2e23,2e13,2e12)
      // the following line is correct; the shear strain vector should not be set to 0.5
    initstrn[initstrntype-1]=1.0;

  }

  //for(int ii=0;ii<6;ii++){cout<<"in linelas3d.cpp inistrn[ "<<ii<<"]=="<<initstrn[ii]<<endl;}
  //for ( int ii = 0; ii < 24; ii++){cout<<__FILE__<<__LINE__<<" erhs, ieqn "<< erhs[ii]<<" "<<ii<<endl;}

  // at each integration point calculate B^{T}D
  for ( int iinte = 0; iinte < shpstrct->tinte; iinte++){

    double jacodet = shpstrct->jdet[iinte];
    double weight  = shpstrct->wts[iinte];

    // set up the D array at this point
    // first interpolate the properties
    double lambda = 0.0; double mu = 0.0;
    
    for(int inode = 0;inode<shpstrct->numnodes;inode++){
      lambda += shpstrct->shape[iinte][inode]*prop[inode][0];
      mu     += shpstrct->shape[iinte][inode]*prop[inode][1];
    }// close interpolation of \lambda & \mu

    double DD[6][6];
    
    // zero DD
    for(int ii= 0;ii<6;ii++){
      for(int jj=0;jj<6;jj++){
	DD[ii][jj] = 0.0;
      }
    }
     
    DD[0][0] = lambda  + 2*mu;
    DD[0][1] = lambda;
    DD[0][2] = lambda;
    DD[1][0] = lambda;
    DD[1][1] = lambda + 2*mu;
    DD[1][2] = lambda;
    DD[2][0] = lambda;
    DD[2][1] = lambda;
    DD[2][2] = lambda + 2*mu;
    DD[3][3] = mu;
    DD[4][4] = mu;
    DD[5][5] = mu;
    

    for ( int inode = 0; inode < shpstrct->numnodes; inode++){
      // set up  the B array at this point 
      double BA[6][3];
      //zero it 
      for(int ii = 0;ii<6;ii++){
	for(int jj =0;jj<3;jj++){
	  BA[ii][jj] = 0.0;
	}
      }
      
      BA[0][0] = shpstrct->glbder[iinte][inode][0];
      BA[1][1] = shpstrct->glbder[iinte][inode][1];
      BA[2][2] = shpstrct->glbder[iinte][inode][2];

      BA[3][1] = shpstrct->glbder[iinte][inode][2];
      BA[3][2] = shpstrct->glbder[iinte][inode][1];

      BA[4][0] = shpstrct->glbder[iinte][inode][2];
      BA[4][2] = shpstrct->glbder[iinte][inode][0];

      BA[5][0] = shpstrct->glbder[iinte][inode][1];
      BA[5][1] = shpstrct->glbder[iinte][inode][0];

      for ( int idofn = 0; idofn < dofnatnode[inode]; idofn++){
	int ieqn = (dofnatnode[inode]*(inode) + idofn);

	for (int ii = 0; ii < 6; ii++){
	  for(int jj = 0; jj < 6; jj++){
	    erhs[ieqn] += BA[ii][idofn]*DD[ii][jj]*initstrn[jj]*jacodet*weight;
	  }
	}

      } // idofn
    } // inode

  }   //iinte
  
  //for ( int ii = 0; ii < 24; ii++){cout<<__FILE__<<__LINE__<<" erhs, ieqn "<< erhs[ii]<<" "<<ii<<endl;}

}
//=============================================================
void linelas3d::homogenizehassanihinton(){
  zero4D(shpstrct->ndime,shpstrct->ndime,shpstrct->ndime,shpstrct->ndime,this->BASE_EHIJKL);

  // interpolate the properties lambda and mu
  int kk = 0;
  int ll = 0;  
  if (  !( ( initstrntype <= 3*(shpstrct->ndime-1) ) && (initstrntype >=0)  ) ){
    cout<<endl<<"INVALID VALUE FOR initstrntype in linelas3d.cpp, initstrntype = "<<initstrntype<<endl;
    exit(1);
  }
  
  // strain vector ordering (11,22,33,23,13,12)
  if ( initstrntype == 1 ){ kk = 0; ll = 0; }
  if ( initstrntype == 2 ){ kk = 1; ll = 1; }
  if ( initstrntype == 3 ){ kk = 2; ll = 2; }
  if ( initstrntype == 4 ){ kk = 1; ll = 2; }
  if ( initstrntype == 5 ){ kk = 0; ll = 2; }
  if ( initstrntype == 6 ){ kk = 0; ll = 1; }

  for( int iinte = 0; iinte < shpstrct->tinte; iinte++){
    double jacodet = shpstrct->jdet[iinte];
    double weight  = shpstrct->wts[iinte];
    
    // calculate lambda and mu at this integration point
    double lambda = 0; double mu = 0;
    for(int inode = 0;inode<shpstrct->numnodes;inode++){
      lambda += shpstrct->shape[iinte][inode]*prop[inode][0];
      mu     += shpstrct->shape[iinte][inode]*prop[inode][1];
    }

    // cout<<__FILE__<<__LINE__<<" "<<lambda<<" "<<mu<<endl;
    
    for (int ii = 0; ii < shpstrct->ndime; ii++){
      for (int jj = 0; jj < shpstrct->ndime; jj++){
	
	// calculate the homogenized tensor
	double tempfourth1 = ( delta(ii,jj)*delta(kk,ll) );
	double tempfourth2 = ( delta(ii,kk)*delta(jj,ll) + delta(jj,kk)*delta(ii,ll));
	
	// we've multiplied by the jacodet here
	double firstcontrib  = ( lambda*tempfourth1 + mu*tempfourth2 )*jacodet*weight;

	double secondcontrib = 0.0;
	for (int pp = 0; pp < shpstrct->ndime; pp++){
	  for (int qq = 0; qq < shpstrct->ndime; qq++){
	    // the derivative of pth component of khi^{kl} wrt q 
	    double khi_kl_pq = 0.0;

	    for ( int inode = 0; inode < shpstrct->numnodes; inode++){
	      // the kl index implicitly choses the primal -- that's how it was localized
	      khi_kl_pq += shpstrct->glbder[iinte][inode][qq]*this->primal[inode][pp];
	    } //inode 
	    double eijpq   = lambda*(delta(ii,jj)*delta(pp,qq)) + mu*( delta(ii,pp)*delta(jj,qq) + delta(jj,pp)*delta(ii,qq));
	    secondcontrib += eijpq*khi_kl_pq*jacodet*weight;

	  } //qq
	}//pp

	BASE_EHIJKL[ii][jj][kk][ll]     += (firstcontrib-secondcontrib)/(this->unitcellvolume);
	
	if ( initstrntype > 3 ) {
	  // make sure we calculate for the switched indices as well
	  BASE_EHIJKL[ii][jj][ll][kk] += (firstcontrib-secondcontrib)/(this->unitcellvolume);     
	}
      }//jj
    }//ii
  }//iinte
  

    
  
}

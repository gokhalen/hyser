#ifndef HAVE_SHAPESTRUCT_HPP
#include "shapestruct.hpp"
#endif

#include <new>
using namespace std;

//====================================================================================
// implementation of shapealloc
 void shapestruct::shapealloc(int *cinte){
   
   
   // ninte     = allocvec(ndime,ninte);
   
   tinte = 1;
   
   for(int ii=0;ii<selfdime;ii++){ninte[ii]=cinte[ii];tinte*=ninte[ii];} // calculate total integration points, set integration points
   
   jdet      = allocvec (tinte,jdet);
   jacomat   = alloc3D  (tinte,ndime,ndime,jacomat);
   jacoinv   = alloc3D  (tinte,ndime,ndime,jacoinv);
   coord     = allocmat (numnodes,ndime,coord);   // allocate memory for **coord (numnodes,ndime); 
   intpts    = allocmat (tinte,ndime,intpts);    // need tinte here, constructor provides it 
   glbcrdint = allocmat (tinte,ndime,glbcrdint);
   shape     = allocmat (tinte,numnodes,glbcrdint);
   glbder    = alloc3D  (tinte,numnodes,ndime,glbder);
   natder    = alloc3D  (tinte,numnodes,ndime,natder);
   wts       = allocvec (tinte,wts);

   // wilson's shape functions
   wishape  = allocmat(tinte,ndime,wishape);
   winatder = alloc3D(tinte,ndime,ndime,winatder);
   wiglbder = alloc3D(tinte,ndime,ndime,wiglbder);
   

} // shapealloc ends


//====================================================================================
//implementation of default constructor
shapestruct::shapestruct(){
}
   

//===================================================================================
//implementation of the second constructor
shapestruct::shapestruct(int cndime,int cselfdime,int *cinte,int cnumnodes,int celtype,long int celemno){

  ndime    = cndime;
  selfdime = cselfdime;
  ninte    = allocvec(ndime,ninte); // allocate memory for ninte 

  for(int idime=0;idime<ndime;idime++){
    ninte[idime]    = cinte[idime]; //ninte has cinte's value
  }

   numnodes = cnumnodes;
   eltype   = celtype;
   elemno   = celemno;
   
   shapealloc(cinte);
  
}
//====================================================================================
// implementation of the destructor
shapestruct::~shapestruct(){
  
 // dellocate memory for all the things allocated in shapealloc
  deallocvec (ninte);
  deallocvec (jdet);
  dealloc3D  (tinte,ndime,jacomat);
  dealloc3D  (tinte,ndime,jacoinv);
  deallocmat (numnodes,coord); 
  deallocmat (tinte,intpts);
  deallocmat (tinte,glbcrdint);
  deallocmat (tinte, shape);
  dealloc3D  (tinte,numnodes,glbder);
  dealloc3D  (tinte,numnodes,natder);
  deallocvec (wts);

  deallocmat(tinte,wishape);
  dealloc3D (tinte,ndime,winatder);
  dealloc3D (tinte,ndime,wiglbder);
   
}



//==================================================================================
//implementation of viewgauss function: views the data set by the quadrature routines
void shapestruct:: viewgauss(){
 
  cout<<endl; //next line before diagnostic output
  for(int iinte=0;iinte<tinte;iinte++){
    cout<<" integration point number "<<iinte<<" weight "<<wts[iinte]<<" location ";
    for(int idime =0;idime<ndime;idime++){
      cout<<"\t"<<intpts[iinte][idime]<<" ";
    }
    cout<<endl;
  }
  


}

//====================================================================================
//implementation of viewshape function: views the data set by the shape function lib

void shapestruct::viewshape(){
  // view the jacodian determinant
  cout<<"\nThe jacobian determinant == ";
  for(int iinte = 0;iinte < tinte;iinte++){
    cout<<" "<<jdet[iinte]<<" ";
  }
  cout<<endl;

  //view the shape functions
  for(int inode=0;inode<numnodes;inode++){
    for(int iinte=0;iinte<tinte;iinte++){
      cout<<"shape function of node "<<inode<<" at integ. point "<<iinte<<" is "<<shape[iinte][inode]<<endl;
    }
  }

  //view shape function natural & global derivatives
  for(int idime = 0;idime<ndime;idime++){
    for(int inode =0;inode <numnodes;inode++){
      for(int iinte = 0;iinte<tinte;iinte++){
	cout<<"glbder "<<glbder[iinte][inode][idime]<<" natder "<<natder[iinte][inode][idime]<<" iinte == "<<iinte<<" inode "<<inode<<" idime "<<idime<<endl;
      }
    }
  }
  
  //global locations of integration points
  for(int iinte = 0;iinte < tinte; iinte++){
    cout << "global coord of inte point "<<iinte;
    for(int idime = 0;idime< ndime; idime++){
      cout<<" "<<glbcrdint[iinte][idime];
    }
    cout<<endl;
  }
	
  cout<< " Wilsons' Shape Functions" <<endl;

  // wilsons's shape functions
    //view the shape functions
  for(int inode=0;inode<2;inode++){
    for(int iinte=0;iinte<tinte;iinte++){
      cout<<"shape function of node "<<inode<<" at integ. point "<<iinte<<" is "<<wishape[iinte][inode]<<endl;
    }
  }

  //view shape function natural & global derivatives
  for(int idime = 0;idime<ndime;idime++){
    for(int inode =0;inode <2;inode++){
      for(int iinte = 0;iinte<tinte;iinte++){
	cout<<"glbder "<<wiglbder[iinte][inode][idime]<<" natder "<<winatder[iinte][inode][idime]<<" iinte == "<<iinte<<" inode "<<inode<<" idime "<<idime<<endl;
      }
    }
  }
  
}

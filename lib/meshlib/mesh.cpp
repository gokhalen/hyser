/*****************************************************************************
* Nachiket Gokhale gokhalen@bu.edu 
* Implementations of the methods for the mesh class. 
* This class includes the following functions 
*
* 1) Reading in the mesh/part of the mesh on this processor 
* 2) Allocating memory 
* 3) Making equation numbers, and the fillin 
* 4) It also contains a CSR Matrix object             
*  
* This class/or a similar class will eventually hold PETSc DAs and objects for
* 1) Coordinates
* 2) Connectivity
* 3) Material properties
* 4) Global to local mappings
******************************************************************************/

#ifndef HAVE_MESH_HPP
#include "mesh.hpp"
#endif

#ifndef HAVE_MATUTILS_HPP
#include "../utillib/matutils.hpp"
#endif

#ifndef HAVE_GLOBALS_H
#include "../../headers/globals.h"
#endif

#include <math.h>
#include <iomanip>

// implement constructor 
mesh::mesh(){
  // data initaliztion

  this->rhsnorm    = 0;
  this->functional = 0;
  
  // this is because of the hack in make_fillin, which used localizer.cpp, which uses
  // mesh->which_field

  this->which_field      = 0;  
  this->which_disperpt   = 0;
  this->which_disperseg  = 0;
  this->inverseiteration = 1;
  this->iloadits         = 0;
  this->inewtonits       = 0;
  this->localizelast     = false;
  this->init_primal_goodguess_flag=false;
  this->propdiscontinflag = 0;

  read_mesh();    // read amd alllocate memory
  make_ien();     // ien(element num, local node num) = glb node number
  make_eqn_num(); // makes equation num, starting at zero
}
//============================================================================
// implement the destructor
mesh::~mesh(){
  deallocate_arrays();
}
//============================================================================
void mesh::accumulatenormmeas(double measnorm){
  // this calculates the L2 norm(s) of the measured displacement field(s)
  // and writes norms(s) to the file measnorm{$which_field$}.out 
  // normmeasdata is zeroed after allocation

  this->normmeasdata[this->which_field] += measnorm;
}
//=============================================================================
void mesh::applybcstosol(double **solrhsarg){
  int ifield = this->which_field;
  for(long int inode = 0; inode < this->num_nodes ; inode++){
    for(long int idofn = 0; idofn < this->mdofn ; idofn++){
      
      if(isdirichlet(inode,idofn,ifield)){
	
	// this stores the solution
	this->solvec[inode][idofn][ifield] = getdirichlet(inode,idofn,iloadits,ifield);
	  //this->dirich[this->isbc[inode][idofn][ifield]-1][iloadits][ifield];
	// this stores the last update to the solution
	this->lastupdate[inode][idofn][ifield] = 0.0;
      }
      else{
	if(this->linearflag==0){
	  // update \b{U}
	  this->solvec[inode][idofn][ifield] += solrhsarg[this->id[inode][idofn][ifield]-1][ifield];
	  // store \Delta\b{U} needed for computing enhanced parameters
	  this->lastupdate[inode][idofn][ifield] = solrhsarg[this->id[inode][idofn][ifield]-1][ifield]; 
	}
	else{
	  this->solvec[inode][idofn][ifield] = solrhsarg[this->id[inode][idofn][ifield]-1][ifield];
	}
	
      }
    }
  }
}
//=============================================================================
void mesh::read_mesh(){
  /**
     read_mesh is a repackaging of read_mesh.c 
     This has to be changed at some point of time
     
  */   
  FILE *fp;

  volatile long int litemp,dirctr,tracctr, perctr;
  int isdirbc,istracbc;
  int pointctr; int pointbcflag;

  volatile double dtemp,dvalue;

  long int startelem,endelem;
  int contin,elemnodes,elemtype,elemopt;

  char fbuffer[MAX_BUFFER];
  // open input data file
  fp = fopen("data.in","r"); 
  
  // check if file correctly opened
  if(!fp){printf("Error opening data.in  in mesh.cpp\n");exit(1);} 

  //printf("%s",fbuffer); // read and  print title of problem, useless
  fgets(fbuffer,MAX_BUFFER,fp); 
  
  fscanf(fp,"%d %ld %ld %d %d %d %d %d %d %ld %ld %ld %ld %d %d %d %d %d\n",&nfields,&num_elem,&num_nodes,&ndime,&mnode,&mdofn,&nprop,&mregparm, &ndispersionseg, &dirdofn,&tracdofn,&perdofn, &pdofn,&band,&loadits,&newtonits,&intelligentnewton,&conseqnno); /*read in data*/

  this->nregprop = 2;  // only lambda and mu are regularized for inverse problems
  this->physics = PHYSICS_ELAS;
  if (( this->mdofn == 3 ) && ( this->nprop == 8 )){
    this->physics = PHYSICS_PIEZO;
  }

  this->allocmisc();

  // read: ndirich
  fgets(fbuffer,MAX_BUFFER,fp); 
  for (int ifield = 0; ifield < nfields; ifield++){
    fscanf(fp,"%ld %ld %ld %ld %ld %ld %ld %ld\n",&ndirich[ifield],&ndirdofn[ifield], &ntrac[ifield],&ntracdofn[ifield],&nperiodic[ifield],&nperdofn[ifield],&npforc[ifield],&npforcdofn[ifield]);
    printf(" ndirich == %d ndirdofn= = %d\n",ndirich[ifield], ndirdofn[ifield]);
    printf(" ntrac == %d ntracdofn = %d\n",ntrac[ifield], ntracdofn[ifield]);
    printf(" nperiodic == %d nperdofn = %d\n",nperiodic[ifield], nperdofn[ifield]);
  }

  
 

  // read options: 
  fgets(fbuffer,MAX_BUFFER,fp);//printf("should be options %s",fbuffer);// throw away "options:fdgrad"
  fscanf(fp,"%d %lf %d %d %d %d\n",&fdgradflag,&fdstep,&taoiterates,&ctangflag,&primalguessflg,&firstguessflag);
  
  fgets(fbuffer,MAX_BUFFER,fp); // throw away linear flag etc 
  fscanf(fp,"%d %d %d\n",&linearflag,&printflag,&propdiscontinflag);
  printf("%d,%d %d\n",linearflag,printflag,propdiscontinflag);
  
  if (propdiscontinflag == 1){
    num_prop_nodes = num_elem;
  }
  else{
    num_prop_nodes = num_nodes;
  }

  fgets(fbuffer,MAX_BUFFER,fp); // throw away solvetype
  fscanf(fp,"%d %d %ld \n",&solvetype, &which_mass, &num_eigen);

  this->computenumeqns();  // depends on solvetype

  if ( this->issolvedispersion()) {
    cout<<"Dispersion K->W Seg = "<<ndispersionseg<<endl;
    maxdispersionpts = 0;
    for (int iseg =0; iseg < ndispersionseg ; iseg++){
      
      if ( this->solvetype == SOLVE_TYPE_DISPERSIONK2W ){
	for(int ipt = 0; ipt < 2; ipt++){
	  for(int idime=0;idime < ndime; idime++){
	    double temp;
	    fscanf(fp,"%lf",&temp);
	    ibzcoord[iseg][ipt][idime] = temp;
	  }
	}
	fscanf(fp,"%d\n",&ndispersionpts[iseg]);
	
	if ( ndispersionpts[iseg] > maxdispersionpts ){maxdispersionpts = ndispersionpts[iseg];}
	if ( ndispersionpts[iseg] < 2 ){ cout<<"Segment "<<iseg<<" has number of dispersion pts less than 2 ....Exiting."<<endl; exit(1);}

	cout<<"coordibz start =="<<ibzcoord[iseg][0][0]<<" "<<ibzcoord[iseg][0][1]<<endl;
	cout<<"coordibz end   =="<<ibzcoord[iseg][1][0]<<" "<<ibzcoord[iseg][1][1]<<endl;
	cout<<"ndispersionpts=="<<ndispersionpts[iseg]<<endl;
      } // dispersion k->w
      
      if ( this->solvetype == SOLVE_TYPE_DISPERSIONW2K){
	fscanf(fp,"%lf %lf %lf %lf %lf %lf %d\n",&w2karray[iseg][0],&w2karray[iseg][1],&w2karray[iseg][2],&w2karray[iseg][3],&w2karray[iseg][4],&w2karray[iseg][5],&ndispersionpts[iseg]);
	
	{	  
	  double kx = w2karray[iseg][4]; double ky = w2karray[iseg][5]; 
	  double knorm    = sqrt((kx*kx) + (ky*ky));
	  double residual =  knorm - 1.0;
	  if ( fabs(residual - knorm) > 1.0e-6  ){
	    cout<<endl<<"Renormalizing wave-vector direction"<<endl;
	    kx = kx / knorm;  ky = ky / knorm;
	    w2karray[iseg][4] = kx;  w2karray[iseg][5] = ky;
	      
	  }
	}
	
	cout<<"fmin== "<<w2karray[iseg][0]<<" fmax== "<<w2karray[iseg][1]<<" alpha== "<<w2karray[iseg][2]<<" beta == "<<w2karray[iseg][3]<<" kx = "<<w2karray[iseg][4]<<" ky = "<<w2karray[iseg][5]<<" ndispersionpts= "<<ndispersionpts[iseg]<<endl;
	if ( ndispersionpts[iseg] > maxdispersionpts ){maxdispersionpts = ndispersionpts[iseg];}
	if ( ndispersionpts[iseg] < 2 ){ cout<<"Segment "<<iseg<<" has number of dispersion pts less than 2 ....Exiting."<<endl; exit(1);}
      }
      
    }// iseg
  fscanf(fp,"\n");
  }

  // read the optimization flag ( 0=(all params are free), i=ith parameter is the variable and everything else deps on it)
  fgets(fbuffer,MAX_BUFFER,fp);printf("should be opttype %s",fbuffer);
  fscanf(fp,"%d %lf\n",&opttype, &poissons);

  if ( poissons == 0.5 ){
    cout<<"Poisson's ratio = 0.5 not allowed in opttype"<<endl;
    exit(1);
  }
			  

  this->mutolam = (2*poissons/(1-2*poissons));
  this->lamtomu = 1.0/this->mutolam;
  
  // read the movie flag
  fgets(fbuffer,MAX_BUFFER,fp);  // throw away "options: movie"
  fscanf(fp,"%d \n",&movieflag); 

  //exit(1);
  fgets(fbuffer,MAX_BUFFER,fp);//printf("%s",fbuffer);
  allocate_arrays();   //allocate coord and other FEM arrays*/

  // can be called after arrays are allocated
  if ( this->solvetype == SOLVE_TYPE_DISPERSIONK2W) {this->computedisperkvecs();}
  if ( this->solvetype == SOLVE_TYPE_DISPERSIONW2K) {this->computedisperwvecs();}
  
  //read element type arrays
  do{
    fscanf(fp,"%ld %ld %d %d %d %d\n",&startelem,&endelem,&elemtype,&elemnodes,&elemopt,&contin);
    printf("%ld %ld %d %d %d\n",startelem,endelem,elemtype,elemnodes,contin);

    assert(startelem <= endelem);
    assert(endelem   <= num_elem);

    for(long int ielem = startelem;ielem<=endelem;ielem++){
      lmref[ielem-1][0] = elemtype; //element type
      lmref[ielem-1][1] = elemnodes; // element nodes
      lmref[ielem-1][2] = elemopt;
      //printf("lmref %ld %d %d \n",ielem-1,this->lmref[ielem-1][0],this->lmref[ielem-1][1]);
    }
    
  }while(contin ==1);
  
  

  // read the element arrays
  fgets(fbuffer,MAX_BUFFER,fp);
  fscanf(fp,"%lf\n",&this->unitcellvolume);
  if ( this->unitcellvolume < 0 ){
    printf("\n ... Negative unit cell volume ... exiting\n");
    fflush(stdout);
    exit(1);
  }
  printf("-------------------\n");
  printf("Volume is %lf\n",this->unitcellvolume);
  

  fgets(fbuffer,MAX_BUFFER,fp);printf("%s",fbuffer);

  //read coordinate field
  int ictr = 0;
  for(long int inode = 0;inode<this->num_nodes;inode++){
    fscanf(fp,"%ld",&litemp);
    //printf("%ld\n ",litemp);
    for(int idime=0;idime<ndime;idime++){
      fscanf(fp,"%lf",&dtemp);
      // this is a hack - there seem to be mysterious issues with the first format 
      // apparent loss of synchronization for large input sizes
      // coord[litemp-1][idime] = dtemp;
      coord[litemp-1][idime] = dtemp;
      if ( ictr == 0 ){
	coordmin[idime] = dtemp;
	coordmax[idime] = dtemp;
      }
      else{
	if ( coordmin[idime] > dtemp ){ coordmin[idime] = dtemp;}
	if ( coordmax[idime] < dtemp ){ coordmax[idime] = dtemp;}
      }
    }
    fscanf(fp,"\n");
    ictr++;
  }
  printf("\n");

  // write the coordinates
  this->writecoord();
  double tempvolume=1.0;
  for(int idime = 0; idime< ndime; idime++){
    tempvolume *=(coordmax[idime] -coordmin[idime]);
    meshdim[idime]=coordmax[idime] - coordmin[idime];
  } 

  this->unitcellvolume = tempvolume;
  
  printf("\n-----------------------------------------------------------\n");
  printf("Overriding unit Cell Volume with internally computed volume\n");
  printf("-----------------------------------------------------------\n");

  //get "conn"
  cout<<"Last point="<<coord[this->num_nodes-1][0]<<" " <<coord[this->num_nodes-1][1]<<endl;
  fgets(fbuffer,MAX_BUFFER,fp);printf("%s",fbuffer);
  //exit(1);
  for(long int ielem = 0;ielem<num_elem;ielem++){
    fscanf(fp,"%ld",&litemp); //printf("%ld ",litemp);
    
    for(int inode = 0;inode<lmref[litemp-1][1];inode++){
      long int nodetemp;
      fscanf(fp,"%ld",&nodetemp);
      if ( nodetemp > this->num_nodes ){
	printf("Error in Conn table node number greater than num_nodes \n");
	exit(1);
      }
      conn[litemp-1][inode]=nodetemp;
      //printf("%ld ",conn[litemp-1][inode]);
      this->elemsattached[nodetemp-1] +=1;
    }

    // for each element increment the node counter by 1 

    fscanf(fp,"%ld \n",&litemp);//printf("\n");
  }

  // find max connectivity
  long int maxattachnodeid  = 0;
  long int maxelemsattached = 0;
  for ( long int inode = 0; inode < this->num_nodes ; inode++){
    if ( elemsattached[inode] > maxelemsattached ){
      maxattachnodeid  = inode+1;
      maxelemsattached = elemsattached[inode];
    }
  }
  
  // 4 nodes per elements, worst case 3 independent 1 shared.
  int approxband = (maxelemsattached)*(this->mnode)*this->mdofn + 3;
  
  cout<<endl<<"Maximum elems attached = "<< maxelemsattached << " at node "<<maxattachnodeid<<" approximate bandwidth  = "<<approxband<<endl;
  
  if ( this->band < approxband ){
    cout<<" ------------- WARNING ------------- "<<endl;
    cout<<"Bandwidth in data.in is "<< band <<" - must be greater than approximate bandwidth "<<approxband<<endl;
    cout<<" ------------- WARNING ------------- "<<endl;
  }


  // get the prop field
  fgets(fbuffer,MAX_BUFFER,fp);//printf("%s",fbuffer);
  for(long int inode = 0;inode<num_prop_nodes;inode++){
    fscanf(fp,"%ld",&litemp);//printf("%ld ",litemp);
    // get the property type -
    // 1 - (\lambda, \mu) 
    // 2 - (E, \nu) 
    int iproptype;  double propfactor; 
    fscanf(fp,"%d",&iproptype);

    // lame parameters, do nothing
    if ( iproptype == 1) {
      for(int iprop = 0;iprop<nprop;iprop++){
	fscanf(fp,"%lf",&prop[litemp-1][iprop]);
	//printf(" %lf ",prop[litemp-1][iprop]);
      }
      //printf("\n");
      //double rho;
      //fscanf(fp,"%lf ",&rho);
      double rho=prop[litemp-1][2];
      density[litemp-1][0]=rho;
    }
    
    if ( iproptype == 2){
      double EE, nu,rho;
      fscanf(fp,"%lf %lf %lf",&EE,&nu, &rho);
      if ( nu == 0.5 ) {
	cout<<"Poisson's ratio = 0.5 in element/node = "<<inode<<" exiting.."<<endl;
	exit(1);
      }
      prop[litemp-1][0]=(EE*nu)/((1.0 + nu)*(1.0 -2.0*nu));
      prop[litemp-1][1]=EE/(2.0*(1.0 + nu));
      prop[litemp-1][2]=rho;
      density[litemp-1][0]=rho;
    }
    fscanf(fp,"\n");//printf("\n");

  }
  
  // get the constraints on the property values
  fgets(fbuffer, MAX_BUFFER,fp);
  for(int inode = 0; inode<num_prop_nodes;inode++){
    int tempnode;
    fscanf(fp,"%d\n",&tempnode);

    for(int iprop = 0; iprop < this->nregprop ; iprop++){
      for(int ibound = 0; ibound < 2; ibound++){
	
	fscanf(fp,"%lf",&propcons[inode][iprop][ibound]);

      }//ibound
      fscanf(fp,"\n"); //move to next line
    }//iprop
  }//inode

  
  // the dirichlet boundary conditions*/


  getbcs(isbc,dirich,ndirich,BC_TYPE_DIRICHLET,fp);
  
  // the initial strain forcing
  {
    int itemp;
    fgets(fbuffer,MAX_BUFFER,fp); printf("%s\n",fbuffer);
    for ( int ifield = 0; ifield < nfields; ifield++){
      fscanf(fp,"%d\n",&initstrn[ifield]);
      printf("initstrn[%d]=%d",ifield,initstrn[ifield]);
    }
    // scan the zero
    fscanf(fp,"%d\n",&itemp);
    printf("\n");
  }

  getbcs(isbc,trac,ntrac,BC_TYPE_TRACTION,fp);
  getbcs(isperbc,perbc,nperiodic,BC_TYPE_PERIODIC,fp);
  getbcs(ispointbc,pforc,npforc,BC_TYPE_POINTFORC,fp);
  fgets(fbuffer,MAX_BUFFER,fp);printf("%s",fbuffer); // ninte
  fscanf(fp,"%d\n",&ninte);printf("%d\n",this->ninte);

  
  // get the field title "goodguess"
  fgets(fbuffer,MAX_BUFFER,fp);
  {
    int itemp;
    for(int ifield = 0; ifield < this->nfields; ifield++){
      for (int inode = 0; inode < this->num_nodes; inode++){
	for (int idofn = 0; idofn < this->mdofn; idofn++){
	  fscanf(fp,"%lf",&this->gooddispguess[inode][idofn][ifield]);
	}
	fscanf(fp,"\n");	  
      }
      // get the terminator '0'
      fscanf(fp,"%ld\n",&itemp);
    }
  }


  // get the field title regfunctionaltype
  fgets(fbuffer,MAX_BUFFER,fp);
  {
    fscanf(fp,"%d %d",&this->functionaltype,&this->regfunctionaltype);
    fscanf(fp,"\n");
  }

  // get the field title regparm
  fgets(fbuffer,MAX_BUFFER,fp);
  {
    for(int iprop = 0; iprop < this->nregprop ;iprop++){
      for(int iparam = 0; iparam < this->mregparm ; iparam++){
	fscanf(fp,"%lf",&this->regparm[iprop][iparam]);
      }
      fscanf(fp,"\n");
    }
    fscanf(fp,"\n");
  }


  // get the tensor T which is used to modify the simple L_2 norm 
  // of (u-u_meas)
  // get the field title lintnsr
  // the tensor is stored in a row format
  // A tensor 
  //       [1  0   2]
  //  T  = [0  0   0]
  //       [0  0   3]
  //  is inputted (for lack of a better word) as 
  //
  //       1  0  2  0  0  0 0 0 3 
 
  fgets(fbuffer,MAX_BUFFER,fp);
  {
    for(int idime = 0; idime < this->ndime; idime++){
      for(int jdime = 0; jdime < this->ndime; jdime++){
	fscanf(fp,"%lf",&functens[idime][jdime]);
      }
    }
    fscanf(fp,"\n");
  }

  // get the weight for the functional 
  // no getting around that
  fgets(fbuffer,MAX_BUFFER,fp);
  {
    for(int ifield = 0; ifield  < this->nfields; ifield++){
      fscanf(fp,"%lf",&fieldfuncweight[ifield]);
    }
    fscanf(fp,"\n");
  }

  // get the target pentamode properties
  fgets(fbuffer,MAX_BUFFER,fp);printf("%s",fbuffer);
  {
    double D11,D12,D22,D66;
    fscanf(fp,"%lf %lf %lf %lf\n",&D11,&D12,&D22,&D66);
    printf("D11=%lf,D12=%lf,D22=%lf,D66=%lf\n",D11,D12,D22,D66);

    // set the properties - the minus one is ignored in the functional, dual, gradient evaluation
    PMIJKL[0][0][0][0] = D11;
    PMIJKL[0][0][0][1] = -10.0;
    PMIJKL[0][0][1][0] = -10.0;
    PMIJKL[0][0][1][1] = D12;

    PMIJKL[0][1][0][0] = -10.0;
    PMIJKL[0][1][0][1] = D66;
    PMIJKL[0][1][1][0] = D66;
    PMIJKL[0][1][1][1] = -10.0;    


    PMIJKL[1][0][0][0] = -10.0;
    PMIJKL[1][0][0][1] = D66;
    PMIJKL[1][0][1][0] = D66;
    PMIJKL[1][0][1][1] = -10.0;        

    PMIJKL[1][1][0][0] = D12;
    PMIJKL[1][1][0][1] = -10.0;
    PMIJKL[1][1][1][0] = -10.0;
    PMIJKL[1][1][1][1] = D22;
    
  }
  
  fgets(fbuffer,MAX_BUFFER,fp);//printf("%s",fbuffer); // end
  fclose(fp);
  
} // closes read-mesh
  



// define allocate arrays

void mesh::allocmisc(){
   // allocate memory
  this->num_equations   = allocvec(nfields,this->num_equations);
  
  // dispersion stuff
  this->ibzcoord  = alloc3D(ndispersionseg,ndime,2,this->ibzcoord);
  this->w2karray  = allocmat(ndispersionseg,6,this->w2karray);

  this->ndispersionpts = allocvec(ndispersionseg,this->ndispersionpts);

  this->ndirich    = allocvec(nfields,this->ndirich);
  this->ndirdofn   = allocvec(nfields,this->ndirdofn);

  this->ntrac      = allocvec(nfields,this->ntrac);
  this->ntracdofn  = allocvec(nfields,this->ntracdofn);

  this->nperiodic  = allocvec(nfields,this->nperiodic);
  this->nperdofn   = allocvec(nfields,this->nperdofn);

  this->npforc     = allocvec(nfields,this->npforc);
  this->npforcdofn = allocvec(nfields,this->npforcdofn);

  this->elemsattached      = allocvec(this->num_nodes,this->elemsattached);
  zerovec(this->num_nodes,this->elemsattached);
}

void mesh::allocate_arrays(){

  /**
     Allocates memory for the arrays
  */
    
  coord    = allocmat(num_nodes,ndime,coord);
  coordmin = allocvec(ndime,coordmin); zerovec(ndime,coordmin);
  coordmax = allocvec(ndime,coordmax); zerovec(ndime,coordmax);
  meshdim  = allocvec(ndime,meshdim); zerovec(ndime,meshdim);
  conn     = allocmat(num_elem,mnode,conn);
  prop     = allocmat(num_prop_nodes,nprop,prop);
  propcons = alloc3D(num_prop_nodes,nprop,2,propcons);
  density  = allocmat(num_prop_nodes,1,density);
  lmref    = allocmat(num_elem,3,lmref);

  id       = alloc3D(num_nodes,mdofn,nfields,id);
  zero3D(num_nodes,mdofn,nfields,id);

  ien      = allocmat(num_elem,mnode,ien);
  zeromat(num_elem,mnode,ien);

  rhsvec = allocmat(num_equations_max,nfields,rhsvec);
  zeromat(num_equations_max,nfields,rhsvec);
  
  solrhsreal = allocmat(num_equations_max,nfields,solrhsreal);
  zeromat(num_equations_max,nfields,solrhsreal);

  solrhsimag = allocmat(num_equations_max,nfields,solrhsimag);
  zeromat(num_equations_max,nfields,solrhsimag);

  initstrn = allocvec(nfields,initstrn);
    
  /*allocate dirchlet bc array*/
  if(dirdofn !=0){
    dirich   = alloc3D(dirdofn,loadits,nfields,dirich);
  }
  
  /*allocate traction bc array*/
  if(tracdofn !=0){
    trac   = alloc3D(tracdofn,loadits,nfields,trac);
  }

  if(pdofn != 0 ){
    pforc = alloc3D(pdofn,loadits,nfields,pforc);
  }
  
  if(perdofn != 0 ){
    perbc  = alloc3D(perdofn,mdofn,nfields,perbc);
  }
  
  /*allocate isbc (is boundary condition ) array*/
  // isbc = allocmat(num_nodes,mdofn,isbc);
  isbc  = alloc3D(num_nodes, mdofn, nfields, isbc);
  zero3D(num_nodes,mdofn,nfields,isbc);

  ispointbc = alloc3D(num_nodes,mdofn,nfields,ispointbc);
  zero3D(num_nodes,mdofn,nfields,ispointbc);

  isperbc  = alloc3D(num_nodes, mdofn, nfields, isperbc);
  zero3D(num_nodes,mdofn,nfields,isperbc);

  solvec=alloc3D(num_nodes,mdofn,nfields,solvec);
  zero3D(num_nodes,mdofn,nfields,solvec);

  if ( solvetype != SOLVE_TYPE_BVP ) {
    // allocate eigenvalue and vectors
    eigenvecreal = alloc3D(num_nodes,mdofn,nfields,eigenvecreal);
    eigenvecimag = alloc3D(num_nodes,mdofn,nfields,eigenvecimag);
    
    eigenvalreal = allocmat(num_equations[this->which_field],nfields,eigenvalreal);
    eigenvalimag = allocmat(num_equations[this->which_field],nfields,eigenvalimag);
    eigenerror   = allocmat(num_equations[this->which_field],nfields,eigenerror);
  }

  // allocate dispersion stuff 
  eigdisperreal  = alloc3D(num_eigen,maxdispersionpts,ndispersionseg,eigdisperreal);
  eigdisperimag  = alloc3D(num_eigen,maxdispersionpts,ndispersionseg,eigdisperimag);
  eigdispererror = alloc3D(num_eigen,maxdispersionpts,ndispersionseg,eigdispererror);
  disperkvecs    = alloc3D(ndime,maxdispersionpts,ndispersionseg,disperkvecs);
  disperwvecs    = allocmat(maxdispersionpts,ndispersionseg,disperwvecs);

  zero3D(num_eigen,maxdispersionpts,ndispersionseg,eigdisperreal);
  zero3D(num_eigen,maxdispersionpts,ndispersionseg,eigdisperimag);
  zero3D(ndime,maxdispersionpts,ndispersionseg,disperkvecs);
  zero3D(ndime,maxdispersionpts,ndispersionseg,eigdispererror);

  lastupdate = alloc3D(num_nodes,mdofn,nfields,lastupdate);
  zero3D(num_nodes,mdofn,nfields,lastupdate);

  enhvec = alloc3D(num_elem,4,nfields,enhvec);
  zero3D(num_elem,4,nfields,enhvec);
  
  for(int ii = 0; ii < num_elem; ii++){
    for(int ienh =0; ienh < 4; ienh++){
      enhvec[ii][ienh][this->which_field] = 0.0;
    }
  }
  
  // allocate fillin
  fillin    = allocmat(num_equations_max,band,fillin);
  fillinno  = allocvec(num_equations_max,fillinno); 

  col_ptr      = allocvec(num_equations_max*band,col_ptr); // the column indices
  row_ptr        = allocvec(num_equations_max+1,row_ptr); //row pointers


  // allocate arrays needed for inversion
  functarray     = allocvec(nfields,functarray);
  regfunctarray  = allocvec(nfields,regfunctarray);

  zerovec(nfields,functarray);
  zerovec(nfields,regfunctarray);
  
  measdata       = alloc3D(num_nodes,mdofn,nfields,measdata);

  // fdgradgrad stores the gradient calculated by finite differences 
  // by function fdgrad
  fdgradgrad     = alloc3D(num_prop_nodes,nprop,nfields,fdgradgrad);
  zero3D(num_prop_nodes,nprop,nfields,fdgradgrad);
  
  // to store the dual(adjoint) field 
  dualfield      = alloc3D(num_nodes,mdofn,nfields,dualfield);
  zero3D(num_nodes,mdofn,nfields,dualfield);
  
  // to store the gradient by adjoint method
  // we store the gradient of each function w.r.t the properties
  // then we add it all up and feed it to TAO
  gradient = alloc3D(num_prop_nodes,nprop,nfields,gradient);
  
  // this combines all the fields in the gradient array into another array
  combgrad   = allocmat(this->num_prop_nodes,this->nprop,combgrad);
  combfdgrad = allocmat(this->num_prop_nodes,this->nprop,combfdgrad);

  zero3D(num_prop_nodes,nprop,nfields,gradient);

  // finally allocate the good guess array
  gooddispguess = alloc3D(num_nodes,mdofn,nfields,gooddispguess);

  // allocate the regparm array
  regparm = allocmat(this->nprop,this->mregparm,regparm);

  // allocate the tensor T in |T(U) - T(U^{m}|
  functens  = allocmat(this->ndime, this->ndime, functens);

  // allocate the weight vector
  fieldfuncweight = allocvec(this->nfields,fieldfuncweight);

  // allocate the vector for the norm of the measured data
  normmeasdata = allocvec(this->nfields,normmeasdata);
  zerovec(this->nfields,normmeasdata);

  // allocate data for armero hack
  //stiffoneelem = allocmat(8,8,stiffoneelem);
  //zeromat(8,8,stiffoneelem);

  EHIJKL = alloc4D(this->ndime,this->ndime,this->ndime,this->ndime,EHIJKL);
  EHAVIJKL = alloc4D(this->ndime,this->ndime,this->ndime,this->ndime,EHAVIJKL);
  PMIJKL = alloc4D(this->ndime,this->ndime,this->ndime,this->ndime,PMIJKL);
  GIJKL  = alloc4D(this->ndime,this->ndime,this->ndime,this->ndime,GIJKL);
  
  // alloc pz and directric tensors
  PZHKIJ   = alloc3D(this->ndime,this->ndime,this->ndime,PZHKIJ);
  PZHAVKIJ = alloc3D(this->ndime,this->ndime,this->ndime,PZHAVKIJ);

  KHIJ   = allocmat(this->ndime,this->ndime,KHIJ);
  KHAVIJ = allocmat(this->ndime,this->ndime,KHAVIJ);


  zero4D(this->ndime,this->ndime,this->ndime,this->ndime,EHIJKL);
  zero4D(this->ndime,this->ndime,this->ndime,this->ndime,EHAVIJKL);
  zero4D(this->ndime,this->ndime,this->ndime,this->ndime,PMIJKL);
  zero4D(this->ndime,this->ndime,this->ndime,this->ndime,GIJKL);

  zero3D(this->ndime,this->ndime,this->ndime,PZHKIJ);
  zero3D(this->ndime,this->ndime,this->ndime,PZHAVKIJ);

  zeromat(this->ndime,this->ndime,KHIJ);
  zeromat(this->ndime,this->ndime,KHAVIJ);


}
//==========================================================================
void mesh::combine_gradients(){
  // first zero the combined gradient vector
  zeromat(this->num_prop_nodes,this->nprop,combgrad);
  
  for (int inode = 0; inode < this->num_prop_nodes ; inode ++){
    for (int ifield = 0; ifield <this->nfields; ifield++){

      if ( this->opttype == 0 ) {
	for (int iprop = 0; iprop < this->nprop; iprop++){
	  this->combgrad[inode][iprop] += this->gradient[inode][iprop][ifield];
	}
      }
      else if ( this->opttype==1){
	this->combgrad[inode][0] += this->gradient[inode][0][ifield];
      }
      else if ( this->opttype==2){
	this->combgrad[inode][1] += this->gradient[inode][1][ifield];
      }

    }
  }

}
//==========================================================================
void mesh::combine_fdgradients(){
  // first zero the combined finite differences gradient vector
  zeromat(this->num_prop_nodes,this->nprop,combfdgrad);
  
  for(int inode = 0; inode < this->num_prop_nodes ; inode++ ){
    for (int ifield = 0; ifield < this->nfields; ifield++){

      if ( this->opttype == 0 ) {
	for (int iprop = 0; iprop < this->nprop; iprop++){
	  this->combfdgrad[inode][iprop] += this->fdgradgrad[inode][iprop][ifield];
	}
      }
      else if (this->opttype == 1 ){
	this->combfdgrad[inode][0] += this->fdgradgrad[inode][0][ifield];
      }
      else if (this->opttype == 2 ){
	this->combfdgrad[inode][1] += this->fdgradgrad[inode][1][ifield];
      }

    }
  }
}

//==========================================================================

void mesh::delete_stiff(){  
  // this is not called by the destructor - must be called explicitly by the user in forward.cpp
  delete stiff;  
  if (this->solvetype == SOLVE_TYPE_DISPERSIONK2W ){
    delete stiffimag;
  }
  
  if (this->solvetype == SOLVE_TYPE_DISPERSIONW2K){
    delete quadkr1; delete quadkr2; delete quadki;
  }

}      
//==========================================================================
void mesh::delete_mass(){  
  // this is not called by the destructor - must be called explicitly by the user in forward.cpp
  delete mass;    
  if (this->solvetype == SOLVE_TYPE_DISPERSIONK2W ){
    delete massimag;
  }
  
  if (this->solvetype == SOLVE_TYPE_DISPERSIONW2K){
    delete quadmr;
  }
  
}
//==========================================================================
void mesh::delete_damping(){
  // this is not called by the destructor - must be called explicitly by the user in forward.cpp
  if ( this->solvetype == SOLVE_TYPE_DISPERSIONW2K ){
    delete quadcreal;
    delete quadcimag;
  }
}

//==========================================================================
void mesh::deallocate_arrays(){
  // stiff, mass and damping csr mats are not destroyed here. Must be explicitly destroyed by the user in forward.cpp

  deallocmat(num_nodes,coord); //coord
  deallocvec(coordmin);
  deallocvec(coordmax);
  deallocvec(meshdim);
  deallocmat(num_elem,conn);   //conn
  deallocmat(num_prop_nodes,prop);  //prop
  deallocmat(num_prop_nodes,density);
  dealloc3D(num_prop_nodes,nprop,propcons);
  deallocmat(num_elem,lmref);  //lmref
  dealloc3D(num_nodes,mdofn,id);    //ld
  deallocmat(num_elem,ien);    //ien
  
  
  deallocmat(num_equations_max,rhsvec);         //rhsvec
  deallocmat(num_equations_max,solrhsreal);     //solrhs
  deallocmat(num_equations_max,solrhsimag);    

  deallocvec(initstrn);

  if(dirdofn !=0){
    // dellocate the vector containing the dirichlet boundary
    dealloc3D(dirdofn,loadits,dirich); 
  }

  if(tracdofn !=0){
    dealloc3D(tracdofn,loadits,trac);
  }

  if(pdofn !=0){
    dealloc3D(pdofn,loadits,pforc);
  }

  if (perdofn !=0){
    dealloc3D(perdofn,mdofn,perbc);
  }
  
  dealloc3D(num_nodes,mdofn,isbc);
  dealloc3D(num_nodes,mdofn,ispointbc);
  dealloc3D(num_nodes,mdofn,isperbc);

  dealloc3D(num_nodes,mdofn,solvec);

  if ( solvetype != SOLVE_TYPE_BVP ) {
    // allocate eigenvalue and vectors
    dealloc3D(num_nodes,mdofn,eigenvecreal);
    dealloc3D(num_nodes,mdofn,eigenvecimag);
    deallocmat(num_equations[this->which_field],eigenvalreal);
    deallocmat(num_equations[this->which_field],eigenvalimag);
    deallocmat(num_equations[this->which_field],eigenerror);
  }

  dealloc3D(num_nodes,mdofn,lastupdate);
  dealloc3D(num_elem,4,enhvec);

  deallocmat(num_equations_max,fillin);
  deallocvec(fillinno);
  deallocvec(col_ptr);
  deallocvec(row_ptr);
  
  // inversion stuff
  deallocvec(functarray);
  deallocvec(regfunctarray);
  dealloc3D(num_nodes,mdofn,measdata);
  dealloc3D(num_prop_nodes,nprop,fdgradgrad);
  dealloc3D(num_nodes,mdofn,dualfield);
  dealloc3D(num_prop_nodes,nprop,gradient);
  deallocmat(num_prop_nodes,combgrad);
  deallocmat(num_prop_nodes,combfdgrad);


  dealloc3D(num_nodes,mdofn,gooddispguess);
  deallocmat(this->nprop,regparm); 
  

  // dispersion stuff 
  dealloc3D(ndispersionseg,ndime,ibzcoord);
  deallocmat(ndispersionseg,w2karray);

  deallocvec(ndispersionpts);
  dealloc3D(num_eigen,maxdispersionpts,eigdisperreal);
  dealloc3D(num_eigen,maxdispersionpts,eigdisperimag);
  dealloc3D(num_eigen,maxdispersionpts,eigdispererror);
  dealloc3D(ndime,maxdispersionpts,disperkvecs);
  deallocmat(maxdispersionpts,disperwvecs);

  
  deallocmat(ndime,functens);
  deallocvec(fieldfuncweight);
  deallocvec(normmeasdata);
  deallocvec(num_equations);
  deallocvec(elemsattached);
  deallocvec(ndirich);  deallocvec(ntrac);  deallocvec(nperiodic);deallocvec(npforc);
  deallocvec(ndirdofn);deallocvec(ntracdofn);deallocvec(nperdofn);deallocvec(npforcdofn);

  // delete stiff; // delete the csr stiffness matrix
  //  this->delete_stiff();
  //deallocmat(8,stiffoneelem);

  dealloc4D(this->ndime,this->ndime,this->ndime,this->ndime,EHIJKL);
  dealloc4D(this->ndime,this->ndime,this->ndime,this->ndime,EHAVIJKL);

  dealloc4D(this->ndime,this->ndime,this->ndime,this->ndime,PMIJKL);
  dealloc4D(this->ndime,this->ndime,this->ndime,this->ndime,GIJKL);

  dealloc3D(this->ndime,this->ndime,PZHKIJ);
  deallocmat(this->ndime,KHIJ);
  dealloc3D(this->ndime,this->ndime,PZHAVKIJ);
  deallocmat(this->ndime,KHAVIJ);

  
}

//============================================================================
void mesh::make_ien(){

/**
makes the array ien in mesh
ien(element number, local node number) = global node number
Nachiket Gokhale gokhalen@bu.edu 15 Feb 2005,
Nachiket Gokhale gokhalen@bu.edu, converted into a class function, from a function
*/

 for(int ielem = 0; ielem<num_elem; ielem++){
    for(int inode = 0;inode <lmref[ielem][1];inode++){
      ien[ielem][inode] = conn[ielem][inode];
    }
  }

} 

//===================================================================

void mesh::make_eqn_num(){
/**
Nachiket Gokhale gokhalen@bu.edu 3rd March 2005, made class method
sets the (global node number,degree of freedom) = eqn num array 
*/

 for (int ifield = 0; ifield < nfields; ifield++){
   long int eqnno=1;

   // set equation numbers for unconstrained nodes
   for(long int inode = 0; inode < num_nodes;inode++){
     for(int idofn = 0;idofn < mdofn;idofn++){
       if( (!isdirichlet(inode,idofn,ifield)) && ( !isslave(inode,idofn,ifield)) ) {
	 id[inode][idofn][ifield]  = eqnno;
	 eqnno++;
       }
     }
   }

   // set equation numbers of slave nodes to those of their masters
   for ( long int inode = 0; inode < num_nodes; inode++){
     for ( int idofn = 0; idofn < mdofn; idofn++){
       if ( isslave(inode,idofn,ifield)){
	 long int masternode = getmasterperiodic(inode,idofn,ifield);
	 // force equation number
	 id[inode][idofn][ifield] = id[masternode-1][idofn][ifield];  
       } 
     }
   }

 }// ifield 

}

//===================================================================================== 
void mesh::make_fillin(){

  /*
   *  Makes the fillin for the stiffness matrix
   *  Nachiket Gokhale gokhalen@bu.edu Mar 3 2005
   *  This is done by looping over the elements 
   *  And has the unnecessary steps of localization of data into elements
   *  Nevertheless it provdides a clear procedure to make CSR arrays 
   *  
   *  This has the following arrays:
   *  
   *  1) fillinno[dofnnumber]  = number of degrees of freedom it is conn to
   *  2) fillin[ieqn][iband]   = the sorted list of dofns ieqn is connected to
   *
   *  The algorithm works this way:
   *
   *  1) Initialize: fillin[ieqn][iband] = ieqn+1;  //each equation sees itself
   *  2)             fillinno[ieqn]      = 1; //each equation sees itself
   *
   *  3) Loop over elements: (ielem)
   *  4)     Loop over equation numbers (ieqn) (check for non-dirichlet)
   *  5)          Loop over equation number (jeqn) (check non-dirchlet)
   *  6)             if entry exists then break; else add at the end
   */

  
  nnz = 0;
  // Initialize the fillin arrays
  for (long int ieqn = 0;ieqn<num_equations[this->which_field];ieqn++){

    for(int iband = 0;iband<band;iband++){
      fillin[ieqn][iband]=ieqn+1; // ieqn starts at zero
      
    }
    // because each equation sees (couples to) itself
    fillinno[ieqn] = 1; 
    nnz++;
  }
  
  
  // 
  for(long int ielem = 0;ielem<num_elem;ielem++){
    // this is deleted at the end of this loop
    elem_base *felem;
    // cout<<"caling getelemref"<<endl;
    felem = getelemref(this,ielem,felem);

    // dummy value for nonlin iteration
    // interested in localizing degree of freedom etc ( ) in here

    felem = localizer(this,ielem,felem,CURRENT_LOAD_ITERATION);

    for(int inode = 0;inode<felem->shpstrct->numnodes;inode++){
      for(int idofn = 0;idofn<felem->dofnatnode[inode];idofn++){

	int ieqn = felem->ideqn[inode][idofn];
	
	if(ieqn != 0){
	  for(int jnode = 0;jnode<felem->shpstrct->numnodes;jnode++){
	    for(int jdofn =0;jdofn<felem->dofnatnode[jnode];jdofn++){
	      int jeqn = felem->ideqn[jnode][jdofn];

	      if(jeqn !=0){
		int fillflag = 0; // assume the entry does not exist
		for(int ifill = 0; ifill<fillinno[ieqn-1]; ifill++){
		  if(fillin[ieqn-1][ifill] == jeqn){
		    fillflag =1; // the dofn has been counted
		  break;
		  }
		}
		
		if(fillflag == 0){
		  int oldnum = fillinno[ieqn-1]; // number of entries dofn sees
		  fillinno[ieqn-1]++;
		  fillin[ieqn-1][oldnum] = jeqn; // add new entry to the end.
		  nnz++; // increment nnz
		  // sorting step to get things in correct order
		  // the band width is small. therefor not using qsort
		  // stupidest sort is used // Numerical Recipies, 8.1 pg330.
		  // could be moved out of the ielem loop
		}
	      } // jeqn !=0
	    }//jdofn
	  }//jnode
	} // ieqn !=0
	
      }//idofn
    }//inode
    delete felem; // delete memory allocated. 
  }//ielem

  // Feb 10 - 2011 - Find max of fillinno
  {
    max_nnz_per_row = 0;
    for ( int ieqn = 0; ieqn < this->num_equations[this->which_field];ieqn++){
      if ( max_nnz_per_row < fillinno[ieqn] ){
	max_nnz_per_row = fillinno[ieqn];
      }
    }
    if ( max_nnz_per_row > this->band ){
      cout<<" Bandwidth in data.in is "<< band <<" - must be greater than exact bandwidth "<<max_nnz_per_row<<endl;
      cout<<" - Please - Fix - "<<endl;
      exit(1);
    }
  }

  // convert the the fillin to row_ptr

  row_ptr[0] = 0; // initialize

  for(int irow = 0;irow<num_equations[this->which_field];irow++){
    // points to the start of the next row
    // the last pointer is set correctly 
    row_ptr[irow+1] = row_ptr[irow] + fillinno[irow];
  }
  

  // Feb 9 - 2011 - This appears to be the sort 
  for(int ieqn = 0;ieqn<num_equations[this->which_field];ieqn++){
    int srtlength = fillinno[ieqn];
    for(int isort=1;isort<srtlength;isort++){
      int numout = fillin[ieqn][isort];
      int jsort  = isort -1;
      while((jsort >=0) && (fillin[ieqn][jsort] > numout)){
	fillin[ieqn][jsort+1] = fillin[ieqn][jsort];
	jsort--;
      }
      fillin[ieqn][jsort+1] = numout;
    }
  }

  // done the matrix is in csr format
  int icolpos = 0;
  for(int ieqn =0;ieqn<num_equations[this->which_field];ieqn++){
    for(int ifill =0;ifill<fillinno[ieqn];ifill++){
      // account for zero based indexing
      col_ptr[icolpos] = fillin[ieqn][ifill]-1; 
      
      icolpos++;
    }
  }
  
  assert(icolpos == nnz);
  
  // row and column pointers are constant - the value vector is allocated internally and stored
  stiff = (MATCSR *)new MATCSR(num_equations[this->which_field],band,row_ptr,col_ptr);
  mass  = (MATCSR *)new MATCSR(num_equations[this->which_field],band,row_ptr,col_ptr);

  if ( this->solvetype == SOLVE_TYPE_DISPERSIONK2W ){
    stiffimag = (MATCSR *)new MATCSR(num_equations[this->which_field],band,row_ptr,col_ptr);
    massimag = (MATCSR *)new MATCSR(num_equations[this->which_field],band,row_ptr,col_ptr);
    massimag->zeromat();
  }
  
  if ( this->solvetype == SOLVE_TYPE_DISPERSIONW2K ){
    quadmr    = (MATCSR *)new MATCSR(num_equations[this->which_field],band,row_ptr,col_ptr); 
    quadkr1   = (MATCSR *)new MATCSR(num_equations[this->which_field],band,row_ptr,col_ptr);
    quadkr2   = (MATCSR *)new MATCSR(num_equations[this->which_field],band,row_ptr,col_ptr);
    quadki    = (MATCSR *)new MATCSR(num_equations[this->which_field],band,row_ptr,col_ptr);
    quadcreal = (MATCSR *)new MATCSR(num_equations[this->which_field],band,row_ptr,col_ptr);
    quadcimag = (MATCSR *)new MATCSR(num_equations[this->which_field],band,row_ptr,col_ptr);
    
  }


}
//======================================================================================
// view sol
void mesh::view_and_store_sol(){

  // prints solution (displacement/primal) to sol.out

  // updateflag decides whether to increment solution or set solution
  // if we solve for the increment (as we do in non-linear elasticity)
  // the we need to increment the solution vector 
  // if we solve for the solution (a linear problem solved over and over) then we need to set sol
  
  int ifield = this->which_field;

  this->applybcstosol(this->solrhsreal);

  char outfilename[100];
  sprintf(outfilename,"sol%d-%d-%d.out",ifield,iloadits,inewtonits);
  ofstream outfile (outfilename);
  outfile.precision(DEFAULT_DOUBLE_PRECISION);
  if (outfile.is_open())
    {
     
	for(long int inode = 0;inode<this->num_nodes;inode++){
	  if ( this->printflag == OUTPUT_GNUPLOT){
	    for ( int idime = 0; idime < ndime; idime++){
	      outfile<<this->coord[inode][idime]<<" ";
	    }
	  }

	  for(long int idofn = 0;idofn<this->mdofn;idofn++){
  
	    outfile<<this->solvec[inode][idofn][ifield]<<"  ";
	  }
	  outfile<<endl;
	}
    }
     
  outfile.close();

}
//=======================================================================
void mesh::init_primal(){
  
  int ifield = this->which_field;

  // important: fdgradflag should be accounted for here.
  // while checking the finite difference gradient not accounting for this
  // causes the second primal guess to fail with negative jacobian determinant

  if(this->init_primal_goodguess_flag == false) {
    for(long int inode = 0;inode<this->num_nodes;inode++){
      for(long int idofn = 0;idofn<this->mdofn;idofn++){
	if(this->isbc[inode][idofn][ifield]>0){
	  int dirlocation = this->isbc[inode][idofn][ifield];
	  this->solvec[inode][idofn][ifield] = this->dirich[dirlocation-1][iloadits][ifield];
	}
	else{
	  this->solvec[inode][idofn][ifield] = 0;
	}
      }
    }
  } // if(this->primalflag)
  else {
      for(long int inode = 0;inode<this->num_nodes;inode++){
	for(long int idofn = 0;idofn<this->mdofn;idofn++){
	  this->solvec[inode][idofn][ifield] = this->gooddispguess[inode][idofn][ifield];
	}
      }
  }

}
//======================================================================
void mesh::zeroenhanced(){

  for(int ifield = 0; ifield < nfields; ifield++){
    for(int ielem = 0; ielem < this->num_elem; ielem++){
      for(int ienh = 0; ienh < 4; ienh++){
	this->enhvec[ielem][ienh][ifield] = 0.0;
      }
    }
  }
  
}
//======================================================================
void mesh::zerorhs(){
    for(int ieqn = 0;ieqn<this->num_equations[this->which_field];ieqn++){
      rhsvec[ieqn][this->which_field] = 0;
    }
}
//=============================================================================================
void mesh::printrhsnorm(){
  // print file: rhsnorm${which_field}.out
  
  
  this->rhsnorm = 0;
  char filename[100];
  sprintf(filename, "rhsnorm%d.out", this->which_field); 
  for(int ieqn = 0;ieqn<this->num_equations[this->which_field];ieqn++){
    this->rhsnorm += pow(rhsvec[ieqn][this->which_field],2);
  }
  
  this->rhsnorm = this->rhsnorm/this->num_equations[this->which_field];
  this->rhsnorm = sqrt(this->rhsnorm);

  // rhsnorm.out: stores norm of rhs (residual) at every iteration
  ofstream outfile (filename,ios::app);
  outfile.precision(DEFAULT_DOUBLE_PRECISION);
  outfile<<this->rhsnorm<<endl;
  outfile.close(); // close the file rhsnorm
  
}
// ================================================================================================

void mesh::writecoord(){

  ofstream  coordfile("coords.out");
  coordfile.precision(DEFAULT_DOUBLE_PRECISION);     

  for ( long int inode = 0; inode < num_nodes; inode++){
    for ( int idime = 0; idime < ndime ; idime++){
      coordfile<<coord[inode][idime]<<" ";
    }
    coordfile<<endl;
  }
  coordfile.close();
  
}
// ================================================================================================
void mesh::writeeigenvalues(){
  char filename[100];
  sprintf(filename,"eigenvalues-fld%d-seg%d-pt%d.out",this->which_field,this->which_disperseg,this->which_disperpt);
  ofstream eigvalfile(filename);
  eigvalfile.precision(DEFAULT_DOUBLE_PRECISION);
  eigvalfile<<"% Real\t\t\t Imag\t\t\t Relative |Ax-kBx|_2 "<<endl;
  eigvalfile<<"%     \t\t\t     \t\t\t    Error    |kBx|_2 "<<endl;  

  for ( int ieigen = 0; ieigen < this->num_eigen; ieigen++){
    eigvalfile<<this->eigenvalreal[ieigen][this->which_field]<<"\t";
    eigvalfile<<this->eigenvalimag[ieigen][this->which_field]<<"\t";
    eigvalfile<<this->eigenerror[ieigen][this->which_field]<<"\t"<<endl;
  }

  eigvalfile.close();
}
//=================================================================================================
void mesh::writeeigenvectorsvtk(int ieigen){

  char outfilename[100];
  sprintf(outfilename,"eigenvector-fld-%d-seg-%d-pt-%d-ieig-%d.vtk",this->which_field,this->which_disperseg,this->which_disperpt,ieigen);
  this->initvtk(&outfilename[0]);
  ofstream outfile(outfilename,ios_base::app);
  outfile.precision(DEFAULT_DOUBLE_PRECISION);
  int formatwidth = ceil(log10(this->num_eigen));
  if( outfile.is_open()){
    // write the eigenvector data - real data

    // vector data eigreal
    outfile<<endl<<"VECTORS VectorBlochReal"<<" double "<<endl;
    this->applybcstosol(this->solrhsreal);
    for(long int inode = 0; inode < this->num_nodes; inode++){
      for(int idofn = 0; idofn < this->mdofn; idofn++){
	outfile<<this->solvec[inode][idofn][this->which_field]<<" ";
      }
      if ( this->mdofn == 2 ){
	outfile<<" 0 ";
      }
      outfile<<endl;
    }

    // scalar data eigreal
    outfile<<endl<<"SCALARS ScalarBlochReal"<<" double "<<this->mdofn<<" "<<endl<<"LOOKUP_TABLE default"<<endl;
    for(long int inode = 0; inode < this->num_nodes; inode++){
      for(int idofn = 0; idofn < this->mdofn; idofn++){
	outfile<<this->solvec[inode][idofn][this->which_field]<<" ";
      }
      outfile<<endl;
    }
    

    // imag
    outfile<<endl<<"VECTORS VectorEigImag"<<" double "<<endl;
    this->applybcstosol(this->solrhsimag);
    for(long int inode = 0; inode < this->num_nodes; inode++){
      for(int idofn = 0; idofn < this->mdofn; idofn++){
	outfile<<this->solvec[inode][idofn][this->which_field]<<" ";
      }
      if ( this->mdofn == 2 ){
	outfile<<" 0 ";
      }
      outfile<<endl;
    }

     // scalar data eigimag
    outfile<<endl<<"SCALARS ScalarBlochImag"<<" double "<<this->mdofn<<" "<<endl<<"LOOKUP_TABLE default"<<endl;
    for(long int inode = 0; inode < this->num_nodes; inode++){
      for(int idofn = 0; idofn < this->mdofn; idofn++){
	outfile<<this->solvec[inode][idofn][this->which_field]<<" ";
      }
      outfile<<endl;
    }
    
    // u = u_0 exp(ikx) =  (ur + i ui)(c(th) + i s(th))     // th = k.x
    // (ur + i ui)*(ct + i st) = (ur ct - ui st )  + i ( ui ct + ur st )

    


  }
  outfile.close();
}
// ================================================================================================
void mesh::writeprimalrhs(){
  char filename[100];
  sprintf(filename,"rhs%d-%d-%d.out",this->which_field,this->iloadits,this->inewtonits);
  ofstream rhsfile(filename);
  rhsfile.precision(DEFAULT_DOUBLE_PRECISION);
  for(int inode = 0;inode<this->num_nodes;inode++){
    for(int idofn=0;idofn<this->mdofn;idofn++){
      long int ipos = this->id[inode][idofn][this->which_field] -1;
      if ( ipos >= 0 ){
	rhsfile<<this->rhsvec[ipos][this->which_field]<<" ";
      }
      else{
	rhsfile<<" bc ";
      }

    }
    rhsfile<<endl;
  }

  rhsfile.close();

}

//======================================================================
void mesh::readmeasdata(){

  for (int ifield = 0;ifield< nfields;ifield++){
    char filename[100];
    sprintf(filename,"meas%d.in",ifield);
    ifstream  measfile(filename);
    if (! measfile.is_open())  {
      cout << "Error opening file "<<filename; exit (1); 
    }
    // read num_nodes lines, mdofn at each line
    for(int inode = 0;inode<this->num_nodes;inode++){
      for(int idofn = 0;idofn<this->mdofn;idofn++){
	measfile>>measdata[inode][idofn][ifield];
	//cout<<measdata[inode][idofn][ifield];
	//measdata[inode][idofn][ifield] = 1.0;
      }
      //cout<<endl;
    }
    // close the file
    measfile.close();
    
  } //ifield
}
//=================================================================================================
void mesh::writefdgradgrad(){

  // open fdgrad${which_field}.out to write finite difference gradient
  char filename[100];
  sprintf(filename,"fdgrad%d.out",this->which_field);

  ofstream fdgradfile(filename);
  fdgradfile.precision(DEFAULT_DOUBLE_PRECISION);
  if(!fdgradfile.is_open()){
    cout<<" Error in opening file fdgrad.out for field "<<this->which_field;exit(1);
  }
  
  for(int inode = 0;inode<this->num_prop_nodes;inode++){
    if ( this->opttype == 0 ){
      for(int iprop = 0;iprop<this->nprop;iprop++){
	fdgradfile<<fdgradgrad[inode][iprop][this->which_field]<<" ";
      }
    }
    else if (this->opttype == 1 ){
      fdgradfile<<fdgradgrad[inode][0][this->which_field]<<" ";
    }
    else if (this->opttype == 2 ){
      fdgradfile<<fdgradgrad[inode][1][this->which_field]<<" ";
    }
    
    
    fdgradfile<<endl;
  }
  fdgradfile.close();
}
//================================================================================================
void mesh::writecombfdgrad(){
  // open fdgrad.out to write finite difference gradient.

  char filename[100];
  sprintf(filename,"fdgrad.out");

  ofstream fdgradfile(filename);
  
  fdgradfile.precision(DEFAULT_DOUBLE_PRECISION);
  if(!fdgradfile.is_open()){
    cout<<" Error in opening file fdgrad.out";exit(1);
  }
  
  for(int inode = 0;inode<this->num_prop_nodes;inode++){
    if ( this->opttype == 0 ){
      for(int iprop = 0;iprop<this->nprop;iprop++){
	fdgradfile<<combfdgrad[inode][iprop]<<" ";
      }
    }
    else if (this->opttype == 1){
      fdgradfile<<combfdgrad[inode][0]<<" ";
    }
    else if (this->opttype == 2){
      fdgradfile<<combfdgrad[inode][1]<<" ";
    }
    fdgradfile<<endl;
  }
  fdgradfile.close();
  
}

//================================================================================================
void mesh::writeonefieldfunctional(){
  // open functional${which_field}.out to write the functional
  char filename[100];
  sprintf(filename,"functional%d.out",this->which_field);
  
  ofstream funcfile(filename,ios::app);
  funcfile.precision(DEFAULT_DOUBLE_PRECISION);
  if(!funcfile.is_open()){
    cout<<" Error in opening file functional.out for field number " <<this->which_field ;exit(1);
  }
  
  funcfile<<this->functarray[this->which_field]<<"\t"<<this->regfunctarray[this->which_field]<<endl;
  funcfile.close();
}

//================================================================================================
void mesh::writefunctional(){
  // writes the combined functional of all fields : f_{1} + f_{2}
  // this is the functional that the optimization code sees.
  // open file
  ofstream funcfile("functional.out",ios::app);
  funcfile.precision(DEFAULT_DOUBLE_PRECISION);

  if(!funcfile.is_open()){
    cout<<" Error opening functional.out ";exit(1);
  }
  funcfile<<this->functional<<endl;
  funcfile.close();
}
//================================================================================================
void mesh::writedualrhs(){
// open file
  char filename[100];
  sprintf(filename,"dualrhs%d.out",this->which_field);
  ofstream dualrhsfile(filename);
  dualrhsfile.precision(DEFAULT_DOUBLE_PRECISION);
  if(!dualrhsfile.is_open()){
    cout<<" Error in opening file dual.out";exit(1);
  }

  int ieqn = 0;
  for(int inode = 0;inode<this->num_nodes;inode++){
    for(int idofn = 0;idofn<this->mdofn;idofn++){
      if(isbc[inode][idofn][this->which_field] <= 0){
	dualrhsfile<<rhsvec[ieqn][this->which_field]<<endl;
	ieqn++;
      }
    }
    
  }
  dualrhsfile.close();
}
//================================================================
void mesh::writedual(){
  // print filenames of the form dual01.out dual02.out etc
  char filename[100];
  sprintf(filename,"dual%d-%d-%d.out",this->which_field,this->iloadits,this->newtonits);

  ofstream dualfile(filename);
  dualfile.precision(DEFAULT_DOUBLE_PRECISION);
  if(!dualfile.is_open()){
    cout<<" Error in opening file dual.out";exit(1);
  }

  for(int inode = 0;inode<this->num_nodes;inode++){
    for(int idofn = 0;idofn<this->mdofn;idofn++){
      dualfile<<dualfield[inode][idofn][this->which_field]<<" ";
    }
    dualfile<<endl;
  }

  dualfile.close();
}
//================================================================================================
void mesh::writeproperties(){
// open file
  char outfilename[100];
  sprintf(outfilename,"prop%d.out",this->inverseiteration);

  ofstream propfile(outfilename);
  propfile.precision(DEFAULT_DOUBLE_PRECISION);
  if(!propfile.is_open()){
    cout<<" Error in opening file prop.out";exit(1);
  }

  
  for(int inode = 0;inode<this->num_prop_nodes;inode++){

    if ( this->printflag == OUTPUT_GNUPLOT){
      if (this->propdiscontinflag == 0 ){
	// for node based properties output node coords
	for ( int idime = 0; idime < ndime; idime++){
	  propfile<<this->coord[inode][idime]<<" ";
	}
      }
      else{
	// for element based properties output centroid
	for ( int idime = 0; idime < ndime; idime++){
	  double ave = 0.0;
	  for ( int ielemnode = 0; ielemnode < lmref[inode][1]; ielemnode++){
	    int inodetemp = conn[inode][ielemnode];
	    ave+=this->coord[inodetemp-1][idime];
	  }
	  ave=ave/lmref[inode][1];
	  propfile<<ave<<" ";
	} // idime
      }
    }
    
    for(int iprop = 0;iprop<this->nprop;iprop++){
      propfile<<prop[inode][iprop]<<" ";
    }
    propfile<<endl;
  }
  propfile.close();
}
//================================================================================================
void mesh::writeehijkl(){
  
  ofstream outfile("ehijkl.out");
  outfile.precision(DEFAULT_DOUBLE_PRECISION);
  if(!outfile.is_open()){
    cout<<" Error in opening file ehijkl.out";exit(1);
  }

  outfile <<"=== Writing the CIJKL tensor ===  "<<endl<<endl;

   for (int ii = 0; ii < this->ndime; ii++){
      for (int  jj = 0; jj < this->ndime; jj++){
	for (int kk = 0; kk < this->ndime; kk++){
	  for (int ll = 0; ll < this->ndime; ll++){
	    outfile<<" ii jj kk ll "<<ii+1 <<" "<<jj+1<<" "<<kk+1 <<" " << ll+1<< " "<<EHIJKL[ii][jj][kk][ll]<<endl;
	  } // ll
	} // kk 
      } // jj 
   } //ii


   outfile <<"=== Writing the D^{H} tensor in TWO dimensions "<<endl<<endl;
   
   double D11 = EHIJKL[0][0][0][0];
   double D12 = EHIJKL[0][0][1][1];
   double D13 = EHIJKL[0][0][0][1];
   double D22 = EHIJKL[1][1][1][1];
   double D23 = EHIJKL[1][1][0][1];
   double D66 = EHIJKL[0][1][0][1];

   double D11_PM = PMIJKL[0][0][0][0];
   double D12_PM = PMIJKL[0][0][1][1];
   double D13_PM = 0.0;
   double D22_PM = PMIJKL[1][1][1][1];
   double D23_PM = 0.0;
   double D66_PM = PMIJKL[0][1][0][1];

   double **DD,**EE;
   DD = allocmat(3,3,DD); EE = allocmat(3,3,EE);

   DD[0][0] = D11; DD[0][1] = D12; DD[0][2] = 0.0;
   DD[1][0] = D12; DD[1][1] = D22; DD[1][2] = 0.0; 
   DD[2][0] = 0.0; DD[2][1] = 0.0; DD[2][2] = D66;

   zeromat(3,3,EE);

   double determ = determinant(3,DD);

   if ( determ != 0 ){
     matinverse(3,DD,EE);
   }
   
   outfile << " D11 is "<<D11<<" D11_Target is "<<D11_PM<<endl;
   outfile << " D12 is "<<D12<<" D12_Target is "<<D12_PM<<endl;
   outfile << " D13 is "<<D13<<" D13_Target is "<<D13_PM<<endl;
   outfile << " D22 is "<<D22<<" D22_Target is "<<D22_PM<<endl;
   outfile << " D23 is "<<D23<<" D23_Target is "<<D23_PM<<endl;
   outfile << " D66 is "<<D66<<" D66_Target is "<<D66_PM<<endl;

   outfile<<endl;
   outfile <<"sqrt(D11*D22) == "<<sqrt(D11*D22)<<endl;
   outfile <<"sqrt(D11*D22)/D12 == "<<sqrt(D11*D22)/(D12)<<endl; 
   outfile <<"D11/D22 == "<<(D11/D22)<<endl;
   
   outfile<<endl;
   outfile<<"----------------------------"<<endl;
   outfile<<"Writing Difference between Target and Homogenized"<<endl;

   outfile<< "D11-Homogenized - D11_Target is "<<(D11 - D11_PM)<<endl;
   outfile<< "D12-Homogenized - D12_Target is "<<(D12 - D12_PM)<<endl;
   outfile<< "D13-Homogenized - D13_Target is "<<(D13 - D13_PM)<<endl;
   outfile<< "D22-Homogenized - D22_Target is "<<(D22 - D22_PM)<<endl;
   outfile<< "D23-Homogenized - D23_Target is "<<(D23 - D23_PM)<<endl;
   outfile<< "D66-Homogenized - D66_Target is "<<(D66 - D66_PM)<<endl;

   outfile<<endl;
   outfile<<"----------------------------"<<endl;
   outfile<<"Elastic Modulus Tensor is"<<endl;

   outfile << " E11 is "<<1.0/EE[0][0]<<endl;
   outfile << " E22 is "<<1.0/EE[1][1]<<endl;
   outfile << " E66 is "<<EE[2][2]<<endl;

   
   outfile<<endl;
   
   if ( this->physics == PHYSICS_PIEZO ) {

     for (int ii = 0; ii < this->ndime; ii++){
       for (int  jj = 0; jj < this->ndime; jj++){
	 for (int kk = 0; kk < this->ndime; kk++){
	   outfile<<"PZH kk ii jj "<<kk<<" "<<ii<<" " <<jj<<" "<<this->PZHKIJ[kk][ii][jj]<<endl;
	 }
       }
     }
     
     outfile<<endl;

     outfile<<"K11 = "<<this->KHIJ[0][0]<<endl;
     outfile<<"K22 = "<<this->KHIJ[1][1]<<endl;
     outfile<<"K12 = "<<this->KHIJ[0][1]<<endl;
     outfile<<"K21 = "<<this->KHIJ[1][0]<<endl;

   }

   outfile<<endl;

   outfile << " Total Mass is " <<this->totalmass<<endl;
   outfile << " Ave. Density is "<<(this->totalmass/this->unitcellvolume)<<endl;
   outfile << " Elem. Volume is "<<this->totalvolume<<endl;
   outfile << " ------------------------------------------------------------ " <<endl;
   
   // write out volume averaged properties
   double DAV11 = EHAVIJKL[0][0][0][0];
   double DAV12 = EHAVIJKL[0][0][1][1];
   double DAV13 = EHAVIJKL[0][0][0][1];
   double DAV22 = EHAVIJKL[1][1][1][1];
   double DAV23 = EHAVIJKL[1][1][0][1];
   double DAV66 = EHAVIJKL[0][1][0][1];

   outfile << " Writing volume averaged properties " << endl;
   outfile << " D11 (Vol Av) = "<<DAV11<<endl;
   outfile << " D12 (Vol Av) = "<<DAV12<<endl;
   outfile << " D13 (Vol Av) = "<<DAV13<<endl;
   outfile << " D22 (Vol Av) = "<<DAV22<<endl;
   outfile << " D23 (Vol Av) = "<<DAV23<<endl;
   outfile << " D66 (Vol Av) = "<<DAV66<<endl;

   
   if ( this->physics == PHYSICS_PIEZO ){
    double ee31 = this->PZHAVKIJ[1][0][0] ;
    double ee33 = this->PZHAVKIJ[1][1][1] ;
    double ee15 = this->PZHAVKIJ[0][0][1] ;
    double kk11 = this->KHAVIJ[0][0];
    double kk22 = this->KHAVIJ[1][1];

    outfile << "EE31 (Vol Av) = "<<ee31<<endl;
    outfile << "EE33 (Vol Av) = "<<ee33<<endl;
    outfile << "EE15 (Vol Av) = "<<ee15<<endl;
    outfile << "KK11 (Vol Av) = "<<kk11<<endl;
    outfile << "KK22 (Vol Av) = "<<kk22<<endl;
   }

   
   deallocmat(3,DD); deallocmat(3,EE);
   outfile.close();
}
//================================================================================================
void mesh::writeehijkl3d(){
  double **CC;
  CC = allocmat(6,6,CC);
  zeromat(6,6,CC);

  ofstream outfile("ehijkl.out");
  outfile.precision(DEFAULT_DOUBLE_PRECISION);
  if(!outfile.is_open()){
    cout<<" Error in opening file ehijkl.out";exit(1);
  }
  
  // \sigma = C \epsilon , where \sigma,\epsilon are 6-vectors,  
  // first three rows
  // see equation on Page 319 Sec 5.46 of Lai Rubin Krempl
  for ( int ii = 0; ii < 3; ii++){
    CC[ii][0] = EHIJKL[ii][ii][0][0]; 
    CC[ii][1] = EHIJKL[ii][ii][1][1];
    CC[ii][2] = EHIJKL[ii][ii][2][2];

    CC[ii][3] = EHIJKL[ii][ii][1][2];
    CC[ii][4] = EHIJKL[ii][ii][0][2];
    CC[ii][5] = EHIJKL[ii][ii][0][1];
  }
  
  // the next three rows of C
  for ( int ii = 3; ii < 6 ; ii++){
    int jj,kk,ll;    
    // strain vector ordering
    if (ii==3) { kk=2;ll=3;}
    if (ii==4) { kk=1;ll=3;}
    if (ii==5) { kk=1;ll=2;}

    CC[ii][0] = EHIJKL[0][0][kk-1][ll-1];
    CC[ii][1] = EHIJKL[1][1][kk-1][ll-1];
    CC[ii][2] = EHIJKL[2][2][kk-1][ll-1];
    CC[ii][3] = EHIJKL[1][2][kk-1][ll-1];
    CC[ii][4] = EHIJKL[0][2][kk-1][ll-1];
    CC[ii][5] = EHIJKL[0][1][kk-1][ll-1];
  }

  // write the C-tensor
  outfile<<"C == "<<endl;
  for ( int ii = 0; ii < 6; ii++){
    for ( int jj = 0; jj < 6; jj++){
      outfile<<CC[ii][jj]<<" ";
    }
    outfile<<endl;
  }
  outfile<<endl;

  deallocmat(6,CC);
  outfile.close();
  
}

//================================================================================================
void mesh::storedual(){
  // takes the dualfield stored in this->solrhs and puts it in this->dualfield

  int ieqn = 0;
  for(int inode = 0;inode<num_nodes;inode++){
    for(int idofn = 0;idofn<mdofn;idofn++){
      
      if(isbc[inode][idofn][this->which_field] <=0){ // check if not dirichlet
	dualfield[inode][idofn][this->which_field] = solrhsreal[ieqn][this->which_field];
	ieqn++;
      }
    }
  }

  assert(ieqn==this->num_equations[this->which_field]);
}


//=============================================================================
void mesh::writeonefieldgradient(){
  // writes the gradient to 'grad.out'
  
  char filename[100]; 
  char filename2[100];
  sprintf(filename,"grad%d.out",this->which_field);
  sprintf(filename2,"grad%d-%d.out",this->which_field,this->inverseiteration);
   
  ofstream gradfile(filename);
  ofstream gradfile2(filename2);

  gradfile.precision(DEFAULT_DOUBLE_PRECISION);
  gradfile2.precision(DEFAULT_DOUBLE_PRECISION);

  if(!gradfile.is_open()){
    cout<<"Error opening file grad.out for field number "<<this->which_field;exit(1);
  }

   if(!gradfile2.is_open()){
     cout<<"Error opening file "<<gradfile2<<" for field number "<<this->which_field;exit(1);
  }

  for(int inode = 0;inode<this->num_prop_nodes;inode++){
    if ( this->opttype == 0 ){
      for(int iprop= 0;iprop<this->nprop;iprop++){
	gradfile<<this->gradient[inode][iprop][this->which_field]<<" ";
	gradfile2<<this->gradient[inode][iprop][this->which_field]<<" ";
      }
    }
    else if ( this->opttype == 1 ){
      gradfile<<this->gradient[inode][0][this->which_field]<<" ";
      gradfile2<<this->gradient[inode][0][this->which_field]<<" ";
    }
    else if ( this->opttype == 2 ){
      gradfile<<this->gradient[inode][1][this->which_field]<<" ";
      gradfile2<<this->gradient[inode][1][this->which_field]<<" ";
    }
    gradfile<<endl;
    gradfile2<<endl;
  }

  gradfile.close();
  gradfile2.close();
}
//==========================================================================
void mesh::writegradient(){
  // writes the gradient to 'grad.out'
  char filename[100];
  sprintf(filename,"grad.out");
   
  ofstream gradfile(filename);
  gradfile.precision(DEFAULT_DOUBLE_PRECISION);
  
  if(!gradfile.is_open()){
    cout<<"Error opening file grad.out";exit(1);
  }

  for(int inode = 0;inode<this->num_prop_nodes;inode++){
    for(int iprop= 0;iprop<this->nprop;iprop++){
      gradfile<<this->combgrad[inode][iprop]<<" ";
    }
    gradfile<<endl;
  }

  gradfile.close();
}
//===========================================================================
void mesh::writemeasnorm(){
  // write the gradient to measnorm%d.out

  char filename[100];
  sprintf(filename,"measnorm%d.out",this->which_field);
  
  ofstream measnormfile(filename);
  measnormfile.precision(DEFAULT_DOUBLE_PRECISION);
  
  if(!measnormfile.is_open()){
    cout<<"Error opening file measnorm.out ";exit(1);
  }

  measnormfile<<sqrt(this->normmeasdata[this->which_field]);

}
//============================================================================
void mesh::writeenhancedparm(){
    char filename[100];
    sprintf(filename,"enhparm%d-%d-%d.out",this->which_field,this->iloadits,this->inewtonits);
    
    ofstream outfilename(filename);
    outfilename.precision(DEFAULT_DOUBLE_PRECISION);
    
    if(!outfilename.is_open()){
      cout<<"Error opening file enhparm.out"<<endl;exit(1);
    }
    
    // this will print enhanced parameters for all elements
    // if the mesh has traction elements then they don't have enhancements 
    // associated with them 
    // this will print out 0 for them
    for(int ielem = 0; ielem < this->num_elem;ielem++){
      for(int ienh = 0; ienh < 4; ienh++){
	outfilename<<this->enhvec[ielem][ienh][this->which_field]<<" ";
      }
      outfilename<<endl;
    }

    outfilename.close();
     
}
//============================================================================
void mesh::initvtk(char *outfilename){
  //char outfilename[100];
  //sprintf(outfilename,"mesh.vtk");
  ofstream outfile(outfilename);
  outfile.precision(DEFAULT_DOUBLE_PRECISION);
  if( outfile.is_open()){
    outfile<<"# vtk DataFile Version 3.0"<<endl<<"HyserIII VTK Outputfile"<<endl;
    outfile<<"ASCII"<<endl<<endl<<"DATASET UNSTRUCTURED_GRID"<<endl<<"POINTS "<<this->num_nodes<<" double"<<endl;

    // this will not work for meshes in which element types are mixed
    for(long int inode = 0; inode < this->num_nodes; inode++){

      for(int idime = 0; idime < this->ndime; idime++){
	outfile<<this->coord[inode][idime]<<" ";
      }

      if ( this->ndime == 3 ){
	outfile<<" "<<endl;
      }
      else{
	outfile<<" 0.0 "<<endl;
      }

    }
    
    int nodes_per_elem = 4;

    

    if ( this->ndime == 3){
      nodes_per_elem = 8; 
    }
    
    //outfile<<endl<<"CELLS "<<this->num_elem<<" "<<( (this->num_elem*nodes_per_elem) + this->num_nodes)<<endl;

    outfile<<endl<<"CELLS "<<this->num_elem<<" "<<( this->num_elem*(nodes_per_elem+1))<<endl;

    for(long int ielem = 0; ielem < this->num_elem; ielem++){
      outfile<<nodes_per_elem<<" ";
      for ( int inode = 0; inode < nodes_per_elem; inode++){
      outfile<<this->conn[ielem][inode]-1<<" ";
      }
      outfile<<endl;
    }

    // now the cell types
    outfile<<endl<<"CELL_TYPES "<<this->num_elem<<endl;
    for(long int ielem = 0; ielem < this->num_elem; ielem++){
      
      if( this->lmref[ielem][0] < 300) {
	outfile<<"9 "<<endl;
      }
      else if ( this->lmref[ielem][0] > 300 ){
	outfile<<"12  "<<endl;
      }
    }
  } // outfile is_open()
  // now the scalar data - first the \displacements
  outfile<<endl<<"POINT_DATA "<<this->num_nodes<<endl;

      
  outfile.close();

}
//============================================================================
void mesh::writevtk(){
  // write the final displacement/properties in VTK format - only quads supported
  char outfilename[100];
  sprintf(outfilename,"mesh.vtk");
  ofstream outfile(outfilename,ios_base::app);
  outfile.precision(DEFAULT_DOUBLE_PRECISION);
  if( outfile.is_open()){
    
        for ( int ifield = 0; ifield < this->nfields; ifield++){
      outfile<<endl<<"SCALARS Disp"<<ifield<<" double "<< this->mdofn<<" "<<endl<<"LOOKUP_TABLE default"<<endl;

      for(long int inode = 0; inode < this->num_nodes; inode++){
	for(int idofn = 0; idofn < this->mdofn; idofn++){
	  outfile<<this->solvec[inode][idofn][ifield]<<" ";
	}
	outfile<<endl;
      }

    }
    
    // now the material properties - this needs to accomodate "Cell Data" for element wise properties
    if ( this->propdiscontinflag == 1 ){outfile<<endl<<"CELL_DATA "<<this->num_prop_nodes<<endl;}
    outfile<<endl<<"SCALARS MaterialProperties"<<" double "<<this->nprop<<" "<<endl<<"LOOKUP_TABLE default"<<endl;
    for(long int inode =0; inode<this->num_prop_nodes; inode++){
      for(int iprop=0; iprop < this->nprop; iprop++){
	outfile<<this->prop[inode][iprop]<<" ";
      }
      outfile<<endl;
    }
  } // outfile.is_open()
      
  outfile.close();

}
//===========================================================================
void mesh::zerofunctional(){
  // zeroes "functional" and "functarray"
  functional = 0;
  zerovec(this->nfields,this->functarray);
  zerovec(this->nfields,this->regfunctarray);
}

//===========================================================================
void mesh::zerogradient(){
  // zeroes the gradient vector
  zero3D(this->num_prop_nodes,this->nprop,this->nfields,this->gradient);
  
}
//===========================================================================
void mesh::calculate_optvars(){
  // this will change and become more complicated 
  // when the ability to code which properties to 
  // optimize with respect to is picked.
  if ( this->opttype == 0 ) {
    this->numoptvars  = this->nprop*this->num_prop_nodes;
  }
  else{
    this->numoptvars = this->num_prop_nodes;
  }
}
//===========================================================================
void mesh::updategoodguess(){

  for (int inode = 0; inode < this->num_nodes; inode++){
    for(int idofn =0; idofn < this->mdofn; idofn++){
      this->gooddispguess[inode][idofn][which_field] = this->solvec[inode][idofn][which_field];
    }
  }

}
//============================================================================
void mesh::resetprops(){
  if ( this->opttype !=0 ){
    for ( int inode = 0; inode < num_prop_nodes; inode++){
      
      if ( this->opttype == 1 ){
	this->prop[inode][1] = this->lamtomu*prop[inode][0];
      }
      else if (this->opttype == 2){
	this->prop[inode][0] = this->mutolam*prop[inode][1];
      }

    }
  }
}
//============================================================================
bool  mesh::isenhanced(int ielem){
  int itype = this->lmref[ielem][0];
  
  if ( (itype == NONLINELAS_PILTNERTAYLOR_IJNME99_2D ) || (itype == NONLINELASARMEROGLASER2D) || (itype == NONLINELASWRIGGERSREESE2D) ){
    return true;
  }
  else{
    return false;
  }

  // never get here
  return false;
}
//==============================================================================
void mesh::computeaveprops(){
  // bounding box volume must be correct for this to work

  for (int ielem = 0 ; ielem < this->num_elem; ielem++){
      elem_base *felement;
      felement = getelemref(this,ielem,felement);
      felement = localizer(this,ielem,felement,LAST_LOAD_ITERATION);    
      // get integration points
      felement->getinteg();
      
      // initialize shape functions 
      felement->getshape();
      felement->computeaveprops();
      
      if ((this->physics == PHYSICS_ELAS )|| (this->physics==PHYSICS_PIEZO)){
	for (int ii = 0; ii < this->ndime; ii++){
	  for (int  jj = 0; jj < this->ndime; jj++){
	    for (int kk = 0; kk < this->ndime; kk++){
	      for (int ll = 0; ll < this->ndime; ll++){ 
		this->EHAVIJKL[ii][jj][kk][ll] += (felement->BASE_EHAVIJKL[ii][jj][kk][ll]/this->unitcellvolume);
	      }//ll
	    }//kk
	  }//jj
	}//ii
      }// elas or piezo

      if ( this->physics == PHYSICS_PIEZO){

	for (int kk = 0; kk < this->ndime; kk++){
	  for (int ii = 0; ii < this->ndime; ii++){
	    for (int  jj = 0; jj < this->ndime; jj++){	
	      this->PZHAVKIJ[kk][ii][jj] += (felement->BASE_PZHAVKIJ[kk][ii][jj]/this->unitcellvolume);
	    }
	  }
	}

	for (int ii = 0; ii < this->ndime; ii++){
	  for (int  jj = 0; jj < this->ndime; jj++){	
	    this->KHAVIJ[ii][jj] += (felement->BASE_KHAVIJ[ii][jj]/this->unitcellvolume);
	  }
	}
	
      }// only piezo
    delete felement;

  } //ielem
  
  
}
//================================================================================
void mesh::computevolumemass(){
  this->totalmass = 0.0;
  this->totalvolume = 0.0;
  
  for (int ielem = 0 ; ielem < this->num_elem; ielem++){
    elem_base *felement;
    felement = getelemref(this,ielem,felement);
    felement = localizer(this,ielem,felement,LAST_LOAD_ITERATION);    

    // get integration points
    felement->getinteg();

    // initialize shape functions 
    felement->getshape();

    felement->computevolumemass();
    this->totalmass   += felement->elemmass;
    this->totalvolume += felement->elemvolume;

    delete felement;
  }
}
//===============================================================================================================================
void mesh::getbcs(long int ***isbc_arg, double  ***dirich_arg, long int *nbcs_arg, int bctype_arg,FILE *fp){
  int itemp, isdirbc;
  char fbuffer[MAX_BUFFER];
  double dvalue;
  long int litemp;
  fgets(fbuffer,MAX_BUFFER,fp);printf("%s",fbuffer);

  int bcfactor = 1;
  if ( bctype_arg == BC_TYPE_TRACTION ) {bcfactor = -1;}
 

  for ( int ifield = 0; ifield < nfields; ifield++){
    long int dirctr=1;

    if  ( nbcs_arg[ifield] != 0 ){
      for (int inode = 0; inode < nbcs_arg[ifield]; inode++){
	int dnode;
	fscanf(fp,"%d\n",&dnode); 
	if ( dnode < 0 ) {
	  fgets(fbuffer,MAX_BUFFER,fp);printf("%s",fbuffer);
	  fscanf(fp,"%d\n",&dnode); 
	}

	if (( bctype_arg == BC_TYPE_DIRICHLET ) || ( bctype_arg == BC_TYPE_TRACTION )){
	    for (int idofn =0; idofn < mdofn; idofn++){
	      fscanf(fp,"%d ", &isdirbc); 
	      printf("\n dnode= %d, idofn == %d isdirbc== %d",dnode,idofn, isdirbc); 

	      for(int ireadload = 0; ireadload < this->loadits; ireadload++){
		fscanf(fp,"%lf",&dvalue);
		printf(" dvalue == %lf\n",dvalue);
		if( isdirbc == 1 ){
		  if ( ireadload == 0 ){
		    isbc_arg[dnode-1][idofn][ifield] = bcfactor*dirctr;
		    dirctr++;
		  }
		  int dirlocation = isbc_arg[dnode-1][idofn][ifield]-1;
		  dirich_arg[dirlocation][ireadload][ifield] = dvalue;
		}
		fscanf(fp,"\n");
	      } // ireadload
	    }//idofn
	}

	if (bctype_arg == BC_TYPE_PERIODIC){
	  long int masternode, slavenode;
	  masternode = dnode;
	  fscanf(fp,"%ld ",&slavenode);
	  for (int idofn =0; idofn < mdofn; idofn++){
	  fscanf(fp,"%d ",&isdirbc);

	  // check if multiple constraints
	  if ( this->isbc[masternode-1][idofn][ifield] > 0 ) { 
	    cout<<"Periodic Node "<<dnode<<" dofn "<<idofn<<" has Dirichlet constraints"<<endl; 
	  }
	  else if ( this->isbc[masternode-1][idofn][ifield] < 0 ){
	    cout<<"Periodic Node "<<dnode<<" dofn "<<idofn<<" has Traction constraints"<<endl; 
	  }
	  else{
	    // There the slavenodes are tied back to the masters
	    // The masters are treated as normal free finite element nodes, they have no constraints on them
	    isbc_arg[slavenode-1][idofn][ifield]   = dirctr;
	    perbc[dirctr-1][idofn][ifield] = masternode;
	    dirctr++;
	    //cout<<"Periodic Node = "<<dnode<<"Slave Node = "<<slavenode<<endl;
	  }
	  } // idofn
	  fscanf(fp,"\n");
	}
	
      } // inode
    }//if nbcs

      // throw away zero
    fscanf(fp,"%ld\n",&litemp); printf("litemp %ld\n",litemp);
    printf("counted %ld dirchlet dofns \n",dirctr-1);fflush(stdout);
  }//ifield
} // ends function getbcs
//==========================================================================================
bool mesh::isdirichlet(long int inode, int idofn, int which_field){
  if ( this->isbc[inode][idofn][which_field] > 0 ) { 
    return true ; 
  }
  return false;
}
//===========================================================================================
bool mesh::istraction(long int inode, int idofn, int which_field){
  if ( this->isbc[inode][idofn][which_field] < 0 ){
    return true;
  }
  return false;
}
//==========================================================================================
bool mesh::isslave(long int inode,int idofn, int which_field){
  if ( this->isperbc[inode][idofn][which_field]>0){
    // slave node in periodic detected
    return true;
  }
  return false;
}
//===================================================================================================
double mesh::getdirichlet(long int inode,int idofn, int which_loadit,int which_field){
  int dirlocation = this->isbc[inode][idofn][which_field];
  return this->dirich[dirlocation-1][which_loadit][which_field];
}
//==================================================================================================
double mesh::gettraction(long int inode, int idofn, int which_loadit, int which_field){
  int tracloc =  -this->isbc[inode][idofn][which_field];
  return this->trac[tracloc-1][which_loadit][which_field];
}
//==================================================================================================
long int mesh::getmasterperiodic(long int inode, int idofn, int which_field){
  int perloc = this->isperbc[inode][idofn][which_field];
  return this->perbc[perloc-1][idofn][this->which_field];
}
//===============================================================================================
void mesh::writedispersion(){
  
  
  if ( this->solvetype == SOLVE_TYPE_DISPERSIONK2W ){
    char outfilename[100];
    sprintf(outfilename,"dispersionk2w.out");
    ofstream outfile(outfilename);
    outfile.precision(DEFAULT_DOUBLE_PRECISION);
    outfile<<"% Number EigReal  EigImag   Kx Ky  Error "<<endl;

    if (outfile.is_open()){
      long int ictr = 1;
      for(int iseg = 0; iseg < ndispersionseg ; iseg++){
	outfile<<"%% Dispersion segment "<<iseg<<endl;
	for ( int ik = 0; ik < ndispersionpts[iseg]; ik++){
	  for ( int ieigen = 0; ieigen < this->num_eigen; ieigen++){
	    
	    outfile<<ictr<<" "; 
	    
	    outfile<<eigdisperreal[ieigen][ik][iseg]<<" ";
	    outfile<<eigdisperimag[ieigen][ik][iseg]<<" ";
	    
	    for (int idime = 0; idime < ndime ; idime++){
	      outfile<< disperkvecs[idime][ik][iseg]<<" "; 
	    }
	    outfile<<eigdispererror[ieigen][ik][iseg]<<" ";
	    outfile<<endl;
	    
	  } // ieigen
	  ictr++;
	  
	} //ik
	
      }// iseg
    
    }//isopen
    outfile.close();

  }//k2w

  if ( this->solvetype == SOLVE_TYPE_DISPERSIONW2K){
    char outfilename[100];
    sprintf(outfilename,"dispersionw2k.out");
    ofstream outfile(outfilename);
    outfile.precision(DEFAULT_DOUBLE_PRECISION);
    outfile<<"% Number EigReal  EigImag   ww nx ny  Error "<<endl;

    if ( outfile.is_open()){
      long int ictr = 1;
      for ( int iseg = 0; iseg < ndispersionseg ; iseg++){
	outfile<<"%% Dispersion segment " <<iseg<<endl;
	for ( int ik = 0; ik < ndispersionpts[iseg]; ik++){
	  for ( int ieigen = 0; ieigen < this->num_eigen; ieigen++){
	    
	    outfile<<ictr<<" "; 
	    
	    outfile<<eigdisperreal[ieigen][ik][iseg]<<" ";
	    outfile<<eigdisperimag[ieigen][ik][iseg]<<" ";
	    

	    outfile<< disperwvecs[ik][iseg]<<" "; 
	    outfile<< w2karray[iseg][4]<<" "; 
	    outfile<< w2karray[iseg][5]<<" "; 

	    outfile<<eigdispererror[ieigen][ik][iseg]<<" ";
	    //outfile<<"-3.14150000"<<" ";
	    outfile<<endl;
	    
	  } // ieigen
	  ictr++;
	  
	} //ik
      }//iseg
    }//isopen
  } // w2k


}
//===============================================================================================
void mesh::computedisperkvecs(){
  // compute the kvecs for the dispersion curve

  for ( int iseg = 0;  iseg < ndispersionseg ; iseg++ ){
    double dk[3];
    
    for ( int idime = 0; idime < ndime ; idime++){
      dk[idime] = (ibzcoord[iseg][1][idime] - ibzcoord[iseg][0][idime]);
      dk[idime]  = dk[idime]/(ndispersionpts[iseg]-1);
    }
    
    for ( int ipt = 0; ipt < ndispersionpts[iseg]; ipt++){
      for ( int idime = 0; idime < ndime ; idime++){
	disperkvecs[idime][ipt][iseg] = ibzcoord[iseg][0][idime] + (ipt)*dk[idime];
      }
    }//ipt

  } //iseg

}
//====================================================================================
void mesh::computedisperwvecs(){
  // compute the wvecs for the dispersion curve

  for ( int iseg = 0;  iseg < ndispersionseg ; iseg++ ){
    double df;
    df = ( w2karray[iseg][1] - w2karray[iseg][0] ) / (ndispersionpts[iseg]-1);
    for ( int ipt = 0; ipt < ndispersionpts[iseg]; ipt++){
      disperwvecs[ipt][iseg] = w2karray[iseg][0] + (ipt)*df;
      disperwvecs[ipt][iseg] = disperwvecs[ipt][iseg]*2*M_PI;
      //cout<<"disperwvec= "<<disperwvecs[ipt][iseg]<<" ipt"<<" iseg ="<<endl;
    }
  } //iseg

}
//==============================================================================================
void mesh::seteigendisper(int ieigen,double kreal,double kimag, double error){
  this->eigdisperreal[ieigen][this->which_disperpt][this->which_disperseg] = kreal;
  this->eigdisperimag[ieigen][this->which_disperpt][this->which_disperseg] = kimag;
  this->eigdispererror[ieigen][this->which_disperpt][this->which_disperseg] = error;

}
//==============================================================================================
void mesh::computenumeqns(){
  this->num_equations_max  = 0; 
  for ( int ifield = 0 ; ifield < this->nfields; ifield++){
    
    num_equations[ifield] =  mdofn*num_nodes - ndirdofn[ifield] - nperdofn[ifield];

    if ( this->solvetype == SOLVE_TYPE_DISPERSIONK2W ){
      num_equations[ifield] = num_equations[ifield];
    }
    
    if ( num_equations_max < num_equations[ifield] ){
      num_equations_max     =  num_equations[ifield];
    } // calculate number of equations
  }

}
//==============================================================================================
bool mesh::issolvedispersion(){
  if (( this->solvetype == SOLVE_TYPE_DISPERSIONK2W ) || ( this->solvetype == SOLVE_TYPE_DISPERSIONW2K)){
    return true;
  }
  else{
    return false;
  }
}
//===============================================================================================
bool mesh::issolveeigen(){
  bool output = false;

  if ( this->solvetype == SOLVE_TYPE_DISPERSIONW2K ) {     output = true;   }
  if ( this->solvetype == SOLVE_TYPE_DISPERSIONK2W ) {     output = true;   }
  if ( this->solvetype == SOLVE_TYPE_EIGENK )        {     output = true;   }
  if ( this->solvetype == SOLVE_TYPE_EIGENM )        {     output = true;   }
  if ( this->solvetype == SOLVE_TYPE_EIGEN  )        {     output = true;   }

  return output;
}

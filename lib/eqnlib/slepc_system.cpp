#ifndef HAVE_SLEPC_SYSTEM_HPP
#include "slepc_system.hpp"
#endif

// constructors and destructor
slepc_system::slepc_system(){}
slepc_system::slepc_system(int num_equations_arg, int which_field, int inewton, int iload,int iprintflag){
  
  this->num_equations = num_equations_arg;
  this->which_field   = which_field;
  this->inewton       = inewton;
  this->iload         = iload;
  this->printflag     = iprintflag;
}

slepc_system::slepc_system(equation_initdata eqninitdata){
  this->num_equations      = eqninitdata.num_equations;
  this->num_equations_real = eqninitdata.num_equations_real;
  this->num_equations_imag = eqninitdata.num_equations_imag;
  this->max_nnz_per_row    = eqninitdata.max_nnz_per_row;

  this->which_field   = eqninitdata.which_field;
  this->inewton       = eqninitdata.inewton;
  this->iload         = eqninitdata.iload;
  this->printflag     = eqninitdata.printflag;
}

slepc_system::~slepc_system(){
  ierr = VecDestroy(&eigvalreal);  CHKERRXX(ierr);
  ierr = VecDestroy(&eigvalimag);  CHKERRXX(ierr);
  ierr = VecDestroy(&eigenerror);  CHKERRXX(ierr);
  ierr = EPSDestroy(&eps);         CHKERRXX(ierr);

  if ( this->allocstiff == 1 ){ ierr = MatDestroy(&this->stiff); CHKERRXX(ierr);}
  if ( this->allocmass  == 1 ){ ierr = MatDestroy(&this->mass); CHKERRXX(ierr);}
  if ( this->allocstiffimag == 1 ){ 
    ierr = MatDestroy(&this->stiff_real); CHKERRXX(ierr);
    ierr = MatDestroy(&this->stiff_imag); CHKERRXX(ierr);
    ierr = MatDestroy(&this->stiff_complex); CHKERRXX(ierr);
  }
  if (this->allocmassimag == 1 ) {
   ierr = MatDestroy(&this->mass_real); CHKERRXX(ierr);
   ierr = MatDestroy(&this->mass_imag); CHKERRXX(ierr);    
   ierr = MatDestroy(&this->mass_complex); CHKERRXX(ierr);
  }
}
//================================================================================================
void slepc_system::createcsrk(int *row_ptr,int *col_ptr,double *value_arr){
  // create the K matrix 
  char filename[100];
  sprintf(filename,"matk%d%d%d.m",this->which_field,this->inewton,this->iload);

  MatCreateSeqAIJWithArrays(PETSC_COMM_SELF,num_equations,num_equations,row_ptr,col_ptr,value_arr,&stiff);
  this->allocstiff = 1;
  
  if(this->printflag == 1){
    PetscViewerASCIIOpen(PETSC_COMM_SELF,filename,&view_mat);
    PetscViewerSetFormat(view_mat,PETSC_VIEWER_ASCII_MATLAB);
    MatView(stiff,view_mat);
    PetscViewerDestroy(&view_mat);
  }
}
//=======================================================================================
void slepc_system::createcsrkreal(int *row_ptr,int *col_ptr,double *value_arr){
  MatCreateSeqAIJWithArrays(PETSC_COMM_SELF,num_equations_real,num_equations_real,row_ptr,col_ptr,value_arr,&stiff_real);
  this->allocstiffimag = 1;
}
//=======================================================================================
void slepc_system::createcsrkimag(int *row_ptr,int *col_ptr,double *value_arr){
 MatCreateSeqAIJWithArrays(PETSC_COMM_SELF,num_equations_imag,num_equations_imag,row_ptr,col_ptr,value_arr,&stiff_imag);
 this->allocstiffimag = 1;
}

//=================================================================================================
void slepc_system::createcsrm(int *row_ptr,int *col_ptr,double *value_arr){
  char filename[100];
  sprintf(filename,"matm%d%d%d.m",this->which_field,this->inewton,this->iload);
  MatCreateSeqAIJWithArrays(PETSC_COMM_SELF,num_equations,num_equations,row_ptr,col_ptr,value_arr,&mass);
  this->allocmass = 1;

  if(this->printflag == 1){
    PetscViewerASCIIOpen(PETSC_COMM_SELF,filename,&view_mat);
    PetscViewerSetFormat(view_mat,PETSC_VIEWER_ASCII_MATLAB);
    MatView(mass,view_mat);
    PetscViewerDestroy(&view_mat);
  }
    
}
//===============================================================================
void slepc_system::createcsrmreal(int *row_ptr,int *col_ptr,double *value_arr){
  this->allocmassimag = 1;
MatCreateSeqAIJWithArrays(PETSC_COMM_SELF,num_equations_imag,num_equations_imag,row_ptr,col_ptr,value_arr,&mass_real);
}
//===============================================================================
void slepc_system::createcsrmimag(int *row_ptr,int *col_ptr,double *value_arr){
  this->allocmassimag = 1;
MatCreateSeqAIJWithArrays(PETSC_COMM_SELF,num_equations_imag,num_equations_imag,row_ptr,col_ptr,value_arr,&mass_imag);
}
//==================================================
void slepc_system::createrhs(double **value_arr){
  //  not needed for eigenvalue computations
  
}
//====================================================================
int slepc_system::solve(){
  cout<<"Beginning of  Slepc_System::Solve"<<endl;
  PetscReal error, tol, re, im, kerror;
  PetscScalar kr, ki;
  PetscInt    	 nev, maxit, i, its, lits;
  char        	 filename[256];
  PetscViewer 	 viewer;
  PetscBool  	 flg;
  Vec            TempVec;
  PetscScalar    TempScalar;

  nev = this->num_eigen;


  MatGetVecs(stiff,&eigvalreal,PETSC_NULL);
  MatGetVecs(stiff,&eigvalimag,PETSC_NULL);
  MatGetVecs(stiff,&eigenerror,PETSC_NULL);
   

  ierr = EPSCreate(PETSC_COMM_WORLD,&eps);CHKERRQ(ierr);
  ierr = EPSSetWhichEigenpairs(eps,EPS_SMALLEST_REAL); CHKERRQ(ierr);
  ierr = EPSSetTolerances(eps,LINEAR_SYS_REL_TOL,PETSC_DECIDE);CHKERRQ(ierr);
  // override number of eigenvalues to compute 
  ierr = EPSSetDimensions(eps,nev,PETSC_DECIDE,PETSC_DECIDE);CHKERRQ(ierr);

  if (( this->solvetype == SOLVE_TYPE_EIGEN )||(this->solvetype==SOLVE_TYPE_DISPERSIONK2W)) {
    ierr =  EPSSetProblemType(eps,EPS_GHEP);
    ierr = EPSSetOperators(eps,stiff,mass);CHKERRQ(ierr);
    
  }
  else{
    ierr = EPSSetProblemType(eps,EPS_HEP);
    if ( this->solvetype == SOLVE_TYPE_EIGENK) {
      ierr = EPSSetOperators(eps,stiff,PETSC_NULL);CHKERRQ(ierr);
    }
    else{
      ierr = EPSSetOperators(eps,mass,PETSC_NULL);CHKERRQ(ierr);
    }
  }

  ierr = EPSSetFromOptions(eps);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Starting eigensolve  ...\n");CHKERRQ(ierr);
  cout<<"About to Call EPSSolve  Slepc_System::Solve"<<endl;
  ierr = EPSSolve(eps);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"...done with eigensolve.\n");CHKERRQ(ierr);

  //    Optional: Get some information from the solver and display it
  ierr = EPSGetIterationNumber(eps, &its);CHKERRQ(ierr);
  ierr = EPSGetOperationCounters(eps,PETSC_NULL,PETSC_NULL,&lits);CHKERRQ(ierr);
  ierr = EPSGetType(eps,&epstype);CHKERRQ(ierr);
  ierr = EPSGetDimensions(eps,&nev,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
  ierr = EPSGetConverged(eps,&nconv);CHKERRQ(ierr);
  ierr = EPSView(eps,PETSC_VIEWER_STDOUT_SELF);CHKERRQ(ierr);
  
  ierr = EPSGetTolerances(eps,&tol,&maxit);CHKERRQ(ierr);


//ierr = PetscPrintf(PETSC_COMM_WORLD," Number of iterations of the method: %d\n",its);CHKERRQ(ierr);
  //ierr = PetscPrintf(PETSC_COMM_WORLD," Number of linear iterations of the method: %d\n",lits);CHKERRQ(ierr);
  //ierr = PetscPrintf(PETSC_COMM_WORLD," Solution method: %s\n\n",epstype);CHKERRQ(ierr);
//ierr = PetscPrintf(PETSC_COMM_WORLD," Number of requested eigenvalues: %d, number of converged = %d\n",nev,nconv);CHKERRQ(ierr);
  //ierr = PetscPrintf(PETSC_COMM_WORLD," Stopping condition: tol=%.4g, maxit=%d\n",tol,maxit);CHKERRQ(ierr);
  

  for ( int ieigen = 0; ieigen < min(nconv,this->num_eigen); ieigen++){
    ierr = EPSGetEigenpair(eps,ieigen,&kr,&ki,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
    ierr = EPSComputeRelativeError(eps,ieigen,&kerror);
    ierr = VecSetValue(eigvalreal,ieigen,kr,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(eigvalimag,ieigen,ki,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(eigenerror,ieigen,kerror,INSERT_VALUES);CHKERRQ(ierr);
    
  }
  return 0;
}

// =======================================================================================
int slepc_system::solvecomplex(){
  // create and solve the matrix system 
  PetscInt nrows = num_equations_real + num_equations_imag;
  PetscInt ncols = num_equations_real + num_equations_imag;
  PetscInt nnzcomplex = 2*this->max_nnz_per_row;

  MatCreateSeqAIJ(PETSC_COMM_SELF,nrows,ncols,nnzcomplex,PETSC_NULL,&stiff_complex);
  MatCreateSeqAIJ(PETSC_COMM_SELF,nrows,ncols,nnzcomplex,PETSC_NULL,&mass_complex);
    
  assemblecomplexmat(&stiff_real,&stiff_imag,&stiff_complex);
  assemblecomplexmat(&mass_real,&mass_imag,&mass_complex);

  ierr=MatAssemblyBegin(stiff_complex,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr=MatAssemblyEnd(stiff_complex,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

  ierr=MatAssemblyBegin(mass_complex,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);  
  ierr=MatAssemblyEnd(mass_complex,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

  cout<<"Beginning of  Complex Slepc_System::Solve"<<endl;
  PetscReal error, tol, re, im, kerror;
  PetscScalar kr, ki;
  PetscInt    	 nev, maxit, i, its, lits;
  char        	 filename[256];
  PetscViewer 	 viewer;
  PetscBool  	 flg;
  Vec            TempVec;
  PetscScalar    TempScalar;

  nev = this->num_eigen;


  MatGetVecs(stiff_complex,&eigvalreal,PETSC_NULL);
  MatGetVecs(stiff_complex,&eigvalimag,PETSC_NULL);
  MatGetVecs(stiff_complex,&eigenerror,PETSC_NULL);
  
  ierr = EPSCreate(PETSC_COMM_WORLD,&eps);CHKERRQ(ierr);
  ierr = EPSSetWhichEigenpairs(eps,EPS_SMALLEST_REAL); CHKERRQ(ierr);
  ierr = EPSSetTolerances(eps,LINEAR_SYS_REL_TOL,PETSC_DECIDE);CHKERRQ(ierr);
  ierr = EPSSetDimensions(eps,nev,PETSC_DECIDE,PETSC_DECIDE);CHKERRQ(ierr);
  ierr = EPSSetOperators(eps,stiff_complex,mass_complex);CHKERRQ(ierr);
  ierr = EPSSetFromOptions(eps);CHKERRQ(ierr);


  ierr = PetscPrintf(PETSC_COMM_WORLD,"Starting eigensolve  ...\n");CHKERRQ(ierr);
  ierr = EPSSolve(eps);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"...done with eigensolve.\n");CHKERRQ(ierr);

  ierr = EPSGetIterationNumber(eps, &its);CHKERRQ(ierr);
  ierr = EPSGetOperationCounters(eps,PETSC_NULL,PETSC_NULL,&lits);CHKERRQ(ierr);
  ierr = EPSGetType(eps,&epstype);CHKERRQ(ierr);
  ierr = EPSGetDimensions(eps,&nev,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
  ierr = EPSGetConverged(eps,&nconv);CHKERRQ(ierr);
  ierr = EPSView(eps,PETSC_VIEWER_STDOUT_SELF);CHKERRQ(ierr);
  ierr = EPSGetTolerances(eps,&tol,&maxit);CHKERRQ(ierr);


  for ( int ieigen = 0; ieigen < min(nconv,this->num_eigen); ieigen++){
    ierr = EPSGetEigenpair(eps,ieigen,&kr,&ki,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
    ierr = EPSComputeRelativeError(eps,ieigen,&kerror);
    ierr = VecSetValue(eigvalreal,ieigen,kr,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(eigvalimag,ieigen,ki,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(eigenerror,ieigen,kerror,INSERT_VALUES);CHKERRQ(ierr);
    
  }


  return 0;
}
//========================================================================================
int slepc_system::assemblecomplexmat(Mat *matreal, Mat *matimag, Mat *matcomplex){

  for ( int irow = 0 ; irow < this->num_equations_real; irow++){
    PetscInt rowcommon     = irow;
    PetscInt rowlower      = rowcommon + num_equations_real;
    PetscInt rownnzreal, rownnzimag;

    const PetscInt *colidxreal,*colidximag;
    const PetscScalar *colvaluesreal, *colvaluesimag;

    MatGetRow(*matreal,rowcommon,&rownnzreal,&colidxreal,&colvaluesreal);
    MatGetRow(*matimag,rowcommon,&rownnzimag,&colidximag,&colvaluesimag);

    PetscInt *colidxreallower     = new PetscInt[rownnzreal];
    PetscInt *colidximagupper     = new PetscInt[rownnzimag];
    PetscScalar *colvaluesimagneg = new PetscScalar[rownnzimag];
    
    
    for ( int icol = 0; icol < rownnzreal; icol++){
      colidxreallower[icol] = colidxreal[icol] + num_equations_real;
    }

    for ( int icol = 0; icol < rownnzimag; icol++){
      colidximagupper[icol]  =  colidximag[icol] + num_equations_real;
      colvaluesimagneg[icol] = -colvaluesimag[icol]; 
    }

    // real stiffness  - upper then lower
    MatSetValues(*matcomplex,1,&rowcommon,rownnzreal,colidxreal,colvaluesreal,INSERT_VALUES);
    MatSetValues(*matcomplex,1,&rowlower,rownnzreal,colidxreallower,colvaluesreal,INSERT_VALUES);

    // imaginary parts - upper, then lower
    MatSetValues(*matcomplex,1,&rowcommon,rownnzimag,colidximagupper,colvaluesimagneg,INSERT_VALUES);
    MatSetValues(*matcomplex,1,&rowlower,rownnzimag,colidximag,colvaluesimag,INSERT_VALUES);
    
    MatRestoreRow(*matreal,rowcommon,&rownnzreal,&colidxreal,&colvaluesreal);
    MatRestoreRow(*matimag,rowcommon,&rownnzimag,&colidximag,&colvaluesimag);

    delete [] colidxreallower;
    delete [] colidximagupper;

  } // irow

}

//========================================================================================
void slepc_system::getsol(double **solvecin){
  PetscScalar *avec;
  VecGetArray(solvec,&avec);
  for(int ieqn = 0;ieqn<this->num_equations;ieqn++){
    solvecin[ieqn][which_field] = avec[ieqn];
  }
  // this frees up memory...
  VecRestoreArray(solvec,&avec); 
}
//========================================================================================
int slepc_system::geteigenvalue(double *kreal, double *kimag, double *eigerror, int ieigen){
  PetscScalar *krealvec,*kimagvec,*kerrorvec;
  VecGetArray(eigvalreal,&krealvec);
  VecGetArray(eigvalimag,&kimagvec);
  VecGetArray(eigenerror,&kerrorvec);

  *kreal    = krealvec[ieigen];
  *kimag    = kimagvec[ieigen];
  *eigerror = kerrorvec[ieigen];
  

  VecRestoreArray(eigvalreal,&krealvec);
  VecRestoreArray(eigvalimag,&kimagvec);
  VecRestoreArray(eigenerror,&kerrorvec);

  return 0;
}
//==========================================================================================
int slepc_system::geteigenvector(double **kreal,double **kimag,int ieigen){
  PetscScalar *epsvecreal, *epsvecimag;
  Vec RightVecReal, RightVecImag;

  MatGetVecs(stiff,&RightVecReal,PETSC_NULL);
  MatGetVecs(stiff,&RightVecImag,PETSC_NULL);


  ierr = EPSGetEigenvector(eps,ieigen,RightVecReal,RightVecImag);CHKERRQ(ierr);
  VecGetArray(RightVecReal,&epsvecreal); 
  VecGetArray(RightVecImag,&epsvecimag); 

  for(int ieqn = 0;ieqn<this->num_equations;ieqn++){
    kreal[ieqn][which_field] = epsvecreal[ieqn];
    kimag[ieqn][which_field] = epsvecimag[ieqn];
  }

  ierr = VecRestoreArray(RightVecReal,&epsvecreal); CHKERRQ(ierr);
  ierr = VecRestoreArray(RightVecImag,&epsvecimag); CHKERRQ(ierr);

  // Destroy vecs
  VecDestroy(&RightVecReal);
  VecDestroy(&RightVecImag);
  return 0;
}
//==============================================================================================
int slepc_system::geteigenvectordispersion(double **kreal, double **kimag,int ieigen){
  // gets the ieigen-th eigenvector for the dispersion system
  // this assumes that u = [ur;ui]
  PetscScalar *epsvec;
  Vec RightVec;
  MatGetVecs(stiff_complex,&RightVec,PETSC_NULL);
  ierr = EPSGetEigenvector(eps,ieigen,RightVec,PETSC_NULL);CHKERRQ(ierr);
  VecGetArray(RightVec,&epsvec);
  
  int ictr = 0;
  for ( int ieqn = 0; ieqn < this->num_equations_real; ieqn++){
    kreal[ictr][this->which_field]=epsvec[ieqn];
    ictr++;
  } 

  ictr = 0;
  for ( int ieqn = this->num_equations_real; ieqn < this->num_equations_imag; ieqn++){
    kimag[ictr][this->which_field]=epsvec[ieqn];
    ictr++;
  }

  ierr = VecRestoreArray(RightVec,&epsvec); CHKERRQ(ierr);
  ierr = VecDestroy(&RightVec); CHKERRQ(ierr);

  return 0;
}
//===================================================================================
void slepc_system::setsolvetype(int stypearg){
  solvetype = stypearg;
}
//====================================================================================
int slepc_system::getconveigen(){
  return this->nconv;
}
//===================================================================================

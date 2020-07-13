#ifndef HAVE_SLEPC_QUAD_SYSTEM_HPP
#include "slepc_quad_system.hpp"
#endif

// constructor and destructor
//==================================================================================================================
slepc_quad_system::slepc_quad_system(){}
//===================================================================================================================
slepc_quad_system::slepc_quad_system(int num_equations_arg, int which_field, int inewton, int iload,int iprintflag){
  this->num_equations = num_equations_arg;
  this->which_field   = which_field;
  this->inewton       = inewton;
  this->iload         = iload;
  this->printflag     = iprintflag;
}

// ==============================================================================================================
slepc_quad_system::slepc_quad_system(equation_initdata eqninitdata){
  this->num_equations      = eqninitdata.num_equations;
  this->num_equations_real = eqninitdata.num_equations_real;
  this->num_equations_imag = eqninitdata.num_equations_imag;
  this->max_nnz_per_row    = eqninitdata.max_nnz_per_row;

  this->which_field = eqninitdata.which_field;
  this->inewton     = eqninitdata.inewton;
  this->iload       = eqninitdata.iload;
  this->printflag   = eqninitdata.printflag;
  this->which_pt    = eqninitdata.which_pt;
  this->which_seg   = eqninitdata.which_seg;
}
// ===========================================================================================================
slepc_quad_system::~slepc_quad_system(){
  ierr = MatDestroy(&this->KK);   CHKERRXX(ierr);
  ierr = MatDestroy(&this->CC);   CHKERRXX(ierr);
  ierr = MatDestroy(&this->MM);   CHKERRXX(ierr);
  ierr = MatDestroy(&this->KKr1); CHKERRXX(ierr);
  ierr = MatDestroy(&this->KKr2); CHKERRXX(ierr);
  ierr = MatDestroy(&this->KKi);  CHKERRXX(ierr);
  ierr = MatDestroy(&this->MMr);  CHKERRXX(ierr);
  ierr = VecDestroy(&eigvalreal); CHKERRXX(ierr);
  ierr = VecDestroy(&eigvalimag); CHKERRXX(ierr);
  ierr = VecDestroy(&eigenerror); CHKERRXX(ierr);
  ierr = QEPDestroy(&qep);        CHKERRXX(ierr);
}
//==========================================================================================================
void slepc_quad_system::createcsrkr1(int *row_ptr,int *col_ptr,double *value_arr){
  MatCreateSeqAIJWithArrays(PETSC_COMM_SELF,num_equations_real,num_equations_real,row_ptr,col_ptr,value_arr,&KKr1);
}
//==========================================================================================================
void slepc_quad_system::createcsrkr2(int *row_ptr,int *col_ptr,double *value_arr){
  MatCreateSeqAIJWithArrays(PETSC_COMM_SELF,num_equations_real,num_equations_real,row_ptr,col_ptr,value_arr,&KKr2);
}
//==========================================================================================================
void slepc_quad_system::createcsrki(int *row_ptr,int *col_ptr,double *value_arr){
  MatCreateSeqAIJWithArrays(PETSC_COMM_SELF,num_equations_imag,num_equations_imag,row_ptr,col_ptr,value_arr,&KKi);
}
//==========================================================================================================
void slepc_quad_system::createcsrcreal(int *row_ptr,int *col_ptr,double *value_arr){
  MatCreateSeqAIJWithArrays(PETSC_COMM_SELF,num_equations_real,num_equations_real,row_ptr,col_ptr,value_arr,&CCr);
}
//==========================================================================================================
void slepc_quad_system::createcsrcimag(int *row_ptr,int *col_ptr,double *value_arr){
  MatCreateSeqAIJWithArrays(PETSC_COMM_SELF,num_equations_imag,num_equations_imag,row_ptr,col_ptr,value_arr,&CCi);
}
//==========================================================================================================
void slepc_quad_system::createcsrmr(int *row_ptr,int *col_ptr,double *value_arr){
  MatCreateSeqAIJWithArrays(PETSC_COMM_SELF,num_equations_real,num_equations_real,row_ptr,col_ptr,value_arr,&MMr);
}
//==========================================================================================================
void slepc_quad_system::createcsrmi(int *row_ptr,int *col_ptr,double *value_arr){
  MatCreateSeqAIJWithArrays(PETSC_COMM_SELF,num_equations_real,num_equations_real,row_ptr,col_ptr,value_arr,&MMi);
}
//==========================================================================================================
int slepc_quad_system::solve(){
}
//===========================================================================================================
int slepc_quad_system::solvecomplex(){
  PetscInt nrows = num_equations_real + num_equations_imag;
  PetscInt ncols = num_equations_real + num_equations_imag;
  PetscInt nnzcomplex = 2*this->max_nnz_per_row;
  PetscInt its,nev,maxit;
  PetscReal tol,kerror;
  PetscScalar kr,ki;


  // (k^2)Px  + (k)Qx + Rx = 0   ==> P=M , Q=C; R=K  

  MatCreateSeqAIJ(PETSC_COMM_SELF,nrows,ncols,nnzcomplex,PETSC_NULL,&MM);
  MatCreateSeqAIJ(PETSC_COMM_SELF,nrows,ncols,nnzcomplex,PETSC_NULL,&CC);
  MatCreateSeqAIJ(PETSC_COMM_SELF,nrows,ncols,nnzcomplex,PETSC_NULL,&KK);
  
  MatGetVecs(MM,&eigvalreal,PETSC_NULL);
  MatGetVecs(MM,&eigvalimag,PETSC_NULL);
  MatGetVecs(MM,&eigenerror,PETSC_NULL);

  assemblecomplexPP();    assemblecomplexQQ();   assemblecomplexRR();

  if ( ( this->inewton == 0 ) && ( this->iload == 0 )  && ( this->printflag == 1 ) ) {
    {
      // view Kr1 
      char filename[100];
      sprintf(filename,"matkr1.m",this->which_field,this->which_seg,this->which_pt);
      PetscViewerASCIIOpen(PETSC_COMM_SELF,filename,&view_mat);
      PetscViewerSetFormat(view_mat,PETSC_VIEWER_ASCII_MATLAB);
      PetscObjectSetName((PetscObject)(KKr1),"KKr1");
      MatView(KKr1,view_mat);
      PetscViewerDestroy(&view_mat);
    }

    {
      // view M_r
      char filename[100];
      sprintf(filename,"matmmr.m",this->which_field,this->which_seg,this->which_pt);
      PetscViewerASCIIOpen(PETSC_COMM_SELF,filename,&view_mat);
      PetscViewerSetFormat(view_mat,PETSC_VIEWER_ASCII_MATLAB);
      PetscObjectSetName((PetscObject)(MMr),"MMr");
      MatView(MMr,view_mat);
      PetscViewerDestroy(&view_mat);
    }

    {
      // view KK = [ Kr1-w^2 Mr              0 ]
      //           [ 0           Kr1 - w^2 Mr  ]   
      char filename[100];
      sprintf(filename,"matkfield%dseg%dpt%d.m",this->which_field,this->which_seg,this->which_pt);
      PetscViewerASCIIOpen(PETSC_COMM_SELF,filename,&view_mat);
      PetscViewerSetFormat(view_mat,PETSC_VIEWER_ASCII_MATLAB);
      PetscObjectSetName((PetscObject)(KK),"KK");
      MatView(KK,view_mat);
      PetscViewerDestroy(&view_mat);
    }

    {
      // view MM
      char filename[100];
      sprintf(filename,"matm.m",this->which_field,this->which_seg,this->which_pt);
      PetscViewerASCIIOpen(PETSC_COMM_SELF,filename,&view_mat);
      PetscViewerSetFormat(view_mat,PETSC_VIEWER_ASCII_MATLAB);
      PetscObjectSetName((PetscObject)(MM),"MM");
      MatView(MM,view_mat);
      PetscViewerDestroy(&view_mat);
    }
    {
      // view CC
      char filename[100];
      sprintf(filename,"matc.m",this->which_field,this->which_seg,this->which_pt);
      PetscViewerASCIIOpen(PETSC_COMM_SELF,filename,&view_mat);
      PetscViewerSetFormat(view_mat,PETSC_VIEWER_ASCII_MATLAB);
      PetscObjectSetName((PetscObject)(CC),"CC");
      MatView(CC,view_mat);
      PetscViewerDestroy(&view_mat);
    }
  }

  // create the QEP system
  ierr = QEPCreate(PETSC_COMM_WORLD,&qep); CHKERRQ(ierr);
  ierr = QEPSetOperators(qep,MM,CC,KK);    CHKERRQ(ierr);
  ierr = QEPSetFromOptions(qep);           CHKERRQ(ierr);
  
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Calling QEPSolve\n");CHKERRQ(ierr);
  ierr = QEPSolve(qep);CHKERRQ(ierr);
  ierr = QEPGetIterationNumber(qep,&its);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," Number of iterations of the method: %D\n",its);CHKERRQ(ierr);
  
  ierr = QEPGetType(qep,&type);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," Solution method: %s\n\n",type);CHKERRQ(ierr);
  ierr = QEPGetDimensions(qep,&nev,NULL,NULL);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," Number of requested eigenvalues: %D\n",nev);CHKERRQ(ierr);
  ierr = QEPGetTolerances(qep,&tol,&maxit);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," Stopping condition: tol=%.4G, maxit=%D\n",tol,maxit);CHKERRQ(ierr);
  
  ierr = QEPView(qep,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  ierr = QEPGetConverged(qep,&nconv);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," Number of converged approximate eigenpairs: %d\n\n",nconv);CHKERRQ(ierr);
  ierr = QEPPrintSolution(qep,NULL);CHKERRQ(ierr);

  for ( int ieigen = 0; ieigen < min(nconv,this->num_eigen); ieigen++){
    ierr = QEPGetEigenpair(qep,ieigen,&kr,&ki,PETSC_NULL,PETSC_NULL);CHKERRQ(ierr);
    ierr = QEPComputeRelativeError(qep,ieigen,&kerror);
    ierr = VecSetValue(eigvalreal,ieigen,kr,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(eigvalimag,ieigen,ki,INSERT_VALUES);CHKERRQ(ierr);
    ierr = VecSetValue(eigenerror,ieigen,kerror,INSERT_VALUES);CHKERRQ(ierr);
  }
  
  cout<<"Destroyed Quad System"<<endl;

}
//===========================================================================================================
int slepc_quad_system::assemblecomplexPP(){

  // (k^2)Px  + (k)Qx + Rx = 0   ==> P=M , Q=C; R=K  
  // P = [ Kr2  0  ] 
  //     [ 0    Kr2]

  for ( int irow = 0; irow < this->num_equations_real; irow++){
    PetscInt rowupper  = irow;
    PetscInt rowlower  = irow + num_equations_real;
    PetscInt rownnz;
    const PetscInt *colidxupper;
    const PetscScalar *colvalues;
    MatGetRow(KKr2,irow,&rownnz,&colidxupper,&colvalues);
   
    PetscInt *colidxlower = new PetscInt[rownnz];
    for ( int icol = 0; icol < rownnz ; icol++){
      // shift lower col indices
      colidxlower[icol] = colidxupper[icol] + (num_equations_real);
    }

    MatSetValues(MM,1,&rowupper, rownnz,colidxupper,colvalues,INSERT_VALUES);
    MatSetValues(MM,1,&rowlower, rownnz,colidxlower,colvalues,INSERT_VALUES);
    
    MatRestoreRow(KKr2,rowupper,&rownnz,&colidxupper,&colvalues);
    delete [] colidxlower;
  } // irow

  ierr=MatAssemblyBegin(MM,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr=MatAssemblyEnd(MM,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  return 0;
}
//===========================================================================================================
int slepc_quad_system::assemblecomplexRR(){

  // (k^2)Px  + (k)Qx + Rx = 0   ==> P=M , Q=C; R=K  
  //  R = [ Kr1-w^2M_r            0     ]
  //      [ 0              Kr1 - w^2M_r ]

  // allocate temporary matrix 
  Mat AA[2],RR;
  ierr = MatDuplicate(KKr1,MAT_COPY_VALUES,&AA[0]); CHKERRQ(ierr);
  ierr = MatDuplicate(MMr ,MAT_COPY_VALUES,&AA[1]); CHKERRQ(ierr);
  ierr = MatScale(AA[1],-ww*ww);                    CHKERRQ(ierr);
  ierr = MatCreateComposite(PETSC_COMM_WORLD,2,AA,&RR); CHKERRQ(ierr);
  ierr = MatCompositeMerge(RR);                     CHKERRQ(ierr);

   for ( int irow = 0; irow < this->num_equations_real; irow++){
    PetscInt rowupper  = irow;
    PetscInt rowlower  = irow + num_equations_real;
    PetscInt rownnz;
    const PetscInt *colidxupper;
    const PetscScalar *colvalues;
    MatGetRow(RR,irow,&rownnz,&colidxupper,&colvalues);

    PetscInt *colidxlower = new PetscInt[rownnz];
    for ( int icol = 0; icol < rownnz ; icol++){
      // shift lower col indices
      colidxlower[icol] = colidxupper[icol] + num_equations_real;
    }

    MatSetValues(KK,1,&irow,    rownnz,colidxupper,colvalues,INSERT_VALUES);
    MatSetValues(KK,1,&rowlower,rownnz,colidxlower,colvalues,INSERT_VALUES);

    MatRestoreRow(RR,irow,&rownnz,&colidxupper,&colvalues);
    delete [] colidxlower;
   }

  ierr=MatAssemblyBegin(KK,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr=MatAssemblyEnd(KK,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr=MatDestroy(&AA[0]); CHKERRQ(ierr);
  ierr=MatDestroy(&AA[1]); CHKERRQ(ierr);
  ierr=MatDestroy(&RR); CHKERRQ(ierr);

  return 0;
}
//===========================================================================================================
int slepc_quad_system::assemblecomplexQQ(){
  // (k^2)Px  + (k)Qx + Rx = 0   ==> P=M , Q=C; R=K  

  // Q = [ 0   -Ki ] 
  //     [ Ki   0  ]

  for ( int irow = 0; irow < this->num_equations_real; irow++){
    PetscInt rowupper  = irow;
    PetscInt rowlower  = irow + num_equations_real;
    PetscInt rownnz;
    const PetscInt *colidxlower;
    const PetscScalar *colvalues;
    MatGetRow(KKi,irow,&rownnz,&colidxlower,&colvalues);
   
    PetscInt    *colidxupper  = new PetscInt[rownnz];
    PetscScalar *colvaluesneg = new PetscScalar[rownnz];

    for ( int icol = 0; icol < rownnz ; icol++){
      // shift upper col indices
      colidxupper[icol]  =  colidxlower[icol] + (num_equations_real);
      colvaluesneg[icol] = -colvalues[icol];
    }

    MatSetValues(CC,1,&rowupper, rownnz,colidxupper,colvaluesneg,INSERT_VALUES);
    MatSetValues(CC,1,&rowlower, rownnz,colidxlower,colvalues,INSERT_VALUES);
    
    MatRestoreRow(KKr2,rowupper,&rownnz,&colidxlower,&colvalues);

    delete [] colidxupper;
    delete [] colvaluesneg;
 
  } // irow

  ierr=MatAssemblyBegin(CC,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr=MatAssemblyEnd(CC,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

  return 0;
}
//===========================================================================================================
void slepc_quad_system::getsol(double **solvecin){
  PetscScalar *avec;
  VecGetArray(solvec,&avec);
  for(int ieqn = 0;ieqn<this->num_equations;ieqn++){
    solvecin[ieqn][which_field] = avec[ieqn];
  }
  // this frees up memory...
  VecRestoreArray(solvec,&avec); 
}
//============================================================================================================
int slepc_quad_system::geteigenvalue(double *kreal, double *kimag, double *eigerror, int ieigen){

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
//===========================================================================================================
int slepc_quad_system::geteigenvectordispersion(double **kreal, double **kimag,int ieigen){
  // gets the ieigen-th eigenvector for the dispersion system
  // this assumes that u = [ur;ui]
  PetscScalar *qepvec;
  Vec RightVec,RightVecImag;
  MatGetVecs(MM,&RightVec,PETSC_NULL);


  ierr = QEPGetEigenpair(qep,ieigen,PETSC_NULL,PETSC_NULL,RightVec,PETSC_NULL);CHKERRQ(ierr);
  VecGetArray(RightVec,&qepvec);
  
  int ictr = 0;
  for ( int ieqn = 0; ieqn < this->num_equations_real; ieqn++){
    kreal[ictr][this->which_field]=qepvec[ieqn];
    ictr++;
  } 

  ictr = 0;
  for ( int ieqn = this->num_equations_real; ieqn < this->num_equations_imag; ieqn++){
    kimag[ictr][this->which_field]=qepvec[ieqn];
    ictr++;
  }

  ierr = VecRestoreArray(RightVec,&qepvec); CHKERRQ(ierr);
  ierr = VecDestroy(&RightVec); CHKERRQ(ierr);

  return 0;
}
//================================================================================================================
int slepc_quad_system::getconveigen(){
  return this->nconv;
}
//=====================================================================================================

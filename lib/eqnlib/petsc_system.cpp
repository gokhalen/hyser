#ifndef HAVE_PETSC_SYSTEM_HPP
#include "petsc_system.hpp"
#endif


petsc_system::petsc_system(){}

petsc_system::petsc_system(int num_equations_arg, int which_field, int inewton, int iload,int iprintflag){
  
  this->num_equations = num_equations_arg;
  this->which_field   = which_field;
  this->inewton       = inewton;
  this->iload         = iload;
  this->printflag     = iprintflag;
}

petsc_system::petsc_system(equation_initdata eqninitdata){
  
  this->num_equations = eqninitdata.num_equations;
  this->which_field   = eqninitdata.which_field;
  this->inewton       = eqninitdata.inewton;
  this->iload         = eqninitdata.iload;
  this->printflag     = eqninitdata.printflag;
}


petsc_system::~petsc_system(){
  ierr = MatDestroy(&this->stiff);CHKERRXX(ierr);
  ierr = VecDestroy(&this->rhsvec);CHKERRXX(ierr);
  ierr = VecDestroy(&this->solvec);CHKERRXX(ierr);
}

void petsc_system::createcsrk(int *row_ptr,int *col_ptr,double *value_arr){

  char filename[100];
  sprintf(filename,"matk%d%d%d.m",this->which_field,this->inewton,this->iload);

  ierr = MatCreateSeqAIJWithArrays(PETSC_COMM_SELF,num_equations,num_equations,row_ptr,col_ptr,value_arr,&stiff);CHKERRXX(ierr);
   
  if(this->printflag == 1){
    ierr = PetscViewerASCIIOpen(PETSC_COMM_SELF,filename,&view_mat);CHKERRXX(ierr);
    ierr = PetscViewerSetFormat(view_mat,PETSC_VIEWER_ASCII_MATLAB);CHKERRXX(ierr);
    ierr = MatView(stiff,view_mat);CHKERRXX(ierr);
    ierr = PetscViewerDestroy(&view_mat);CHKERRXX(ierr);
  }


}

void petsc_system::createrhs(double **value_arr){

  // we do not use VecCreateWithArrays for the rhs
  // This is because we need to pass a pointer to a contiguous array
  // With the current array declaration [num_equations][nfields] this is not possible
  // It is possible if [nfields][num_equations]. However, size/time is not an issue here
  // Therefore copying the array over is justified.
  
  
  ierr  = VecCreate(PETSC_COMM_WORLD,&rhsvec);CHKERRXX(ierr);
  ierr  = VecSetSizes(rhsvec,PETSC_DECIDE,num_equations);CHKERRXX(ierr);
  ierr  = VecSetFromOptions(rhsvec);CHKERRXX(ierr);

  for (int ieqn = 0; ieqn<this->num_equations;ieqn++){
    ierr = VecSetValue(rhsvec,ieqn,value_arr[ieqn][which_field], INSERT_VALUES);CHKERRXX(ierr);
  }
  
  VecAssemblyBegin(rhsvec);
  VecAssemblyEnd(rhsvec);
  
  ierr = VecDuplicate(rhsvec,&solvec);CHKERRXX(ierr);
  ierr = PetscViewerASCIIOpen(PETSC_COMM_SELF,"rhs.out",&view_rhs);CHKERRXX(ierr);
  ierr = VecView(rhsvec,view_rhs);CHKERRXX(ierr);
  ierr = PetscViewerDestroy(&view_rhs);CHKERRXX(ierr);
  
}

int petsc_system::solve(){
  KSP      ksp;      // linear   solver context
  PC       pc;       // linear   solver context
  PetscReal norm;    // norm of  solution error
  Vec uu;
  //PetscReal neg_one = -1.0;
  PetscInt its;
  PetscBool SystemScale;

  // Define a scaled matrix to solve
  Mat MatScal;
  Vec RhsScal,ScalFacs,ScalFacsInv;
  PetscScalar *alpha,*rhsvals; 

  SystemScale = PETSC_FALSE;
  if ( SystemScale ) {
    ierr = PetscPrintf(PETSC_COMM_WORLD,"Scaling Linear System... ");CHKERRQ(ierr);
    ierr = MatDuplicate(stiff,MAT_COPY_VALUES,&MatScal);CHKERRQ(ierr);
    ierr = VecDuplicate(rhsvec,&RhsScal); CHKERRQ(ierr);
    ierr = VecCopy(rhsvec,RhsScal); CHKERRQ(ierr);
    ierr = VecDuplicate(rhsvec,&ScalFacs);CHKERRQ(ierr);
    ierr = VecDuplicate(rhsvec,&ScalFacsInv);CHKERRQ(ierr);
    ierr = MatGetDiagonal(MatScal,ScalFacs);CHKERRQ(ierr);
    
    ierr = VecGetArray(ScalFacs,&alpha);CHKERRQ(ierr);
    ierr = VecGetArray(rhsvec,&rhsvals);CHKERRQ(ierr);

    for ( int ieqn = 0 ; ieqn < this->num_equations ;ieqn++){
      ierr = VecSetValue(RhsScal,ieqn,rhsvals[ieqn]/sqrt(fabs(alpha[ieqn])),INSERT_VALUES);CHKERRQ(ierr);
      ierr = VecSetValue(ScalFacsInv,ieqn,1.0/sqrt(fabs(alpha[ieqn])),INSERT_VALUES);CHKERRQ(ierr);
    }
    
    VecAssemblyBegin(RhsScal);
    VecAssemblyEnd(RhsScal);

    ierr = VecRestoreArray(rhsvec,&rhsvals);CHKERRQ(ierr);
    
    
    ierr = MatDiagonalScale(MatScal,ScalFacsInv,ScalFacsInv);CHKERRQ(ierr);

    
  }

  ierr = VecDuplicate(solvec,&uu); CHKERRQ(ierr);

  ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);CHKERRQ(ierr);
  ierr = KSPSetType(ksp,KSPCG);CHKERRQ(ierr);

  if ( SystemScale ){
    ierr = KSPSetOperators(ksp,MatScal,MatScal,DIFFERENT_NONZERO_PATTERN);CHKERRQ(ierr);
  }
  else{
    ierr = KSPSetOperators(ksp,stiff,stiff,DIFFERENT_NONZERO_PATTERN);CHKERRQ(ierr);
  }
  
  ierr = KSPGetPC(ksp,&pc);CHKERRQ(ierr);
  ierr = PCSetType(pc,PCLU);CHKERRQ(ierr);
  ierr = PCFactorSetLevels(pc,4);CHKERRQ(ierr);
  ierr = KSPSetTolerances(ksp,1.e-12,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT); CHKERRQ(ierr);
  // override the above options if options are set on the command line
  ierr = KSPSetFromOptions(ksp); CHKERRQ(ierr);

  if ( SystemScale ){
    ierr = KSPSolve(ksp,RhsScal,solvec);CHKERRQ(ierr);
    PetscScalar *solvals;
    Vec solvec2;
    ierr = VecGetArray(solvec,&solvals);CHKERRQ(ierr);
    ierr = VecDuplicate(solvec,&solvec2);CHKERRQ(ierr);

    for ( int ieqn = 0 ; ieqn < this->num_equations ;ieqn++){
      ierr = VecSetValue(solvec2,ieqn,solvals[ieqn]/sqrt(fabs(alpha[ieqn])),INSERT_VALUES);CHKERRQ(ierr);
    }

    ierr = VecRestoreArray(solvec,&solvals);CHKERRQ(ierr);
    ierr = VecAssemblyBegin(solvec2);CHKERRQ(ierr);
    ierr = VecAssemblyEnd(solvec2);CHKERRQ(ierr);
    ierr = VecCopy(solvec2,solvec);CHKERRQ(ierr);
    ierr = VecDestroy(&solvec2);CHKERRQ(ierr);

  }
  else{
    ierr = KSPSolve(ksp,rhsvec,solvec); CHKERRQ(ierr);
  }

  ierr = KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
  ierr = KSPGetIterationNumber(ksp,&its); CHKERRQ(ierr);
  ierr = KSPGetResidualNorm(ksp,&norm); CHKERRQ(ierr);
  //MatMult(stiff,solvec,uu);
  
  

  ierr = PetscPrintf(PETSC_COMM_WORLD,"Norm of error %G, Iterations %D\n",
	      norm,its); CHKERRQ(ierr);

  ierr = PetscViewerASCIIOpen(PETSC_COMM_SELF,"solpetsc.out",&view_sol); CHKERRQ(ierr);
  ierr = VecView(solvec,view_sol); CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&view_sol); CHKERRQ(ierr);
  
  if ( SystemScale) {
  ierr = VecRestoreArray(ScalFacs,&alpha);CHKERRQ(ierr);
  }

  ierr = VecDestroy(&uu); CHKERRQ(ierr);
  ierr = KSPDestroy(&ksp); CHKERRQ(ierr);

  return 0;

}

void petsc_system::getsol(double **solvecin){
  PetscScalar *avec;
  ierr = VecGetArray(solvec,&avec);CHKERRXX(ierr);
  for(int ieqn = 0;ieqn<this->num_equations;ieqn++){
    solvecin[ieqn][which_field] = avec[ieqn];
  }
  // this frees up memory...
  ierr = VecRestoreArray(solvec,&avec); CHKERRXX(ierr);
}

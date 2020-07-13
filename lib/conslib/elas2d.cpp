/*---------------------------------------------------------------------------------------
 * Nachiket Gokhale gokhalen@bu.edu
 * 
 * Calculates material stiffness C_IJKL
 * Calculates second Piola-Kirchhoff (2PK)
 * Calculates gradient of 2PK wrt material parameters
 --------------------------------------------------------------------------------------*/  

#ifndef HAVE_CONSLIB_H
#include "conslib.h"
#endif

#ifndef HAVE_UTILS_HPP
#include "../../utils/utils.hpp"
#endif


#include <math.h>

int elas2d(int lawnumber, conselas23d *elasptr){

  // for this law we need rightcauchygreeninverse transpose
    double **rcgreeninvtpose;
    
    
    rcgreeninvtpose = allocmat(elasptr->ndime,elasptr->ndime,rcgreeninvtpose);
    invertsmall(elasptr->ndime, elasptr->rcgreen, rcgreeninvtpose);
    tposeinplace(elasptr->ndime,rcgreeninvtpose);

    zero4D(elasptr->ndime,elasptr->ndime,elasptr->ndime,elasptr->ndime,elasptr->aijklptr);
    zero4D(elasptr->ndime,elasptr->ndime,elasptr->ndime,elasptr->ndime,elasptr->sptan);

  if(lawnumber == 1){
    
    // This is Simo and Hughes 
    // calculate C_IJKL
    double lambda = elasptr->stfparam[0];
    double mu     = elasptr->stfparam[1];
    double dett   = determinant(elasptr->ndime, elasptr->rcgreen);

    for(int ii = 0;ii<elasptr->ndime;ii++){
      for(int jj = 0;jj<elasptr->ndime;jj++){
	for(int kk = 0;kk<elasptr->ndime;kk++){
	  for(int ll = 0;ll<elasptr->ndime;ll++){
	    
	    double symmfourth1  = rcgreeninvtpose[ii][kk]*rcgreeninvtpose[jj][ll];
	    double symmfourth2  = rcgreeninvtpose[jj][kk]*rcgreeninvtpose[ii][ll];
	    double symmfourth   = 0.5*(symmfourth1+symmfourth2);
	    double Cinv_ijCinv_kl  = rcgreeninvtpose[ii][jj]*rcgreeninvtpose[kk][ll];
	    
	      
	    double first  = (mu-0.5*lambda*(dett-1))*symmfourth;
	    double second = 0.5*lambda*dett*Cinv_ijCinv_kl;
	    
	    elasptr->cijklptr[ii][jj][kk][ll] = 2*(first + second);
	  }
	}
      }
    }
    
     
    // calculate second piola kirchhoff
     for(int irow = 0; irow<elasptr->ndime;irow++){
	for(int icol = 0;icol<elasptr->ndime;icol++){
	  
	  double right  = rcgreeninvtpose[irow][icol];
	  double secpiola  = lambda*(dett - 1)*0.5*right + mu*(delta(irow,icol) -right);
	  
	  elasptr->secpk[irow][icol] = secpiola;
	}
      }


    // calculate the gradient contribution
    // check if we are required to calculate it
     if(elasptr->grdflag == 1){
       // the purpose is to set pk2grad

       if(elasptr->grdparm == 0){
	 // taking gradient with first stiffness parameter, \mu
	 
	 for(int ii = 0;ii<elasptr->ndime;ii++){
	   for(int jj =0;jj<elasptr->ndime;jj++){
	     // note: dett is the determinant of C. det(C)= det(J^2)
	     elasptr->pk2grad[ii][jj] = 0.5*(dett-1.0)*rcgreeninvtpose[ii][jj];
	   }
	 }
       }

       if(elasptr->grdparm == 1){
	 // taking gradient with second stiffness paramter, \mu
	 for(int ii = 0; ii<elasptr->ndime;ii++){
	   for(int jj =0; jj<elasptr->ndime;jj++){
	     elasptr->pk2grad[ii][jj] = (delta(ii,jj)-rcgreeninvtpose[ii][jj]);
	   }
	 }
       }
     }
    
     
  }
  

  if(lawnumber == 2){
    // calculate C_IJKL
    double lambda = elasptr->stfparam[0];
    double mu     = elasptr->stfparam[1];
    double dett   = determinant(elasptr->ndime, elasptr->rcgreen);

    for(int ii = 0;ii<elasptr->ndime;ii++){
      for(int jj = 0;jj<elasptr->ndime;jj++){
	for(int kk = 0;kk<elasptr->ndime;kk++){
	  for(int ll = 0;ll<elasptr->ndime;ll++){
	    
	    double symmfourth1  = rcgreeninvtpose[ii][kk]*rcgreeninvtpose[jj][ll];
	    double symmfourth2  = rcgreeninvtpose[jj][kk]*rcgreeninvtpose[ii][ll];
	    double symmfourth   = 0.5*(symmfourth1+symmfourth2);
	    double Cinv_ijCinv_kl  = rcgreeninvtpose[ii][jj]*rcgreeninvtpose[kk][ll];
	    
	      
	    double first  = (mu-0.5*lambda*(dett-1))*symmfourth;
	    double second = 0.5*lambda*dett*Cinv_ijCinv_kl;
	    
	    elasptr->cijklptr[ii][jj][kk][ll] = 2*(first) + 2*(second);
	  }
	}
      }
    }



     
    // calculate second piola kirchhoff
     for(int irow = 0; irow<elasptr->ndime;irow++){
	for(int icol = 0;icol<elasptr->ndime;icol++){
	  
	  double right  = rcgreeninvtpose[irow][icol];
	  double secpiola  = lambda*(dett - 1)*0.5*right + mu*(delta(irow,icol) -right);
	  
	  elasptr->secpk[irow][icol] = secpiola;
	}
      }


    // calculate the gradient contribution
    // check if we are required to calculate it
     if(elasptr->grdflag == 1){
       // the purpose is to set pk2grad

       if(elasptr->grdparm == 0){
	 // taking gradient with first stiffness parameter, \mu
	 
	 for(int ii = 0;ii<elasptr->ndime;ii++){
	   for(int jj =0;jj<elasptr->ndime;jj++){
	     // note: dett is the determinant of C. det(C)= det(J^2)
	     elasptr->pk2grad[ii][jj] = 0.5*(dett-1.0)*rcgreeninvtpose[ii][jj];
	   }
	 }
       }

       if(elasptr->grdparm == 1){
	 // taking gradient with second stiffness paramter, \mu
	 for(int ii = 0; ii<elasptr->ndime;ii++){
	   for(int jj =0; jj<elasptr->ndime;jj++){
	     elasptr->pk2grad[ii][jj] = (delta(ii,jj)-rcgreeninvtpose[ii][jj]);
	   }
	 }
       }
     }
    
  } // lawnumber2


  if(lawnumber == 3){
    double gamma = 1.2;
    double beta   = 1000;
    double mu     = elasptr->stfparam[0];
    double alpha  = elasptr->stfparam[1];
    
    double I1   = tracemat(elasptr->ndime,elasptr->rcgreen);
    double I3   = determinant(elasptr->ndime, elasptr->rcgreen);
    
  

    for(int ii = 0;ii<elasptr->ndime;ii++){
      for(int jj = 0;jj<elasptr->ndime;jj++){
	for(int kk = 0;kk<elasptr->ndime;kk++){
	  for(int ll = 0;ll<elasptr->ndime;ll++){
	    
	    double first = mu*gamma*exp(gamma*(I1-elasptr->ndime))*delta(kk,ll)*delta(ii,jj);
	    double second = alpha*delta(kk,ll)*delta(ii,jj) - 0.5*alpha*(delta(ii,kk)*delta(jj,ll) + delta(jj,kk)*delta(ii,ll));

	    double symm1  = 2*rcgreeninvtpose[kk][ll]*rcgreeninvtpose[ii][jj];
	    double symm2  = rcgreeninvtpose[ii][kk]*rcgreeninvtpose[jj][ll] + rcgreeninvtpose[jj][kk]*rcgreeninvtpose[ii][ll];

	    double third  = 0.5*beta*(I3*symm1 - (I3-1)*symm2);

	   	    
	    elasptr->cijklptr[ii][jj][kk][ll] = 4*(first + second + third);

	  }//ll
	}//kk
      }//jj
    } //ii


   // piola kirchoff
   for( int ii = 0; ii < elasptr->ndime; ii ++ ){
     for(int jj = 0; jj < elasptr->ndime; jj++){

       double first  = (mu*exp(gamma*(I1-elasptr->ndime)) + alpha*I1)*delta(ii,jj);
       double second = -0.5*alpha*(elasptr->rcgreen[ii][jj] + elasptr->rcgreen[jj][ii]);
       double third  =  0.5*beta*(I3-1)*(rcgreeninvtpose[ii][jj] + rcgreeninvtpose[jj][ii]);

       elasptr->secpk[ii][jj] = 2*(first + second + third);
     }
   }


   // gradient of the second piola-kirchoff with respect to stiffness parameters

   // check if required to calculate the gradient
   if (elasptr->grdflag == 1){

     // loop over pk2grad matrix entries
     for(int irow = 0;irow < elasptr->ndime; irow++){
       for(int icol = 0;icol < elasptr->ndime;icol++){

	 if(elasptr->grdparm == 0){
	   // taking the gradient with first stiffness parameter: \mu
	   double first = 2*(exp(gamma*(I1-elasptr->ndime)))*delta(irow,icol);
	   elasptr->pk2grad[irow][icol] = first ;
	   //elasptr->pk2grad[irow][icol] = 0;
	 }
	 
	 if(elasptr->grdparm == 1){
	   // taking the gradient with the second stiffness parameter \alpha
	   double first  = I1*delta(irow,icol);
	   double second = elasptr->rcgreen[irow][icol] + elasptr->rcgreen[icol][irow];
	   //elasptr->pk2grad[irow][icol] = 0;
	   elasptr->pk2grad[irow][icol] = 2*(first - 0.5*(second));
	 }

       }//icol
     }//irow
     
   }//gradflg
   
  } // lawnumber3


  // begin lawnumber 4: Veronda Westmann specialised to 2D.
  // alpha = -mu/2
  // we take gradient with respect to both \mu and \gamma
  // mu and gamma are independent
  
  if (lawnumber == 4){
    
    double mu     = elasptr->stfparam[0];
    double alpha  = -mu/2.0;

    double gamma  = elasptr->stfparam[1];
    double beta   = elasptr->stfparam[2];
    
    
    double I1   = tracemat(elasptr->ndime,elasptr->rcgreen);
    double I3   = determinant(elasptr->ndime, elasptr->rcgreen);

    // CIJKL: a small change in the exponential term from the 3D. VerWestmann

    for(int ii = 0; ii < elasptr->ndime; ii++){
      for(int jj = 0; jj < elasptr->ndime; jj++){
	for(int kk = 0; kk < elasptr->ndime; kk++){
	  for(int ll = 0; ll < elasptr->ndime; ll++){
	    
	    double first = mu*gamma*exp(gamma*(I1-2))*delta(kk,ll)*delta(ii,jj);
	    double second = alpha*delta(kk,ll)*delta(ii,jj) - 0.5*alpha*(delta(ii,kk)*delta(jj,ll) + delta(jj,kk)*delta(ii,ll));

	    double symm1  = 2*rcgreeninvtpose[kk][ll]*rcgreeninvtpose[ii][jj];
	    double symm2  = rcgreeninvtpose[ii][kk]*rcgreeninvtpose[jj][ll] + rcgreeninvtpose[jj][kk]*rcgreeninvtpose[ii][ll];

	    double third  = 0.5*beta*(I3*symm1 - (I3-1)*symm2);

	   	    
	    elasptr->cijklptr[ii][jj][kk][ll] = 4*(first + second + third);

	  }//ll
	}//kk
      }//jj
    } //ii
    
    // Second Piola-Kirchhoff
    for( int ii = 0; ii < elasptr->ndime; ii ++ ){
     for(int jj = 0; jj < elasptr->ndime; jj++){

       double first  = (mu*exp(gamma*(I1-2)) + alpha*I1 + alpha)*delta(ii,jj);
       double second = -0.5*alpha*(elasptr->rcgreen[ii][jj] + elasptr->rcgreen[jj][ii]);
       double third  =  0.5*beta*(I3-1)*(rcgreeninvtpose[ii][jj] + rcgreeninvtpose[jj][ii]);
       
       elasptr->secpk[ii][jj] = 2*(first + second + third);
     }
    }

    // calculate the cauchy stress

    

    // calculate the first PK
    
    // calculate the the spatial tangent modulus
    

    
    // Gradient with respect to material parameters
    // check if required to calculate the gradient
   if (elasptr->grdflag == 1){

     // loop over pk2grad matrix entries
     for(int irow = 0;irow < elasptr->ndime; irow++){
       for(int icol = 0;icol < elasptr->ndime;icol++){

	 if(elasptr->grdparm == 0){
	   // taking the gradient with first stiffness parameter: \mu
	   double first   = (exp(gamma*(I1-2)))*delta(irow,icol);
	   double second  = I1*delta(irow,icol) + delta(irow,icol);
	   double third   = elasptr->rcgreen[irow][icol] + elasptr->rcgreen[icol][irow];
	   elasptr->pk2grad[irow][icol] = 2*(first - 0.5*second + 0.25*(third))
;
	   //elasptr->pk2grad[irow][icol] = 0;
	 }
	 
	 if(elasptr->grdparm == 1){
	   // taking the gradient with the second stiffness parameter \gamma
	   double first = mu*(I1-2)*exp(gamma*(I1-2))*delta(irow,icol);
	   //elasptr->pk2grad[irow][icol] = 0;
	   elasptr->pk2grad[irow][icol] = 2*(first);
	 }

	 if(elasptr->grdparm == 2){
	    double third  =  (I3-1)*(rcgreeninvtpose[irow][icol] + rcgreeninvtpose[icol][irow]);
       
	   elasptr->pk2grad[irow][icol] = third;
	 }

       }//icol
     }//irow
     
   }//gradflg
    
  }

  if (lawnumber == 5 ){

    
    double beta   = 100;
    double mu     = 1;
    double alpha  = -mu;

    double gamma  = elasptr->stfparam[0];
    
    double I1   = tracemat(elasptr->ndime,elasptr->rcgreen);
    double I3   = determinant(elasptr->ndime, elasptr->rcgreen);

    // CIJKL : same as Veronda Westman 3D
    for(int ii = 0;ii<elasptr->ndime;ii++){
      for(int jj = 0;jj<elasptr->ndime;jj++){
	for(int kk = 0;kk<elasptr->ndime;kk++){
	  for(int ll = 0;ll<elasptr->ndime;ll++){
	    
	    double first = mu*gamma*exp(gamma*(I1-2))*delta(kk,ll)*delta(ii,jj);
	    double second = alpha*delta(kk,ll)*delta(ii,jj) - 0.5*alpha*(delta(ii,kk)*delta(jj,ll) + delta(jj,kk)*delta(ii,ll));

	    double symm1  = 2*rcgreeninvtpose[kk][ll]*rcgreeninvtpose[ii][jj];
	    double symm2  = rcgreeninvtpose[ii][kk]*rcgreeninvtpose[jj][ll] + rcgreeninvtpose[jj][kk]*rcgreeninvtpose[ii][ll];

	    double third  = 0.5*beta*(I3*symm1 - (I3-1)*symm2);

	   	    
	    elasptr->cijklptr[ii][jj][kk][ll] = 4*(first + second + third);

	  }//ll
	}//kk
      }//jj
    } //ii
    
    // Second PK
    for( int ii = 0; ii < elasptr->ndime; ii ++ ){
     for(int jj = 0; jj < elasptr->ndime; jj++){

       double first  = (mu*exp(gamma*(I1-2)) + alpha*I1 + alpha)*delta(ii,jj);
       double second = -0.5*alpha*(elasptr->rcgreen[ii][jj] + elasptr->rcgreen[jj][ii]);
       double third  =  0.5*beta*(I3-1)*(rcgreeninvtpose[ii][jj] + rcgreeninvtpose[jj][ii]);
       
       elasptr->secpk[ii][jj] = 2*(first + second + third);
     }
    }
    
    // Gradient with respect to material parameters
    // check if required to calculate the gradient
   if (elasptr->grdflag == 1){

     // loop over pk2grad matrix entries
     for(int irow = 0;irow < elasptr->ndime; irow++){
       for(int icol = 0;icol < elasptr->ndime;icol++){

	 
	 if(elasptr->grdparm == 0){
	   // taking the gradient with the second stiffness parameter \gamma
	   double first = mu*(I1-2)*exp(gamma*(I1-2))*delta(irow,icol);
	   //elasptr->pk2grad[irow][icol] = 0;
	   elasptr->pk2grad[irow][icol] = 2*(first);
	 }

       }//icol
     }//irow
     
   }//gradflg
    
  }


  // this law (lawnumber == 6) ignores \gamma. Justifiable in the limit of small deformations
  
  if (lawnumber == 6 ){
    
    
    double beta   = 100;
    double mu     = elasptr->stfparam[0];
    double alpha  = -mu/2.0;
    
    double gamma  = 0;
    
    double I1   = tracemat(elasptr->ndime,elasptr->rcgreen);
    double I3   = determinant(elasptr->ndime, elasptr->rcgreen);

    // CIJKL : same as Veronda Westman 3D
    for(int ii = 0;ii<elasptr->ndime;ii++){
      for(int jj = 0;jj<elasptr->ndime;jj++){
	for(int kk = 0;kk<elasptr->ndime;kk++){
	  for(int ll = 0;ll<elasptr->ndime;ll++){
	    
	    double first = mu*gamma*exp(gamma*(I1-2))*delta(kk,ll)*delta(ii,jj);
	    double second = alpha*delta(kk,ll)*delta(ii,jj) - 0.5*alpha*(delta(ii,kk)*delta(jj,ll) + delta(jj,kk)*delta(ii,ll));

	    double symm1  = 2*rcgreeninvtpose[kk][ll]*rcgreeninvtpose[ii][jj];
	    double symm2  = rcgreeninvtpose[ii][kk]*rcgreeninvtpose[jj][ll] + rcgreeninvtpose[jj][kk]*rcgreeninvtpose[ii][ll];

	    double third  = 0.5*beta*(I3*symm1 - (I3-1)*symm2);

	   	    
	    elasptr->cijklptr[ii][jj][kk][ll] = 4*(first + second + third);

	  }//ll
	}//kk
      }//jj
    } //ii
    
    // Second PK
    for( int ii = 0; ii < elasptr->ndime; ii ++ ){
     for(int jj = 0; jj < elasptr->ndime; jj++){

       double first  = (mu*exp(gamma*(I1-2)) + alpha*I1 + alpha)*delta(ii,jj);
       double second = -0.5*alpha*(elasptr->rcgreen[ii][jj] + elasptr->rcgreen[jj][ii]);
       double third  =  0.5*beta*(I3-1)*(rcgreeninvtpose[ii][jj] + rcgreeninvtpose[jj][ii]);
       
       elasptr->secpk[ii][jj] = 2*(first + second + third);
     }
    }
    
    // Gradient with respect to material parameters
    // check if required to calculate the gradient
   if (elasptr->grdflag == 1){

     // loop over pk2grad matrix entries
     for(int irow = 0;irow < elasptr->ndime; irow++){
       for(int icol = 0;icol < elasptr->ndime;icol++){

	 
	 if(elasptr->grdparm == 0){
	   // taking gradient with respect to $\mu$ only $\alpha=0.5\mu$
	   double first   = (exp(gamma*(I1-2)))*delta(irow,icol);
	   double second  = I1*delta(irow,icol) + delta(irow,icol);
	   double third   = elasptr->rcgreen[irow][icol] + elasptr->rcgreen[icol][irow];
	   elasptr->pk2grad[irow][icol] = 2*(first - 0.5*second + 0.25*(third));
	     
	 }
	 
       }//icol
     }//irow
     
   }//gradflg
    
  }

  //  only gradient with respect to gamma
  //  gradient with respect to mu is set to zero.
  if (lawnumber == 7 ){
    
    
    double beta   = 100;
    double mu     = elasptr->stfparam[0];
    double alpha  = -mu/2.0;
     
    double gamma  = elasptr->stfparam[1];
    
    double I1   = tracemat(elasptr->ndime,elasptr->rcgreen);
    double I3   = determinant(elasptr->ndime, elasptr->rcgreen);

    // CIJKL : same as Veronda Westman 3D
    for(int ii = 0;ii<elasptr->ndime;ii++){
      for(int jj = 0;jj<elasptr->ndime;jj++){
	for(int kk = 0;kk<elasptr->ndime;kk++){
	  for(int ll = 0;ll<elasptr->ndime;ll++){
	    
	    double first = mu*gamma*exp(gamma*(I1-2))*delta(kk,ll)*delta(ii,jj);
	    double second = alpha*delta(kk,ll)*delta(ii,jj) - 0.5*alpha*(delta(ii,kk)*delta(jj,ll) + delta(jj,kk)*delta(ii,ll));

	    double symm1  = 2*rcgreeninvtpose[kk][ll]*rcgreeninvtpose[ii][jj];
	    double symm2  = rcgreeninvtpose[ii][kk]*rcgreeninvtpose[jj][ll] + rcgreeninvtpose[jj][kk]*rcgreeninvtpose[ii][ll];

	    double third  = 0.5*beta*(I3*symm1 - (I3-1)*symm2);

	   	    
	    elasptr->cijklptr[ii][jj][kk][ll] = 4*(first + second + third);

	  }//ll
	}//kk
      }//jj
    } //ii
    
    // Second PK
    for( int ii = 0; ii < elasptr->ndime; ii ++ ){
     for(int jj = 0; jj < elasptr->ndime; jj++){

       double first  = (mu*exp(gamma*(I1-2)) + alpha*I1 + alpha)*delta(ii,jj);
       double second = -0.5*alpha*(elasptr->rcgreen[ii][jj] + elasptr->rcgreen[jj][ii]);
       double third  =  0.5*beta*(I3-1)*(rcgreeninvtpose[ii][jj] + rcgreeninvtpose[jj][ii]);
       
       elasptr->secpk[ii][jj] = 2*(first + second + third);
     }
    }
    
    // Gradient with respect to material parameters
    // check if required to calculate the gradient
   if (elasptr->grdflag == 1){

     // loop over pk2grad matrix entries
     for(int irow = 0;irow < elasptr->ndime; irow++){
       for(int icol = 0;icol < elasptr->ndime;icol++){

	 
	 if(elasptr->grdparm == 0){
	   // setting \mu gradient to zero
	   elasptr->pk2grad[irow][icol] = 0;
  
	 }
	 
	 if(elasptr->grdparm == 1){
	   // taking the gradient with the second stiffness parameter \gamma
	   double first = mu*(I1-2)*exp(gamma*(I1-2))*delta(irow,icol);
	   //elasptr->pk2grad[irow][icol] = 0;
	   elasptr->pk2grad[irow][icol] = 2*(first);
	 }

	 
       }//icol
     }//irow
     
   }//gradflg
  
  } // closes lawnumber 7


  if(lawnumber == 8){
    double mu = elasptr->stfparam[1];
    double lambda = elasptr->stfparam[0];

    double dett   = determinant(elasptr->ndime, elasptr->rcgreen);

    for(int ii = 0;ii<elasptr->ndime;ii++){
      for(int jj = 0;jj<elasptr->ndime;jj++){
	for(int kk = 0;kk<elasptr->ndime;kk++){
	  for(int ll = 0;ll<elasptr->ndime;ll++){
	    double symmfourth1  = rcgreeninvtpose[ii][kk]*rcgreeninvtpose[jj][ll];
	    double symmfourth2  = rcgreeninvtpose[jj][kk]*rcgreeninvtpose[ii][ll];
	    double symmfourth   = 0.5*(symmfourth1+symmfourth2);
	    double Cinv_ijCinv_kl  = rcgreeninvtpose[ii][jj]*rcgreeninvtpose[kk][ll];

	    double first  = (lambda*0.5*Cinv_ijCinv_kl);
	    double second = (-mu + 0.5*lambda*log(dett))*symmfourth;
	    elasptr->cijklptr[ii][jj][kk][ll] = 2*(first -second);

	  }
	}
      }
    }

    for(int irow = 0; irow<elasptr->ndime;irow++){
	for(int icol = 0;icol<elasptr->ndime;icol++){
	  
	  
	  elasptr->secpk[irow][icol] = mu*delta(irow,icol) + (-mu + lambda*0.5*log(dett))*rcgreeninvtpose[irow][icol];
	  
	}
    }

  }// (lawnumber ==8 from IJNME Piltner Taylor 99)



   if(lawnumber == 9){
    
    // This is from Aremro Glaser: Simo Hughes, different penalty function
    // calculate C_IJKL
    double lambda = elasptr->stfparam[0];
    double mu     = elasptr->stfparam[1];
    double dett   = determinant(elasptr->ndime, elasptr->rcgreen);

    for(int ii = 0;ii<elasptr->ndime;ii++){
      for(int jj = 0;jj<elasptr->ndime;jj++){
	for(int kk = 0;kk<elasptr->ndime;kk++){
	  for(int ll = 0;ll<elasptr->ndime;ll++){
	    
	    double symmfourth1  = rcgreeninvtpose[ii][kk]*rcgreeninvtpose[jj][ll];
	    double symmfourth2  = rcgreeninvtpose[jj][kk]*rcgreeninvtpose[ii][ll];
	    double symmfourth   = 0.5*(symmfourth1+symmfourth2);
	    double Cinv_ijCinv_kl  = rcgreeninvtpose[ii][jj]*rcgreeninvtpose[kk][ll];
	    
	      
	    double first  = (2*mu - lambda*log(dett))*symmfourth;
	    double second = lambda*Cinv_ijCinv_kl;
	    
	    elasptr->cijklptr[ii][jj][kk][ll] = (first + second);
	  }
	}
      }
    }
    
     
    // calculate second piola kirchhoff
     for(int irow = 0; irow<elasptr->ndime;irow++){
	for(int icol = 0;icol<elasptr->ndime;icol++){
	  
	  double right  = rcgreeninvtpose[irow][icol];
	  double secpiola  = lambda*0.5*log(dett)*right + mu*(delta(irow,icol) -right);
	  
	  elasptr->secpk[irow][icol] = secpiola;
	}
      }


    // calculate the gradient contribution
    // check if we are required to calculate it
     if(elasptr->grdflag == 1){
       // the purpose is to set pk2grad

       if(elasptr->grdparm == 0){
	 // taking gradient with first stiffness parameter, \lambda
	 
	 for(int ii = 0;ii<elasptr->ndime;ii++){
	   for(int jj =0;jj<elasptr->ndime;jj++){
	     // note: dett is the determinant of C. det(C)= det(J^2)
	     // elasptr->pk2grad[ii][jj] = 0.5*(dett-1.0)*rcgreeninvtpose[ii][jj];
	     elasptr->pk2grad[ii][jj] = 0.5*log(dett)*rcgreeninvtpose[ii][jj];
	   }
	 }
       }

       if(elasptr->grdparm == 1){
	 // taking gradient with second stiffness paramter, \mu
	 for(int ii = 0; ii<elasptr->ndime;ii++){
	   for(int jj =0; jj<elasptr->ndime;jj++){
	     elasptr->pk2grad[ii][jj] = (delta(ii,jj)-rcgreeninvtpose[ii][jj]);
	   }
	 }
       }
     }
    
     
   }

   if (lawnumber == 10){
     double lambda  = elasptr->stfparam[0];
     double mu      = elasptr->stfparam[1];
     double alpha   = elasptr->stfparam[2];

     double I1   = tracemat(elasptr->ndime,elasptr->rcgreen);
     double dett    = determinant(elasptr->ndime,elasptr->rcgreen);

      for(int ii = 0;ii<elasptr->ndime;ii++){
	for(int jj = 0;jj<elasptr->ndime;jj++){
	  for(int kk = 0;kk<elasptr->ndime;kk++){
	    for(int ll = 0;ll<elasptr->ndime;ll++){
	      double symmfourth1  = rcgreeninvtpose[ii][kk]*rcgreeninvtpose[jj][ll];
	      double symmfourth2  = rcgreeninvtpose[jj][kk]*rcgreeninvtpose[ii][ll];
	      double symmfourth   = 0.5*(symmfourth1+symmfourth2);
	      double Cinv_ijCinv_kl  = rcgreeninvtpose[ii][jj]*rcgreeninvtpose[kk][ll];
	      
	      double first = 4*(lambda*0.25*Cinv_ijCinv_kl + 0.5*mu*(alpha-1.0)*pow(I1,alpha-2.0)*delta(ii,jj)*delta(kk,ll));
	      double second = 4*(mu*0.5 - (lambda/4.0)*log(dett))*symmfourth;
			      
	      elasptr->cijklptr[ii][jj][kk][ll] = first + second;
	      
	    }
	  }
	}
      }

       for(int irow = 0; irow<elasptr->ndime;irow++){
	 for(int icol = 0;icol<elasptr->ndime;icol++){
	   elasptr->secpk[irow][icol] = 2*((0.25*lambda*log(dett) - 0.5*mu)*rcgreeninvtpose[irow][icol] + 0.5*mu*pow(I1,alpha-1.0)*delta(irow,icol));
	}
       }
      

   }// lawnumber = 10

   
   if(lawnumber == 11){
     // plane stress, incompressible veronda westman
     double mu = elasptr->stfparam[0];
     double gamma = elasptr->stfparam[1];
     
     double I1 = tracemat(elasptr->ndime,elasptr->rcgreen);
     double I2 = determinant(elasptr->ndime,elasptr->rcgreen);
     double I11I3 = (I1 + (1.0/I2) -3);
     double GG = gamma*I11I3;


     for(int ii = 0; ii < elasptr->ndime; ii++){
       for(int jj = 0; jj < elasptr->ndime; jj++){
	 for(int kk = 0; kk < elasptr->ndime; kk++){
	   for(int ll = 0; ll < elasptr->ndime; ll++){
	   
	     double D1=    mu*(2*gamma*exp(GG)*(delta(kk,ll)-(rcgreeninvtpose[kk][ll]/I2)) + (rcgreeninvtpose[kk][ll]/I2))*delta(ii,jj);
	     double D2_1 = mu*((delta(kk,ll)/I2) - (I1*rcgreeninvtpose[kk][ll]/I2))*rcgreeninvtpose[ii][jj];
	     double D2_2 = mu*(-I2*rcgreeninvtpose[kk][ll] + (2*exp(GG)*rcgreeninvtpose[kk][ll]/I2))*rcgreeninvtpose[ii][jj];;
	     double D2_3 = mu*(-2*gamma*exp(GG)/I2)*(delta(kk,ll)-(rcgreeninvtpose[kk][ll]/I2))*rcgreeninvtpose[ii][jj];
	     double B2= -mu*((I1/I2) -I2 - (2*exp(GG)/I2));
	     double symmfourth=0.5*(rcgreeninvtpose[ii][kk]*rcgreeninvtpose[jj][ll] + rcgreeninvtpose[jj][kk]*rcgreeninvtpose[ii][ll]);

	     elasptr->cijklptr[ii][jj][kk][ll] = 2*(D1 + D2_1 + D2_2 + D2_3 + B2*symmfourth);
	     
	   }
	 }
       }
     }

     for( int ii = 0; ii < elasptr->ndime; ii ++ ){
       for(int jj = 0; jj < elasptr->ndime; jj++){
	 elasptr->secpk[ii][jj] = mu*(2*exp(GG)-(1.0/I2))*delta(ii,jj) + mu*((I1/I2) -I2 - (2*exp(GG)/I2))*rcgreeninvtpose[ii][jj];
       }
     }


     // now the gradient bit
     // check if we are required to calculate the gradient 
     if(elasptr->grdflag == 1){
       if(elasptr->grdparm == 0){

	 //taking gradient with respect to \mu
	 for(int ii = 0;ii<elasptr->ndime;ii++){
	   for(int jj =0;jj<elasptr->ndime;jj++){
	     elasptr->pk2grad[ii][jj] = (2*exp(GG)-(1/I2))*delta(ii,jj) + ((I1/I2) - (I2)-(2*exp(GG)/I2))*rcgreeninvtpose[ii][jj];
	   }
	 }
       }// grfparm == 0

       if(elasptr->grdparm == 1){
	 for(int ii = 0;ii<elasptr->ndime;ii++){
	   for(int jj =0;jj<elasptr->ndime;jj++){
	     elasptr->pk2grad[ii][jj] = mu*(1.0*delta(ii,jj) - (1.0/I2)*rcgreeninvtpose[ii][jj])*(2*exp(GG))*I11I3;
	   }
	 }
       } // elasptr->grdparm == 1

     }// grdflag

   }  // lawnumber == 11; incompressible plane stress vw
 


   // begin lawnumber == 12; incompressible plane stress neo-hook
   if (lawnumber == 12)
   {
     double mu = elasptr->stfparam[0];

     double I1 = tracemat(elasptr->ndime,elasptr->rcgreen);
     double I2 = determinant(elasptr->ndime,elasptr->rcgreen);
     
     // calculate second pk
     for(int ii = 0; ii < elasptr->ndime; ii++){
       for(int jj = 0; jj < elasptr->ndime; jj++){
	 elasptr->secpk[ii][jj] = mu*(delta(ii,jj) - (rcgreeninvtpose[ii][jj]/I2));
       }
     }

     // calculate cijkl

     for(int ii = 0; ii < elasptr->ndime; ii++){
       for(int jj = 0; jj < elasptr->ndime; jj++){
	 for(int kk = 0; kk < elasptr->ndime; kk++){
	   for(int ll = 0; ll < elasptr->ndime; ll++){
	     
	     double symmfourth1  = rcgreeninvtpose[ii][kk]*rcgreeninvtpose[jj][ll];
	     double symmfourth2  = rcgreeninvtpose[jj][kk]*rcgreeninvtpose[ii][ll];
	     double symmfourth   = 0.5*(symmfourth1+symmfourth2);
	     double Cinv_ijCinv_kl  = rcgreeninvtpose[ii][jj]*rcgreeninvtpose[kk][ll];
	     
	     elasptr->cijklptr[ii][jj][kk][ll] = 2*mu*((Cinv_ijCinv_kl/I2) + (1.0/I2)*(symmfourth));
	   }
	 }
       }
     }

     // calculate the gradient 
     if(elasptr->grdflag==1){
       assert(elasptr->grdparm == 0);  // only one parameter to take the grad

        for(int ii = 0;ii<elasptr->ndime;ii++){
	   for(int jj =0;jj<elasptr->ndime;jj++){
	     // note: dett is the determinant of C. det(C)= det(J^2)
	     elasptr->pk2grad[ii][jj] = (delta(ii,jj)-(1.0*rcgreeninvtpose[ii][jj]/I2));
	   }
	 }
	 
     }
     
     
   } // lawnumber == 12 incomp neo hook
  


  // in here we evaluate the first piola kirchhoff stress tensor
  // then the modulus Tensor A and its matrix representation
  // this is common to all constitutive laws

  // first the first piola-kirchhoff stress tensor
  // this can be written as P_{ij} = F_{ik}S_{kj}

  zeromat(elasptr->ndime,elasptr->ndime,elasptr->firstpk);

  for(int ii = 0 ; ii < elasptr->ndime; ii++){
    for(int jj = 0; jj < elasptr->ndime; jj++){
      for(int kk = 0; kk < elasptr->ndime; kk++){
	elasptr->firstpk[ii][jj] += elasptr->defgrad[ii][kk]*elasptr->secpk[kk][jj];
      }
    }
  }



  // the the modulus tensor A
  // A_{ijkl} = \pdd{P_{ij}}}{F_{kl}} = \pdd{F_{im}S_{mj}}{F_{kl}}
  //          

  for(int ii = 0; ii < elasptr->ndime; ii++){
    for(int jj = 0; jj < elasptr->ndime; jj++){
      for(int kk = 0; kk < elasptr->ndime; kk++){
	for(int ll = 0; ll < elasptr->ndime; ll++){

	  
	  
	  elasptr->aijklptr[ii][jj][kk][ll] = delta(ii,kk)*elasptr->secpk[jj][ll];
	  
	  for(int mm = 0 ; mm < elasptr->ndime; mm++){
	    for(int qq = 0; qq < elasptr->ndime; qq++){
	      elasptr->aijklptr[ii][jj][kk][ll]  += (elasptr->defgrad[ii][mm]*elasptr->cijklptr[qq][ll][jj][mm]*elasptr->defgrad[kk][qq]);
	      
	    }
	  }



	}
      }
    }
  }//ii --close aijklptr

  
  
  
  double **FFF = elasptr->defgrad;

  double **FFINV;
  FFINV = allocmat(elasptr->ndime,elasptr->ndime,FFINV);
  invertsmall(elasptr->ndime, FFF,FFINV);
  
  for(int ii = 0; ii < elasptr->ndime; ii++){
    for(int jj = 0; jj < elasptr->ndime; jj++){
      for(int kk = 0; kk < elasptr->ndime; kk++){
	for(int ll = 0; ll < elasptr->ndime; ll++){
	  elasptr->sptan[ii][jj][kk][ll] = 0;
	  
	  for(int AA = 0; AA < elasptr->ndime; AA++){
	    for(int BB = 0; BB < elasptr->ndime; BB++){
	      for(int CC = 0; CC < elasptr->ndime; CC++){
		for(int DD = 0; DD < elasptr->ndime; DD++){ 
		  elasptr->sptan[ii][jj][kk][ll] += FFF[ii][AA]*FFF[jj][BB]*FFF[kk][CC]*FFF[ll][DD]*elasptr->cijklptr[AA][BB][CC][DD];
		}
	      }
	    }
	  }
	}
      }
    }
  }
	

  for(int ii = 0; ii < 2; ii++){
    for(int jj = 0; jj < 2; jj++){

      elasptr->krchstrs[ii][jj] = 0;

      for(int mm = 0; mm < 2; mm++){
	for(int kk = 0; kk<2; kk++){
	  elasptr->krchstrs[ii][jj] += FFF[ii][mm]*elasptr->secpk[mm][kk]*FFF[jj][kk];
	}
      }

    }
  }

  double newaa[2][2][2][2];
  for(int aa = 0; aa < 2 ; aa++){
    for(int BB = 0; BB < 2; BB++){
      for(int cc = 0; cc < 2; cc++){
	for(int DD = 0; DD < 2; DD++){
	  newaa[aa][BB][cc][DD] = 0;

	  for(int bb = 0; bb < 2 ; bb++){
	    for(int dd = 0; dd < 2; dd++){
	      newaa[aa][BB][cc][DD] +=FFINV[BB][bb]*(elasptr->sptan[aa][bb][cc][dd] + elasptr->krchstrs[aa][cc]*delta(bb,dd))*FFINV[DD][dd];
	    }
	  }
	 
	}
      }
    }
  }


  
  double dett  = determinant(elasptr->ndime,elasptr->rcgreen);
  double A1111 = 1.0 + (1.0/dett)*(1 + 1.5 - 1.5*0.5*log(dett))*elasptr->defgrad[1][1]*elasptr->defgrad[1][1];
  
  elasptr->AIJ[0][0] = elasptr->aijklptr[0][0][0][0];
  elasptr->AIJ[0][1] = elasptr->aijklptr[0][0][1][1];
  elasptr->AIJ[0][2] = elasptr->aijklptr[0][0][0][1];
  elasptr->AIJ[0][3] = elasptr->aijklptr[0][0][1][0];

  elasptr->AIJ[1][0] = elasptr->AIJ[0][1];
  elasptr->AIJ[1][1] = elasptr->aijklptr[1][1][1][1];
  elasptr->AIJ[1][2] = elasptr->aijklptr[1][1][0][1];
  elasptr->AIJ[1][3] = elasptr->aijklptr[1][1][1][0];

  elasptr->AIJ[2][0] = elasptr->AIJ[0][2];
  elasptr->AIJ[2][1] = elasptr->AIJ[1][2];
  elasptr->AIJ[2][2] = elasptr->aijklptr[0][1][0][1];
  elasptr->AIJ[2][3] = elasptr->aijklptr[0][1][1][0];

  elasptr->AIJ[3][0] = elasptr->AIJ[0][3];
  elasptr->AIJ[3][1] = elasptr->AIJ[1][3];
  elasptr->AIJ[3][2] = elasptr->AIJ[2][3];
  elasptr->AIJ[3][3] = elasptr->aijklptr[1][0][1][0];


  // do the same for CIJ
  elasptr->CIJ[0][0] = elasptr->cijklptr[0][0][0][0];
  elasptr->CIJ[0][1] = elasptr->cijklptr[0][0][1][1];
  elasptr->CIJ[0][2] = elasptr->cijklptr[0][0][0][1];
  elasptr->CIJ[0][3] = elasptr->cijklptr[0][0][1][0];

  elasptr->CIJ[1][0] = elasptr->CIJ[0][1];
  elasptr->CIJ[1][1] = elasptr->cijklptr[1][1][1][1];
  elasptr->CIJ[1][2] = elasptr->cijklptr[1][1][0][1];
  elasptr->CIJ[1][3] = elasptr->cijklptr[1][1][1][0];

  elasptr->CIJ[2][0] = elasptr->CIJ[0][2];
  elasptr->CIJ[2][1] = elasptr->CIJ[1][2];
  elasptr->CIJ[2][2] = elasptr->cijklptr[0][1][0][1];
  elasptr->CIJ[2][3] = elasptr->cijklptr[0][1][1][0];

  elasptr->CIJ[3][0] = elasptr->CIJ[0][3];
  elasptr->CIJ[3][1] = elasptr->CIJ[1][3];
  elasptr->CIJ[3][2] = elasptr->CIJ[2][3];
  elasptr->CIJ[3][3] = elasptr->cijklptr[1][0][1][0];
  


  // decallocate the inverse transpose
  deallocmat(elasptr->ndime,rcgreeninvtpose);
  deallocmat(elasptr->ndime,FFINV);
  return 0;
}

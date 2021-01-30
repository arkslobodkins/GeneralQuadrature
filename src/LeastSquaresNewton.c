#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../headers/LeastSquaresNewton.h"
#include "../headers/GetFunction.h"
#include "../headers/GetJacobian.h"
#include "../headers/StoreCopy.h"
#include "../headers/Structs.h"
typedef int boolean;
boolean checkInfAndNan(int k, double *x);  
static void FreeMemory(double *, double *, double *, double *, double *, double *, double *, double *);  

/* LeastSquaresNewton
 * Receives nodes and weights as an initial guess for 
 * specified quadrature parameters. Primarily  solves 
 * underdetermined systems of equations in the least
 * squares sense. Returns success if algorithm converged 
 * and all nodes are inside of the domain.
 *   
 * x: initial nodes
 * w: initial weights
 * k: number of nodes and weights
 * flag: returns 1 if Newton's method succeeded, 0 otherwise 
 * its: number of iterations for Newton's method 
 * functions: structure that contains function handles for the respective domain
 * params: structure that contains quadrature parameters
 */
void LeastSquaresNewton( double *x, double *w,int k, int *flag, int *its, struct functionHandles *functions, struct quadratureParameters *params) {
  int i,j,d;
  int iterations=0; int maxiter=25;          
  int Iteration;
  int p=params->degreeOfPrecision;
  int m=params->numberOfBasisFunctions;
  int dim=params->totalDimension;
  int numRows=(dim+1)*k, numColumns=m;  
  double difference, sum, normError, normErrorPrev;
  double alpha=1.0, tol=pow(10.0,-15);            //damping parameter, convergence criterion
  double *jacobian=malloc(numRows*m*sizeof(double));
  double *xPrev=malloc(numRows*sizeof(double));
  double *xNext=malloc(numRows*sizeof(double));
  double *xDiff=malloc(numRows*sizeof(double));
  double *f=malloc(m*sizeof(double));
  double *leastSqSol=malloc(numRows*sizeof(double)); 
  //Lapack parameters
  int info,nRhs=1, lda=numRows, ldb=(numRows<numColumns)? numColumns: numRows, lwork=numRows+numColumns,\
  nmin= (numRows<numColumns)? numRows: numColumns;
  double *fStore=malloc(ldb*sizeof(double)); 
  double *work=malloc(lwork*sizeof(double)); 
  char trans='T';

  for(j=0;j<k;j++) {
    for(d=0; d<=dim;d++) {
      xPrev[(dim+1)*j+d]=0.0; 
    }
  xNext[(dim+1)*j]=w[j];
    for(d=0; d<dim;d++) {
      xNext[(dim+1)*j+d+1]=x[dim*j+d];  
    }   
  }
  for(j=0;j<numRows;j++) {
    xDiff[j]=xNext[j]-xPrev[j];
  }
  normError=ScaledTwoNorm(numRows, xDiff);
  normErrorPrev=normError;
  while(iterations<maxiter&& normError>tol){
    sum=0.0;
    //function and jacobian of a function to be solved
    CopyNodesAndWeights(k, dim, xPrev, xNext);
    GetJacobian(k, x,w,jacobian, functions, params);  
    GetFunction(k, &x[0],&w[0],&f[0], functions, params);
    //terminate the algorithm if infs of nans entountered
    if((checkInfAndNan( numRows*m, jacobian)==1) || (checkInfAndNan( m, f)==1)) {
      *flag=0;
      FreeMemory( jacobian, xPrev, xNext, xDiff, f, leastSqSol, work, fStore);
      return;
    }
    for(i=0; i<nmin;i++) {
      fStore[i]=f[i];
    }
    for(i=nmin; i<ldb;i++) {
      fStore[i]=0.0;  
    }
    //solve the system using least squares
    dgels_( &trans, &numRows, &numColumns, &nRhs, jacobian, &lda, fStore, &ldb, work, &lwork, &info);
    for(i=0;i<numRows;i++) {
      leastSqSol[i]=fStore[i];
    }   
    if(checkInfAndNan(numRows, leastSqSol)==1) {
      *flag=0; 
      FreeMemory( jacobian, xPrev, xNext, xDiff, f, leastSqSol, work, fStore);
      return;
    }
    for(j=0;j<numRows;j++) {
      xNext[j]=xPrev[j]-alpha*leastSqSol[j];    
    }
    RollNodesAndWeights(k, dim, xNext, x, w);
  
    for(j=0; j<numRows;j++) {
      xDiff[j]=xNext[j]-xPrev[j];
    }
    normError=ScaledTwoNorm(numRows, xDiff);
    iterations++; 
    if(iterations>8&&functions->InDomain(k, x, params->dimension)==0) {
      *flag=0; 
      FreeMemory( jacobian, xPrev, xNext, xDiff, f, leastSqSol, work, fStore);
      return;
    }
    if(iterations>6 && normError> (normErrorPrev+0.5)) {
      *flag=0;
      FreeMemory( jacobian, xPrev, xNext, xDiff, f, leastSqSol, work, fStore);
      return;
    }
    normErrorPrev=normError;
    
  }//end while loop
  FreeMemory( jacobian, xPrev, xNext, xDiff, f, leastSqSol, work, fStore);
  *its=iterations;

  //check whether nodes are inside the domain
  if(functions->InDomain(k, x, params->dimension)==0) {
     *flag=0; 
     return;
  }
  //check whether Newton's method converged in specified number of iterations
  if(iterations<maxiter &&isnan(normError)==0 &&isinf(normError)==0) {    
    *flag=1;
    Iteration=i;
  } else if(isnan(normError)) { 
    *flag=0;
    Iteration=maxiter;
  
  } else if(iterations>=maxiter) {  
    *flag=0;
    Iteration=maxiter;
  } else {
    *flag=0;
  }

} /*end LeastSquaresNewton */

double ScaledTwoNorm(int n, double *x){
  int i;
  double norm=0.0;
  for(i=0;i<n;i++) {    
    norm=norm+x[i]*x[i];
  }
  norm=sqrt(norm)/sqrt(n);    
  return norm;
} /* end ScaledTwoNorm */

double TwoNorm(int n, double *x){
  int i;
  double norm=0.0;

  for(i=0;i<n;i++) {    
    norm=norm+x[i]*x[i];
  } 
  norm=sqrt(norm);
  return norm;
} /* end TwoNorm */

double InfNorm(int n, double *x) {
  int i;
  double norm=0.0;
  double temp;

  for(i=0;i<n;i++) {    
    temp=fabs(x[i]);
    if(temp>=norm) {      
      norm=temp;      
    }    
  }
  return norm;    
} /* end InfNorm */

boolean checkInfAndNan(int k, double *x) {
  int i;
  boolean flag;
  flag=0;
  for(i=0; i<k; i++) {
    if (isnan(x[i])==1 || isinf(x[i])==1) {
      flag=1;
      goto returnValue;
    } 
  }

  returnValue: return flag;
} /*end checkInfAndNan */


static void FreeMemory(double *jacobian, double *xPrev, double *xNext, double *xDiff, double *f, double *leastSqSol, double *work, double *fStore) {
  free(jacobian);
  free(xPrev);
  free(xNext);
  free(xDiff);
  free(f);
  free(leastSqSol);
  free(work);
  free(fStore);
}

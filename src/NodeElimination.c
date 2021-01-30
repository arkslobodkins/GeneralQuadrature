#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../headers/NodeElimination.h"
#include "../headers/GetJacobian.h"
#include "../headers/InsertionSort.h"
#include "../headers/LeastSquaresNewton.h"
#include "../headers/TestIntegral.h"
#include "../headers/StoreCopy.h"
#include "../headers/Structs.h"
#include "../headers/InDomain.h"

/* NodeElimination
 * N: initial number of nodes
 * xInitial: initial nodes
 * wInitial: initial weights 
 * nNodesNew: number of nodes of the new quadrature
 * xFinal: new quadrature nodes
 * wFinal new quadrature weights
 * functions: structure that contains function handles for the respective domain
 * params: structure that contains quadrature parameters
 * 
 * Receives quadrature nodes, weights, domain and quadrature parameters and 
 * uses the new node elimination scheme to eliminate one node to obtain initial guess for 
 * Newton's method. Subsequently, routine calls Newton's method to obtain quadrature rule 
 * with fewer nodes. The procedure is repeated until no more nodes can be eliminated.     
 */
void NodeElimination(int N, double *xInitial,double *wInitial,  int *nNodesNew, double *xFinal, double *wFinal, struct functionHandles *functions, struct quadratureParameters *params, struct eliminationHistory *history) {
  
  int i,j,k,d;  
  int its, flag, counter=0, countSuccessLocal=0;
  int dim=params->totalDimension;
  int m=params->numberOfBasisFunctions;
  history->totalEliminations=0;
  int *totalEliminations= &history->totalEliminations;
  double res, normq2Row, tol=pow(10,-15);
  double *xTemp, *wTemp, *xNew, *wNew, *signIndex;
  int *arrayIndex;
  
  //LAPACK parameters
  int numRows, numCols, qiCols,lda, ldc,lwork, info, lworkqr;
  numRows=(dim+1)*N,  numCols=m; 
  double *jacobian,  *Q, *dZ, *Z;
  char side='L'; char trans='T'; 
  int q2Cols=1;
  double *x, *tau, *work, *workqr, *q2Row;
  xNew=malloc( dim*N*sizeof(double));
  wNew=malloc( N*sizeof(double));
  xTemp=malloc( dim*N*sizeof(double));
  wTemp=malloc( N*sizeof(double));
  arrayIndex=malloc(N*sizeof(int));

  //initialize nodes and weights 
  CopyNodes( N, dim, xNew,xInitial);
  CopyWeights(N, wNew, wInitial);
  CopyNodes(N, dim,  xTemp, xInitial);
  CopyWeights(N, wTemp, wInitial);

  res=TestIntegral(N, xNew,wNew, functions, params);  
  printf("initial residual in NodeElimination=%lf \n", res);
  //run LeastSquaresNewton if initial error contains guess greater than tolerance
  if(fabs(res)>tol) {
    LeastSquaresNewton(&xTemp[0],&wTemp[0],N,&flag, &its, functions, params);
    
    if(flag==1) {
      CopyNodes(N, dim, xNew, xTemp);
      CopyWeights(N, wNew, wTemp);         
    } else {
    printf("Poor initial guess, exiting program");
    exit(-1);
    }
  }
  k=N-1;

  Z=malloc(1*sizeof(double));
  //run Node Elimination Algorithm. Theoretical optimum is reached when (dim+1)*k=m, but a few more iterations are
  //allowed in case under exceptional circumstances theoretical optimum is surpassed. 
  while( (dim+1)*k>m-3) {

    printf("dimension=%i, current number of nodes=%i,  optimal=%i \n", dim, k+1, (int)ceil(1.0*m/(dim+1.0)));
    numRows=(dim+1)*(k+1), numCols=m, lda=numRows, ldc=numRows, lworkqr=numCols, qiCols=k+1;
    //reallocate to consume less memory 
    xTemp=realloc( xTemp, dim*N*sizeof(double));
    wTemp=realloc( wTemp,  N*sizeof(double));
    xNew=realloc( xNew, dim*N*sizeof(double));
    wNew=realloc( wNew,  N*sizeof(double));   
    Z=realloc(Z, k*numRows*sizeof(double));
    jacobian=malloc(numRows*numCols*sizeof(double));      //Jacobian
    Q=malloc( numRows*(k+1)*sizeof(double));
    dZ=malloc( numRows*(k+1)*sizeof(double));
    signIndex=malloc( N*sizeof(double));
    x=malloc(numRows*sizeof(double));
    tau=malloc( numCols*sizeof(double));
    ldc=numCols<numRows? numRows:numCols;
    lworkqr=numCols<numRows? numRows: numCols;
    work=malloc(numRows*sizeof(double));
    workqr=malloc(lworkqr*sizeof(double));
    q2Row=malloc( numRows*sizeof(double));    
    arrayIndex=realloc(arrayIndex, (k+1)*sizeof(int));

    GetJacobian(k+1, xNew,wNew, jacobian, functions, params);
    //construct QR factorization
    dgeqr2_(&numRows, &numCols, jacobian, &lda, tau, work, &info);
    for( i=0; i<k+1;i++) {
      for(j=0; j<numRows;j++) {
        Q[i*numRows+j]=0.0;
      }
      Q[i*numRows+i*(dim+1)]=1.0;
    } 
  
    //obtain Q explicitly
    trans='T';
    dormqr_(&side, &trans, &numRows, &qiCols, &numCols, jacobian, &lda, tau, Q, &ldc, workqr, &lworkqr, &info);
      
    //compute initial guesses for Newton's method and store them in Z
    for(i=0; i<k+1; i++) {

      for(j=0; j<(numRows-m);j++) {
        q2Row[m+j]=Q[i*numRows+m+j];
      } 
      for(j=0; j<m; j++) {
        q2Row[j]=0.0;
      }  
    
      normq2Row=TwoNorm(numRows-m, &q2Row[m]);
      trans='N';
      dormqr_(&side, &trans, &numRows, &q2Cols, &numCols, jacobian, &lda, tau, q2Row, &ldc, workqr, &lworkqr, &info);
      for(j=0; j<numRows;j++) {
        dZ[i*numRows+j]=q2Row[j];     
        dZ[i*numRows+j]=dZ[i*numRows+j]*wNew[i]/(normq2Row*normq2Row);
      } 

      counter=0;
      for(j=0; j<k+1; j++) {
        if(i!=j) {  
          Z[k*(dim+1)*i+(dim+1)*counter]=wNew[j]-dZ[i*numRows+(dim+1)*j];
        }
        if(i!=j) {
          for( d=0; d<dim;d++) {
            Z[k*(dim+1)*i+(dim+1)*counter+1+d]=xNew[j*dim+d]-dZ[i*numRows+(dim+1)*j+d+1];
          }
          counter++;
        } 
      }
      signIndex[i]=TwoNorm(numRows, &dZ[numRows*i]);

    } 
  
    for(j=0; j<k+1;j++) { 
      arrayIndex[j]=j;
    }
    //store indices of signIndex in arrayIndex such that entries are accessed in ascending order 
    InsertionSort(k+1, signIndex ,arrayIndex);
    free(jacobian);
    free(Q);
    free(dZ);   
    free(signIndex);
    free(x);
    free(tau);
    free(work);
    free(workqr);
    free(q2Row);
 
    //extract initial guesses and run Newton's method
    for(i=0;i<k+1;i++) {

      for( j=0;j<k;j++) { 
        for( d=0; d<dim;d++) {          
          xTemp[dim*j+d]=Z[k*(dim+1)*arrayIndex[i]+(dim+1)*j+d+1];
        }
        wTemp[j]=Z[k*(dim+1)*arrayIndex[i]+(dim+1)*j];
      }

      LeastSquaresNewton(&xTemp[0],&wTemp[0],k,&flag, &its, functions, params);   
      
      //store nodes and weights if Newton's method succeeded, update history
      if(flag==1) {
        *nNodesNew=k;
        CopyNodes(k, dim, xNew, xTemp);
        CopyWeights(k, wNew, wTemp);
        history->nodesTotal[*totalEliminations]=k;  
        history->successNodeIndex[*totalEliminations]=i;     
        history->successNewtonIterations[*totalEliminations]=its;
        history->totalEliminations++;     
        break;  //break i loop    
      }  //end if     

    } //end i loop   */ 

    //end NodeElimination if all iterations were unsuccessful
    if(flag==0) {
      printf("all nodes for the last iteration did not converge, return quadrature \n \n");
      break;
    } 
    k--;

  }  //end while loop 

  //save nodes and weights in output array
  CopyNodes(*nNodesNew,dim, xFinal, xNew);
  CopyWeights(*nNodesNew, wFinal, wNew);
  free(Z);
  free(xTemp);
  free(wTemp);
  free(xNew);
  free(wNew);
  free(arrayIndex);
} /*end NodeElimination */

          

/* Arkadijs Slobodkins
 * SMU Mathematics
 * February 2021
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../headers/ComputeDomain.h"
#include "../headers/NodeElimination.h"
#include "../headers/GaussTensor.h"
#include "../headers/AddDimension.h"
#include "Gauss_Lib/Jacobi.h"
#include "../headers/BasisIndices.h"
#include "../headers/TestIntegral.h"
#include "../headers/Output.h"
#include "../headers/SetParams.h"
#include "../headers/SetDomain.h"
#include "../headers/StoreCopy.h"
#include "../headers/Structs.h"
static void EstimateMemory(struct quadratureParameters *);
static void FreeMemory(double *, double *, double *, double *, struct eliminationHistory *);
static void AllocateMemory(int dim, int nNodes,double **, double **, double **, double **, struct eliminationHistory **);
static void ReallocateMemory(int dim, int nNodes,double **, double **, double **, double **, struct eliminationHistory **);

/* ComputeCube
 * domainFunctions: structure that contains function handles for the cube
 * params: structure that contains quadrature parameters
 * 
 * Sets up the initial guess for the Node Elimination algorithm. 
 * Initial guess is computed via recursive scheme that reuses 
 * quadrature rules obtained in lower dimensions. More specifically, 
 * it computes the tensor product of (d-1)-dimensional quadrature 
 * rule over the unit cube and interval to obtain the initial guess for 
 * quadrature over the d-dimensional unit cube. When Node Elimination 
 * for the dimension assigned by a user is finalized, the final quadrature
 *  and analysis is printed to quadRule.txt and results.txt. 
*/
void ComputeCube(struct functionHandles *domainFunctions, struct quadratureParameters *params) {
  EstimateMemory(params);
  int dim=params->totalDimension;  
  int p=params->degreeOfPrecision;
  char flag[0];
  int n=floor(p/2)+1;
  int nNodesPrev, nNodesCur, nNodesNew;
  double res;
  double x[n], w[n];
  double *wInitial, *xInitial, *xNew, *wNew,*zNew;
  struct eliminationHistory *history;
  //generate Gaussian nodes and weights in 1-d on [0,1]
  double p1=0.0, p2=0.0;
  Jacobi(n, p1, p2, x,w);     
  int i,j,d;

  if(dim==1) {
    DumpCubatureRule(x, w, n, params, "CUBE");
    return;
  }

  //implement recursive scheme for computing the initial guess and run Node Elimination algorithm 
  for(d=2; d<=dim;d++) {
    SetParams(d, 0, 0,  p, params);

    if(d==2) {
      //set up parameters for 2-dim cube
      nNodesPrev=n;
      nNodesCur=n*n;
      nNodesNew=nNodesCur;	
      AllocateMemory(d, nNodesCur, &xInitial, &wInitial, &xNew, &wNew, &history);
      AddLine(d, n, n, x,x , w, w, xInitial, wInitial); 
		
    }else if(d>2) {
      //set up parameters for d-dim cube
      nNodesPrev=nNodesNew;
      nNodesCur=nNodesNew*n;	
      nNodesNew=nNodesCur;
      ReallocateMemory(d, nNodesCur, &xInitial, &wInitial, &xNew, &wNew, &history);
      AddLine(d, n, nNodesPrev, x,xNew, w, wNew, xInitial, wInitial); 
    }
    //perform Node Elimination
    NodeElimination(nNodesCur ,xInitial,wInitial, &nNodesNew,  xNew, wNew, domainFunctions, params, history);    
  }
  //print results to results.txt and quadRule.txt
  Output(xInitial, wInitial, nNodesCur,xNew, wNew, nNodesNew, domainFunctions, params,  history, "CUBE"); 
  DumpCubatureRule(xNew, wNew, nNodesNew, params, "CUBE");
  res=TestIntegral(nNodesNew, xNew,wNew, domainFunctions, params);
  printf("reached \n");
  printf("Final residual of the sum of all basis function integrals =%.16f \n", res); 
	
  FreeMemory(xInitial, wInitial, xNew, wNew, history);
} /* end ComputeCube */

/* ComputeSimplex
 * domainFunctions: structure that contains function handles for the simplex
 * params: structure that contains quadrature parameters
 * 
 * Sets up the initial guess for the Node Elimination algorithm. 
 * Initial guess is computed via recursive scheme that reuses 
 * quadrature rules obtained in lower dimensions. More specifically, 
 * it computes the tensor product of (d-1)-dimensional quadrature 
 * rule over the unit simplex and interval and uses variant of Duffy 
 * Transformation to obtain the initial guess for quadrature over the
 *  d-dimensional unit simplex. When Node Elimination for the dimension
 *  assigned by a user is finalized, the final quadrature and analysis 
 * is printed to quadRule.txt and results.txt. 
*/
void ComputeSimplex(struct functionHandles *domainFunctions, struct quadratureParameters *params) {
  EstimateMemory(params);
  char flag[0];
  int dim=params->totalDimension;  
  int p=params->degreeOfPrecision;
  int n;
  int nNodesPrev, nNodesCur, nNodesNew;
  double res;
  double *wInitial, *xInitial, *xNew, *wNew,*zNew;
  struct eliminationHistory *history;
  //input parameters for generating Gaussian rules
  double p1=0.0, p2=0.0;
  int i,j,d;

  for(d=2; d<=dim;d++) {
    SetParams(d, 0, 0,  p, params);

    if(d==2) {
      n=floor(p/2)+1+floor(p/8);
      double x[n], w[n];
      Jacobi(n, p1, p2, x,w);     //generate Gaussian nodes and weights in 1-d on [0,1]
      //set up parameters for 2-dim simplex
      nNodesPrev=n;
      nNodesCur=n*n;
      nNodesNew=nNodesCur;	
      AllocateMemory(d, nNodesCur, &xInitial, &wInitial, &xNew, &wNew, &history);
      AddLineSimplex(d, n, n, x,x , w, w, xInitial, wInitial); 
    } else if(d>2) {
      n=floor(p/2)+2+floor(p/6);
      double x[n], w[n];
      Jacobi(n, p1, p2, x,w);      //generate Gaussian nodes and weights in 1-d on [0,1]
      //set up parameters for d-dim simplex
      nNodesPrev=nNodesNew;
      nNodesCur=nNodesNew*n;	
      nNodesNew=nNodesCur;
      ReallocateMemory(d, nNodesCur, &xInitial, &wInitial, &xNew, &wNew, &history);
      AddLineSimplex(d, n, nNodesPrev, x,xNew, w, wNew, xInitial, wInitial); 
    }

    //perform Node Elimination
    NodeElimination(nNodesCur ,xInitial,wInitial, &nNodesNew,  xNew, wNew, domainFunctions, params, history);   
  }
  //print results to results.txt and quadRule.txt
  Output(xInitial, wInitial, nNodesCur,xNew, wNew, nNodesNew, domainFunctions, params,  history, "SIMPLEX"); 
  DumpCubatureRule(xNew, wNew, nNodesNew, params, "SIMPLEX");
  res=TestIntegral(nNodesNew, xNew,wNew, domainFunctions, params);
  printf("Final residual of the sum of all basis function integrals =%.16f \n", res); 
  
  FreeMemory(xInitial, wInitial, xNew, wNew, history);
} /* end ComputeSimplex */


/* ComputeCubeSimplex
 * domainFunctions: structure that contains function handles 
 * params: structure that contains quadrature parameters
 * 
 * Sets up the initial guess for the Node Elimination algorithm. 
 * Initial guess is computed via recursive scheme that reuses 
 * quadrature rules obtained in lower dimensions. At the first stage, 
 * it computes the tensor product of (d-1)-dimensional quadrature 
 * rule over the unit simplex and interval and uses variant of Duffy 
 * Transformation to obtain the initial guess for quadrature over the
 * d-dimensional unit simplex. Later on, it computes tensor product of 
 * (d-1)-dimensional CUBESIMPLEX and interval to obtain initial guess for quadrature 
 * over d-dimensional CUBESIMPLEX. When Node Elimination for the dimensions
 * assigned by a user is finalized, the final quadrature and analysis 
 * is printed to quadRule.txt and results.txt. 
*/
void ComputeCubeSimplex(struct functionHandles *domainFunctions, struct quadratureParameters *params) {
  EstimateMemory(params);
  int dimCubeMax=params->dimension[0];
  int dimSimplexMax=params->dimension[1];
  char flag[0];
  int dimCube, dimSimplex;
  int dim=params->totalDimension;  
  int p=params->degreeOfPrecision;
  int n;
  int nNodesPrev, nNodesCur, nNodesNew;
  double res;
  double *wInitial, *xInitial, *xNew, *wNew,*zNew;
  struct eliminationHistory *history;
  //input parameters for generating Gaussian rules
  double p1=0.0, p2=0.0;
  int i,j,d;

  //implement recursive scheme for computing the initial guess and run Node Elimination algorithm 
  for(d=2; d<=dim;d++) {
    dimSimplex=(d <= dimSimplexMax) ? (d) : (dimSimplexMax);   //when maximum dimension for simplex is reached, increase dimension for cube
    dimCube=d-dimSimplex;
	
    //set up parameters
    if(dimCube==0) {
      SetSimplex(domainFunctions);
      SetParams(dimSimplex, 0, 0, p, params);
      if(d==2) {
        n=floor(p/2)+1+floor(p/8);
        double x[n], w[n];
        Jacobi(n, p1, p2, x,w);     //generate Gaussian nodes and weights in 1-d on [0,1]
        nNodesPrev=n;
        nNodesCur=n*n;
        nNodesNew=nNodesCur;	
        AllocateMemory(d, nNodesCur, &xInitial, &wInitial, &xNew, &wNew, &history);     

        AddLineSimplex(d, n, n, x,x , w, w, xInitial, wInitial); 
      }else if(d>2) {
        n=floor(p/2)+2+floor(p/6);
        double x[n], w[n];
        Jacobi(n, p1, p2, x,w);           //generate Gaussian nodes and weights in 1-d on [0,1]
        nNodesPrev=nNodesNew;
        nNodesCur=nNodesNew*n;	
        nNodesNew=nNodesCur;
        ReallocateMemory(d, nNodesCur, &xInitial, &wInitial, &xNew, &wNew, &history);
        AddLineSimplex(d, n, nNodesPrev, x,xNew, w, wNew, xInitial, wInitial); 
      }
	
    } else if(dimCube>0) {
      SetCubeSimplex(domainFunctions);
      SetParams(dimCube, dimSimplex, 0, p, params);
      n=floor(p/2)+1;
      double x[n], w[n];
      Jacobi(n, p1, p2, x,w);           //generate Gaussian nodes and weights in 1-d on [0,1]
      nNodesPrev=nNodesNew;
      nNodesCur=nNodesNew*n;	
      nNodesNew=nNodesCur;
      ReallocateMemory(d, nNodesCur, &xInitial, &wInitial, &xNew, &wNew, &history);
      AddLine(d, n, nNodesPrev, x,xNew, w, wNew, xInitial, wInitial); 
    }
    //perform Node Elimination
    NodeElimination(nNodesCur ,xInitial,wInitial, &nNodesNew,  xNew, wNew, domainFunctions, params, history);    
  }
  //print results to results.txt and quadRule.txt
  Output(xInitial, wInitial, nNodesCur,xNew, wNew, nNodesNew, domainFunctions, params,  history, "CUBESIMPLEX"); 
  DumpCubatureRule(xNew, wNew, nNodesNew, params, "CUBESIMPLEX");
  res=TestIntegral(nNodesNew, xNew,wNew, domainFunctions, params);
  printf("Final residual of the sum of all basis function integrals =%.16f \n", res); 
  FreeMemory(xInitial, wInitial, xNew, wNew, history);
} /* end ComputeCubeSimplex */


/* ComputeSimplexSimplex
 * domainFunctions: structure that contains function handles 
 * params: structure that contains quadrature parameters
 * 
 * Sets up the initial guess for the Node Elimination algorithm. 
 * Initial guess is computed via recursive scheme that reuses 
 * quadrature rules obtained in lower dimensions. It computes 
 * the tensor product of (d-1)-dimensional quadrature rule over 
 * the unit simplex and interval and uses variant of Duffy Transformation 
 * to obtain the initial guess for quadrature over the d-dimensional unit simplex. 
 * When quadrature for dim1-dimensional simplex and dim2-dimensional simplices 
 * are computed by Node Elimiantion, the tensor product of dim1 and dim2 simplices 
 * is used as the final initial guess for SIMPLEXSIMPLEX. When Node Elimination for the dimensions assigned by a user
 * is finalized, the final quadrature and analysis  is printed to quadRule.txt and results.txt. 
*/
void ComputeSimplexSimplex(struct functionHandles *domainFunctions, struct quadratureParameters *params) {
  EstimateMemory(params);
  int dimSimplexMax1=params->dimension[0];
  int dimSimplexMax2=params->dimension[1];
  char flag[0];
  int dimSimplex1, dimSimplex2;
  int dim=params->totalDimension;  
  int p=params->degreeOfPrecision;
  int n;
  int nNodesPrevSimplex1, nNodesCurSimplex1, nNodesNewSimplex1;
  int nNodesPrevSimplex2, nNodesCurSimplex2, nNodesNewSimplex2;
  int nNodesCurSimplexSimplex, nNodesNewSimplexSimplex;
  double res;
  double *wInitial, *xInitial, *xNewSimplex1, *wNewSimplex1, *xNewSimplex2, *wNewSimplex2,  \
  *xNewSimplexSimplex, *wNewSimplexSimplex;
  struct eliminationHistory *history;
  //input parameters for generating Gaussian rules
  double p1=0.0, p2=0.0;
  int flip, i,j,d;
  //flip coordinates in a favourable order
  flip=0;
  if(dimSimplexMax1<dimSimplexMax2) {
    int temp=dimSimplexMax2;
    dimSimplexMax2=dimSimplexMax1;
    dimSimplexMax1=temp;
    flip=1;
  }
  //implement recursive scheme for computing the initial guess and run Node Elimination algorithm 
  for(d=2; d<=dimSimplexMax1;d++) {

    //set up parameters
    dimSimplex1=d;
    SetSimplex(domainFunctions);
  	SetParams(dimSimplex1, 0, 0, p, params);
    if(d==2) {
      n=floor(p/2)+1+floor(p/8);
      double x[n], w[n];
      Jacobi(n, p1, p2, x,w);     //generate Gaussian nodes and weights in 1-d on [0,1]
      nNodesPrevSimplex1=n;
      nNodesCurSimplex1=n*n;
      nNodesNewSimplex1=nNodesCurSimplex1;	
      AllocateMemory(d, nNodesCurSimplex1, &xInitial, &wInitial, &xNewSimplex1, &wNewSimplex1, &history); 
      AddLineSimplex(d, n, n, x,x , w, w, xInitial, wInitial); 
    }else if(d>2) {
	  //set up parameters
      n=floor(p/2)+2+floor(p/6);
      double x[n], w[n];
      Jacobi(n, p1, p2, x,w);           //generate Gaussian nodes and weights in 1-d on [0,1]
      nNodesPrevSimplex1=nNodesNewSimplex1;
      nNodesCurSimplex1=nNodesNewSimplex1*n;	
      nNodesNewSimplex1=nNodesCurSimplex1;

      ReallocateMemory(d, nNodesCurSimplex1, &xInitial, &wInitial, &xNewSimplex1, &wNewSimplex1, &history);
      AddLineSimplex(d, n, nNodesPrevSimplex1, x,xNewSimplex1, w, wNewSimplex1, xInitial, wInitial); 
    }
    //perform Node Elimiantion
    NodeElimination(nNodesCurSimplex1 ,xInitial,wInitial, &nNodesNewSimplex1,  xNewSimplex1, wNewSimplex1, domainFunctions, params, history);    
    if(d==dimSimplexMax2) {
      nNodesNewSimplex2=nNodesNewSimplex1;
      xNewSimplex2=malloc(d*nNodesNewSimplex2*sizeof(double));
      wNewSimplex2=malloc( nNodesNewSimplex2*sizeof(double)); 
      CopyNodes(nNodesNewSimplex2, d, xNewSimplex2, xNewSimplex1);
      CopyWeights(nNodesNewSimplex2,  wNewSimplex2, wNewSimplex1);
    }
  }
  //set up parameters for SIMPLEXSIMPLEX
  nNodesCurSimplexSimplex=nNodesNewSimplex1*nNodesNewSimplex2;
  nNodesNewSimplexSimplex=nNodesCurSimplexSimplex;
  xInitial=realloc(xInitial, nNodesCurSimplexSimplex*dim*sizeof(double));
  wInitial=realloc(wInitial, nNodesCurSimplexSimplex*sizeof(double));
  xNewSimplexSimplex=malloc(dim*nNodesCurSimplexSimplex*sizeof(double));
  wNewSimplexSimplex=malloc(nNodesCurSimplexSimplex*sizeof(double));
  SetSimplexSimplex(domainFunctions);
 
  SetParams(dimSimplexMax1, dimSimplexMax2,0, p, params);
  int dim1, dim2, counter;
  //compute tensor product of simplex1 x simplex2 to obtain SIMPLEXSIMPLEX
  for(i=0; i<nNodesNewSimplex1; i++) {
    for(j=0; j<nNodesNewSimplex2; j++) {

      for( dim1=0; dim1<dimSimplexMax1; dim1++) {
        xInitial[i*nNodesNewSimplex2*dim+dim*j+dim1]=xNewSimplex1[i*dimSimplexMax1+dim1];
      }		
      counter=0;
      for(dim2=dimSimplexMax1; dim2<dim; dim2++) {
        xInitial[i*nNodesNewSimplex2*dim+dim*j+dim2]=xNewSimplex2[j*dimSimplexMax2+counter];
        counter++;
      }	
    }
  }
  WeightsTensor(nNodesNewSimplex1, nNodesNewSimplex2, wNewSimplex1, wNewSimplex2, wInitial);
  free(xNewSimplex1);
  free(xNewSimplex2);
  free(wNewSimplex1);
  free(wNewSimplex2);
  history->nodesTotal=realloc(history->nodesTotal, nNodesCurSimplexSimplex*sizeof(int));  
  history->successNodeIndex=realloc(history->successNodeIndex, nNodesCurSimplexSimplex*sizeof(int));    
  history->successNewtonIterations=realloc(history->successNewtonIterations, nNodesCurSimplexSimplex*sizeof(int)); 
  //perform Node Elimination
  NodeElimination(nNodesCurSimplexSimplex ,xInitial,wInitial, &nNodesNewSimplexSimplex,  xNewSimplexSimplex, wNewSimplexSimplex, domainFunctions, params, history);  
  //rearrange coordinates of nodes in order assigned by the user
  if(flip==1) {
    double xTempInitial[dim];
    double xTemp[dim];
    for(i=0; i<nNodesCurSimplexSimplex; i++) {
      for(d=0; d<dim; d++) {
        xTempInitial[d]=xInitial[i*dim+d];
      }
      for(d=0; d<dimSimplexMax1; d++) {
        xInitial[i*dim+d]=xTempInitial[dimSimplexMax2+d];
      }	
      for(d=0; d<dimSimplexMax2; d++) {
        xInitial[i*dim+dimSimplexMax1+d]=xTempInitial[d];
      }	
    }

    for(i=0; i<nNodesNewSimplexSimplex; i++) {
      for(d=0; d<dim; d++) {
        xTemp[d]=xNewSimplexSimplex[i*dim+d];
      }
      for(d=0; d<dimSimplexMax2; d++) {
        xNewSimplexSimplex[i*dim+d]=xTemp[dimSimplexMax1+d];
      }	
      for(d=0; d<dimSimplexMax1; d++) {
        xNewSimplexSimplex[i*dim+dimSimplexMax2+d]=xTemp[d];
      }	
    }
    SetParams(dimSimplexMax2, dimSimplexMax1, 0, p, params);
  }  
  //print results to results.txt and quadRule.txt
  res=TestIntegral(nNodesNewSimplexSimplex, xNewSimplexSimplex,wNewSimplexSimplex, domainFunctions, params);
  printf("Final residual of the sum of all basis function integrals =%.16f \n", res); 
  Output(xInitial, wInitial, nNodesCurSimplexSimplex,xNewSimplexSimplex, wNewSimplexSimplex, nNodesNewSimplexSimplex, domainFunctions, params,  history, "SIMPLEXSIMPLEX"); 
  DumpCubatureRule(xNewSimplexSimplex, wNewSimplexSimplex, nNodesNewSimplexSimplex, params, "SIMPLEXSIMPLEX");
  FreeMemory(xInitial, wInitial, xNewSimplexSimplex, wNewSimplexSimplex, history);
} /* end ComputeSimplexSimplex */

/* ComputeCubeSimplexSimplex
 * domainFunctions: structure that contains function handles 
 * params: structure that contains quadrature parameters
 * 
 * Sets up the initial guess for the Node Elimination algorithm. 
 * Initial guess is computed via recursive scheme that reuses 
 * quadrature rules obtained in lower dimensions. It computes 
 * the tensor product of (d-1)-dimensional quadrature rule over 
 * the unit simplex and interval and uses variant of Duffy Transformation 
 * to obtain the initial guess for quadrature over the d-dimensional unit simplex. 
 * When quadrature for dim1-dimensional simplex and dim2-dimensional simplices 
 * are computed by Node Elimiantion, the tensor product of dim1 and dim2 simplices 
 * is used as the initial guess for SIMPLEXSIMPLEX. Then, tensor product of SIMPLEXSIMPLEX 
 * and interval is computed to obtain initial guess for I x SIMPLEXSIMPLEX, and so on. When Node Elimination
 * for the dimensions assigned by a user is finalized, the final quadrature and analysis 
 * is printed to quadRule.txt and results.txt. 
*/
void ComputeCubeSimplexSimplex(struct functionHandles *domainFunctions, struct quadratureParameters *params){
  EstimateMemory(params);
  int dimCubeMax=params->dimension[0];
  int dimSimplexMax1=params->dimension[1];
  int dimSimplexMax2=params->dimension[2];
  char flag[0];
  int dimCube, dimSimplex1, dimSimplex2;
  int dim=params->totalDimension;  
  int p=params->degreeOfPrecision;
  int n;
  int nNodesPrevSimplex1, nNodesCurSimplex1, nNodesNewSimplex1;
  int nNodesPrevSimplex2, nNodesCurSimplex2, nNodesNewSimplex2;
  int nNodesPrev, nNodesCur, nNodesNew;
  double res;
  double *wInitial, *xInitial, *xNewSimplex1, *wNewSimplex1, *xNewSimplex2, *wNewSimplex2,  \
  *xNew, *wNew;
  struct eliminationHistory *history; 
  //input parameters for generating Gaussian rules
  double p1=0.0, p2=0.0;
  int flip, i,j,d;
  //flip coordinates in a favourable order
  flip=0;
  if(dimSimplexMax1<dimSimplexMax2) {
    int temp=dimSimplexMax2;
    dimSimplexMax2=dimSimplexMax1;
    dimSimplexMax1=temp;
    flip=1;
  }
  //implement recursive scheme for computing the initial guess and run Node Elimination algorithm 
  for(d=2; d<=dimSimplexMax1;d++) {
    //set up parameters
    dimSimplex1=d;
    SetSimplex(domainFunctions);
    SetParams(dimSimplex1, 0, 0, p, params);
    if(d==2) {
      n=floor(p/2)+1+floor(p/8);
      double x[n], w[n];
      Jacobi(n, p1, p2, x,w);     //generate Gaussian nodes and weights in 1-d on [0,1]
      nNodesPrevSimplex1=n;
      nNodesCurSimplex1=n*n;
      nNodesNewSimplex1=nNodesCurSimplex1;	
      AllocateMemory(d, nNodesCurSimplex1, &xInitial, &wInitial, &xNewSimplex1, &wNewSimplex1, &history); 
      AddLineSimplex(d, n, n, x,x , w, w, xInitial, wInitial); 
		
    }else if(d>2) {
      n=floor(p/2)+2+floor(p/6);
      double x[n], w[n];
      Jacobi(n, p1, p2, x,w);           //generate Gaussian nodes and weights in 1-d on [0,1]
      nNodesPrevSimplex1=nNodesNewSimplex1;
      nNodesCurSimplex1=nNodesNewSimplex1*n;	
      nNodesNewSimplex1=nNodesCurSimplex1;
      ReallocateMemory(d, nNodesCurSimplex1, &xInitial, &wInitial, &xNewSimplex1, &wNewSimplex1, &history); 
      AddLineSimplex(d, n, nNodesPrevSimplex1, x,xNewSimplex1, w, wNewSimplex1, xInitial, wInitial); 
    }
	 
    NodeElimination(nNodesCurSimplex1 ,xInitial,wInitial, &nNodesNewSimplex1,  xNewSimplex1, wNewSimplex1, domainFunctions, params, history);    
    if(d==dimSimplexMax2) {
      nNodesNewSimplex2=nNodesNewSimplex1;
      xNewSimplex2=malloc(d*nNodesNewSimplex2*sizeof(double));
      wNewSimplex2=malloc( nNodesNewSimplex2*sizeof(double)); 
      CopyNodes(nNodesNewSimplex2, d, xNewSimplex2, xNewSimplex1);
      CopyWeights(nNodesNewSimplex2,  wNewSimplex2, wNewSimplex1);
    }
  }

  int twoDims=dimSimplexMax1+dimSimplexMax2;
  nNodesCur=nNodesNewSimplex1*nNodesNewSimplex2;
  nNodesNew=nNodesCur;
  xInitial=realloc(xInitial, nNodesCur*twoDims*sizeof(double));
  wInitial=realloc(wInitial, nNodesCur*sizeof(double));
  xNew=malloc(dim*nNodesCur*sizeof(double));
  wNew=malloc(nNodesCur*sizeof(double));
  SetSimplexSimplex(domainFunctions);
  SetParams(dimSimplexMax1, dimSimplexMax2,0, p, params);
  int dim1, dim2, counter;

  for(i=0; i<nNodesNewSimplex1; i++) {
    for(j=0; j<nNodesNewSimplex2; j++) {

      for( dim1=0; dim1<dimSimplexMax1; dim1++) {
        xInitial[i*nNodesNewSimplex2*twoDims+twoDims*j+dim1]=xNewSimplex1[i*dimSimplexMax1+dim1];
      }		
      counter=0;
      for(dim2=dimSimplexMax1; dim2<twoDims; dim2++) {
        xInitial[i*nNodesNewSimplex2*twoDims+twoDims*j+dim2]=xNewSimplex2[j*dimSimplexMax2+counter];
        counter++;
      }	
    }
  }
		
  WeightsTensor(nNodesNewSimplex1, nNodesNewSimplex2, wNewSimplex1, wNewSimplex2, wInitial);
  free(xNewSimplex1);
  free(xNewSimplex2);
  free(wNewSimplex1);
  free(wNewSimplex2);
  history->nodesTotal=realloc(history->nodesTotal, nNodesCur*sizeof(int));   
  history->successNodeIndex=realloc(history->successNodeIndex, nNodesCur*sizeof(int));     
  history->successNewtonIterations=realloc(history->successNewtonIterations, nNodesCur*sizeof(int)); 
  NodeElimination(nNodesCur ,xInitial,wInitial, &nNodesNew,  xNew, wNew, domainFunctions, params, history);    

  //rearrange coordinates of nodes in order assigned by the user	
  if(flip==1) {
    double xTemp[twoDims];

    for(i=0; i<nNodesNew; i++) {
      for(d=0; d<twoDims; d++) {
        xTemp[d]=xNew[i*twoDims+d];
      }
      for(d=0; d<dimSimplexMax2; d++) {
        xNew[i*twoDims+d]=xTemp[dimSimplexMax1+d];
      } 	
      for(d=0; d<dimSimplexMax1; d++) {
        xNew[i*twoDims+dimSimplexMax2+d]=xTemp[d];
      }	
    }
    int temp=dimSimplexMax2;
    dimSimplexMax2=dimSimplexMax1;
    dimSimplexMax1=temp;
  } 

  n=floor(p/2)+1;
  double x[n], w[n];
  Jacobi(n, p1, p2, x,w);     //generate Gaussian nodes and weights in 1-d on [0,1]
  SetCubeSimplexSimplex(domainFunctions);

  //implement recursive scheme for computing the initial guess and run Node Elimination algorithm 			
  for(d=1; d<=dimCubeMax; d++) {
    //set up parameters
    SetParams(dimSimplexMax1, dimSimplexMax2, d,p, params);
    nNodesPrev=nNodesNew;
    nNodesCur=nNodesNew*n;	
    nNodesNew=nNodesCur;
    ReallocateMemory(twoDims+d, nNodesCur, &xInitial, &wInitial, &xNew, &wNew, &history); 
    AddLineLast(twoDims+d, n, nNodesPrev, x,xNew, w, wNew, xInitial, wInitial); 
    //perform Node Elimination 
    NodeElimination(nNodesCur ,xInitial,wInitial, &nNodesNew,  xNew, wNew, domainFunctions, params, history);    
  }

  res=TestIntegral(  nNodesNew,  xNew, wNew, domainFunctions, params);
  printf("Final residual of the sum of all basis function integrals =%.16f \n", res); 
  double xTempInitial[dim];
  double xTemp[dim];
  //rearrange coordinates of nodes in order assigned by the user
  for(i=0; i<nNodesCur; i++) {
    for(d=0; d<dim; d++) {
      xTempInitial[d]=xInitial[i*dim+d];
    }
    for(d=0; d<dimCubeMax; d++) {
      xInitial[i*dim+d]=xTempInitial[twoDims+d];
    }	
    for(d=0; d<twoDims; d++) {
      xInitial[i*dim+dimCubeMax+d]=xTempInitial[d];
    }	
  }

  for(i=0; i<nNodesNew; i++) {
    for(d=0; d<dim; d++) {
      xTemp[d]=xNew[i*dim+d];
    }
    for(d=0; d<dimCubeMax; d++) {
      xNew[i*dim+d]=xTemp[twoDims+d];
    }	
    for(d=0; d<twoDims; d++) {
      xNew[i*dim+dimCubeMax+d]=xTemp[d];
    }	
  } 

  //print results to results.txt and quadRule.txt
  Output(xInitial, wInitial, nNodesCur,xNew, wNew, nNodesNew, domainFunctions, params,  history, "CUBESIMPLEXSIMPLEX"); 
  DumpCubatureRule(xNew, wNew, nNodesNew, params, "CUBESIMPLEXSIMPLEX");
  FreeMemory(xInitial, wInitial, xNew, wNew, history);

} /* end ComputeCubeSimplexSimplex */


/* EstimateMemory 
 * Computes an estimate of memory consumption by
 * computing an approximate number of entries needed
 * by jacobian in Newton's method for the initial guess 
 * of final dimension of interest. The memory consumption 
 * is dominated by jacobian.    
*/	 
static void EstimateMemory(struct quadratureParameters *params) {
  char flag[0];
  unsigned long int memoryEstimate;
  memoryEstimate= 2*8*pow(params->numberOfBasisFunctions, 2)*2;

  if(memoryEstimate> 8*pow(10, 9)){
    printf("Warning: Memory consumption estimated to eventually exceed %i gigabytes. \
    Enter character 'Y' if you would like to continue. Exiting otherwise. \n", (int)(memoryEstimate/(pow(10,9))));
    scanf("%c", flag);
    if(flag[0]!='Y') {
      printf("Exiting program \n");
      exit(-1);
    }
  }

} /* end EstimateMemory */


static void FreeMemory(double *xInitial, double *wInitial, double *xNew, double *wNew, struct eliminationHistory *history){
  free(xInitial);
  free(wInitial);
  free(xNew);
  free(wNew);
  free(history->nodesTotal);  
  free(history->successNodeIndex);      
  free(history->successNewtonIterations); 
  free(history);
} /* end FreeMemory */


static void AllocateMemory(int dim, int nNodes, double **xInitial, double **wInitial, double **xNew, double **wNew, struct eliminationHistory **history){
  *xInitial=malloc(nNodes*dim*sizeof(double));
  *wInitial=malloc(nNodes*sizeof(double));
  *xNew=malloc(dim*nNodes*sizeof(double));
  *wNew=malloc(nNodes*sizeof(double)); 
  (*history)=(struct eliminationHistory*)malloc(sizeof(struct eliminationHistory)); 
  (*history)->nodesTotal=malloc(nNodes*sizeof(int));  
  (*history)->successNodeIndex=malloc(nNodes*sizeof(int));   
  (*history)->successNewtonIterations=malloc(nNodes*sizeof(int));
} /* end ReallocateMemory */


static void ReallocateMemory(int dim, int nNodes, double **xInitial, double **wInitial, double **xNew, double **wNew, struct eliminationHistory **history){
  *xInitial=realloc(*xInitial, nNodes*dim*sizeof(double));
  *wInitial=realloc(*wInitial, nNodes*sizeof(double));
  *xNew=realloc(*xNew, dim*nNodes*sizeof(double));
  *wNew=realloc(*wNew, nNodes*sizeof(double)); 
  (*history)->nodesTotal=realloc((*history)->nodesTotal, nNodes*sizeof(int));  
  (*history)->successNodeIndex=realloc((*history)->successNodeIndex, nNodes*sizeof(int));   
  (*history)->successNewtonIterations=realloc((*history)->successNewtonIterations, nNodes*sizeof(int));
} /* end ReallocateMemory */


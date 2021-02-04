/* Arkadijs Slobodkins
 * SMU Mathematics
 * February 2021
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../headers/TestIntegral.h"
#include "../headers/BasisIndices.h"
#include "../headers/IntegralsOfBasisFunctions.h"
#include "../headers/Structs.h"

/* TestIntegral
 * k: number of nodes
 * x: quadrature nodes
 * w: quadrature weights
 * functions: structure that contains function handles for the respective domain
 * params: structure that contains quadrature parameters
 * 
 * Approximates the difference between the sum of all analytical integrals of polynomial basis functions 
 * and sum of their quadrature approximations. The difference sum (residual) is returned. 
 */
double TestIntegral(const int k, const double *const x, const double *w, struct functionHandles *functions, struct quadratureParameters *params){
  int p=params->degreeOfPrecision;
  int m=params->numberOfBasisFunctions;
  int dim1=params->dimension[0];
  int dim2=params->dimension[1];
  int dim3=params->dimension[2];
  int dim=dim1+dim2+dim3;  
  int i,j;
  double In=0;
  double Ie=0;;
  int basisSize=BasisSize(dim, p);
  double *phi=malloc(basisSize*sizeof(double));
  double *integrals=malloc(basisSize*sizeof(double));
  //approximate integrals of basis functions
  for(i=0;i<k;i++) {
    functions->EvalBasis(&x[dim*i], phi, params);
    for(j=0; j<basisSize;j++) {
      In=In+phi[j]*w[i];        
    }
  }
  functions->IntegralsOfBasisFunctions(integrals, params);
  Ie=integrals[0];    //integral of only the first function is nonzero, because all bases are orthogonal
  free(phi);
  free(integrals);
  double res=Ie-In;
  return res;
} /* end TestIntegral */


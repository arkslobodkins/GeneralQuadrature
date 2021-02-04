/* Arkadijs Slobodkins
 * SMU Mathematics
 * February 2021
 */

#include <stdlib.h>
#include "../headers/GetFunction.h"
#include "../headers/Structs.h"

/* GetFunction
 * numNodes: number of nodes
 * X: array of nodes 
 * w: array of weights
 * f: function to be computed 
 * functions: structure that contains function handles for the respective domain
 * params: structure that contains quadrature parameters
 * 
 * Computes an evaluates function f=b-X*W 
 * which will employ Newton's method.
 */
void GetFunction(int numNodes, double *X, double *w, double *f, struct functionHandles *functions, struct quadratureParameters *params) {      
  int i,j;
  int m=params->numberOfBasisFunctions;
  int dim=params->totalDimension;
  double *xtemp;  
  double *b=malloc(m*sizeof(double));
  double *phi=malloc(m*sizeof(double));

  functions->IntegralsOfBasisFunctions(b, params);  
  for(i=0;i<m;i++){
    f[i]=-1.0*b[i];
  }
  for(j=0;j<numNodes;j++) {
    xtemp=&X[dim*j];
    functions->EvalBasis(xtemp,phi, params);     
    for(i=0;i<m;i++){
      f[i]=f[i]+w[j]*phi[i];        
    }
  }
  free(b);
  free(phi);
} /*end GetFunction */

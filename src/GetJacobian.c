#include <stdlib.h>
#include "../headers/GetJacobian.h"
#include "../headers/Structs.h"

/* GetJacobian
 * numNodes: number of nodes
 * x: array of nodes
 * w: array of weights
 * jacobian: array for computing Jacobian
 * functions: structure that contains function handles for the respective domain
 * params: structure that contains quadrature parameters
 * 
 * Computes an returns Jacobian of size mx3*k 
 * of a vector X*w, where each row of X corresponds to 
 * a distinct polynomial of orthogonal basis, and each 
 * column of X and entries in w correspond to ith node and weight
 * of the quadrature
 */
void GetJacobian( int numNodes , double *x, double *w,  double *jacobian   , struct functionHandles *functions, struct quadratureParameters *params) {
  int i,j,d;            
  int numRows=params->numberOfBasisFunctions;
  int numColumns=numNodes;
  int dim=params->totalDimension;
  double *phi=malloc(numRows*sizeof(double));
  double *phiPrime=malloc(numRows*dim*sizeof(double));
  double *curNode;
  
  //compute analytical derivatives with respect to all weights and node components
  for(j=0;j<numColumns;j++){    
    curNode=&x[dim*j];
    functions->EvalBasis(&curNode[0], &phi[0], params);            
    functions->EvalBasisDerivatives(&curNode[0], &phiPrime[0], params);  
    for(i=0;i<numRows;i++){       
      jacobian[(dim+1)*numColumns*i+(dim+1)*j]=phi[i];
      for(d=0; d<dim;d++) {
        jacobian[(dim+1)*numColumns*i+(dim+1)*j+d+1]=phiPrime[i+d*numRows]*w[j];        
      }
    }
  }
  free(phi);
  free(phiPrime);
} /*end GetJacobian */


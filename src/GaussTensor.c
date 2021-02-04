/* Arkadijs Slobodkins
 * SMU Mathematics
 * February 2021
 */
#include "../headers/GaussTensor.h"

/* NodesTensor
 * n1: number of nodes of x1
 * n2: number of nodes of x2
 * x1: vector used for computing tensor product 
 * x2: vector used for computing tensor product 
 * X: node tensor product of x1 and x2 
 *  
 * Computes and returns 2-d tensor product of nodes 
 */
void NodesTensor(int n1,int n2, double *x1, double *x2, double *X)  {
  int i,j;
  int N=n1*n2;
        
  for(i=0;i<n1;i++) {
    for(j=0;j<n2;j++) {                 
      X[2*(n2*i+j)]=x1[i];
      X[2*(n2*i+j)+1]=x2[j];            
    }
  }
  
} /*end NodesTensor */


/* WeightsTensor
 * n1: size of w1
 * n2: size of w2
 * w1: vector used for computing tensor product
 * w2: vector used for computing tensor product
 * W: weight tensor product of w1 and w2
 *
 * Computes and returns 2-d tensor product of weights
 */
void WeightsTensor(int n1, int n2, double *w1,double *w2, double *W) {
  int i,j;
        
  for(i=0;i<n1;i++) {
    for(j=0;j<n2;j++) {                                                                 
      W[n2*i+j]=w1[i]*w2[j];                                                    
    }
  } 

} /* end WeightsTensor */  


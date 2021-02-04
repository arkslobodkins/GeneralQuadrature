/* Arkadijs Slobodkins
 * SMU Mathematics
 * February 2021
 */

#include "../headers/StoreCopy.h"

void CopyNodes(int n, int dim, double *x, double *y) {
  int i,j;

  for(i=0;i<n;i++) {  
    for(j=0;j<dim;j++){
      x[i*dim+j]=y[i*dim+j];  
    }
  }
  
} /* end CopyNodes */

void CopyWeights(int n, double *w1, double *w2) {
  int i;

  for(i=0; i<n; i++) {
    w1[i]=w2[i];
  }

} /* end CopyWeights */

void CopyNodesAndWeights(int n, int dim, double *z1, double *z2) {
  int i,d;

  for(i=0; i<n; i++) {
    for(d=0; d<dim+1; d++) {
      z1[i*(dim+1)+d]=z2[i*(dim+1)+d];
    }
  }

} /* end CopyNodesAndWeights */

void UnrollNodesAndWeights(int k, int dim, double *xVector, double *x, double *w) {
  int i,j;

  for(j=0; j<k;j++) {
    xVector[j*(dim+1)]=w[j];
    for(i=0; i<dim; i++) {
      xVector[i+j*(dim+1)+1]=x[i+dim*j];
    }
  }

} /* end UnrollNodesAndWeights */

void RollNodesAndWeights(int k, int dim, double *xVector, double *x, double *w) {
  int i,j;

  for(j=0; j<k;j++) {
    w[j]=xVector[j*(dim+1)];
    for(i=0; i<dim; i++) {
      x[i+dim*j]=xVector[i+j*(dim+1)+1];
    }
  }

} /* end RollNodesAndWeights */

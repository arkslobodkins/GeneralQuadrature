#include <math.h>
#include <stdlib.h>
#include "../headers/GeneralizedGaussTensor.h"
#include "../headers/GaussTensor.h"

/* GeneralizedNodesTensor */
void GeneralizedNodesTensor(int n, int dim, double *x1, double *x2, double *X) {
  int i,j,k,s;
  int nFinal=pow(n, dim);
  double *xTemp1=malloc(2*n*n*sizeof(double));
  double *xTemp2=malloc(nFinal*dim*sizeof(double)); 
  double *xTemp3=malloc(nFinal*sizeof(double));
  int counter=pow(n, dim-2);
  
  for(i=0;i<counter;i++) {    
    NodesTensor(n,n, x1,x2, xTemp1);  
    
    for(k=0;k<n*n;k++) {
      X[ dim-2+dim*i*n*n+ dim*k ]= xTemp1[2*k];
      X[ dim-1+dim*i*n*n +dim*k]= xTemp1[2*k+1];
    }
  }

  for(j=dim-3;j>=0;j--) {
    counter=pow(n,j);
    
    for(i=0; i<counter;i++) { 
      for(s=0;s<pow(n,dim-j-1);s++) { 
        xTemp3[s]= X[ s*dim+j+1];
      }
      NodesTensor(n, pow(n,dim-j-1), x1, xTemp3, xTemp2 ) ;
                   
      for(k=0; k<pow(n, dim-j); k++) {
        X[ j+ (i)*(int)pow(n,dim-j)*dim+dim*k ]= xTemp2[2*k];
      }
    }

  }

  free(xTemp1);
  free(xTemp2);
  free(xTemp3);
} 


void GeneralizedWeightsTensor (int n, int dim, double *w1,double *w2, double *W) {
  int i,j,k;
  int nFinal=pow(n,dim);
  double *wTemp1=malloc(nFinal*sizeof(double));
  double *wTemp2=malloc(nFinal*sizeof(double));
  int counter=pow(n, dim-2);  

  for(i=0;i<counter;i++) {    
    WeightsTensor(n,n,w1,w2, wTemp1); 
    for(k=0;k<n*n;k++) {
      W[ i*n*n+k]=wTemp1[k];
    }   
  }
  
  for(j=dim-3;j>=0;j--) {
    counter=pow(n,j);
    for(i=0; i<counter;i++) {   
      WeightsTensor(n, pow(n,dim-j-1), w1, W , wTemp2) ;        
      for(k=0; k<pow(n, dim-j-1)*n; k++) {
        W[i*(int)pow(n,dim-j)+k]= wTemp2[k];
      } 
    } 
  } 

free(wTemp1);
free(wTemp2);
}

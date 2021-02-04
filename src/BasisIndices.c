/* Arkadijs Slobodkins
 * SMU Mathematics
 * February 2021
 */

#include <stdlib.h>
#include <math.h>
#include "../headers/BasisIndices.h"

/* BasisIndices
 * dim: dimension of polynomial basis
 * deg: degree of polynomial basis
 * f: array to store the indices 
 * 
 * Recursively computes an array of power multindices for 
 * generating polynomial basis of dimension 
 * dim and degree deg. 
*/
void BasisIndices(int  dim, int deg, int *f) {
  int i,j,k,d;
  int size;

  int counter=0;
  if(dim!=1) {
    //compute basis indices using recursion
    for(j=2;j<=dim;j++) {	
      for(i=0; i<=deg; i++) {
	
        size=BasisSize(j-1, deg-i);
        int* recursiveF=malloc(size*(j-1)*sizeof(int));
        BasisIndices(j-1, deg-i, recursiveF);
					
        if(j==dim) {			
          for(k=0;k<size;k++) {
            for(d=0; d<j-1;d++) {		
              f[counter+k*dim+d]=recursiveF[k*(j-1)+d];	
            }	
            f[counter+k*dim+j-1]=i;		
          }
          counter=counter+dim*size;
        } //end if		
        free(recursiveF); 
      }
    }
	
  }//end if

  if(dim==1) {
    int index=0;
    for(i=0; i<=deg;i++){	
      f[index]=i;
      index++;	
    }	
  } 

} /* end BasisIndices */


/* BasisSize 
 * Computes the minimum number of basis
 * functions to generate polynomial basis 
 * of dimension dim and degree deg. 
*/
int BasisSize(int dim, int deg) {
  int first=deg+1;
  int last=deg+dim;
  int product=1;
  int counter=1;
  int i;
	
  for(i=first; i<=last;i++) {
    product=product*(i)/counter;	
    counter++;
  }
	
  if(dim==0) {
    product=0;
  }
	
  return product;
} /* end BasisSize */

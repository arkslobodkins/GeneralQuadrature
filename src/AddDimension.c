/* Arkadijs Slobodkins
 * SMU Mathematics
 * February 2021
 */

#include "../headers/AddDimension.h"
#include "../headers/GaussTensor.h"
/* AddLine
 * dim: dimension of the new quadrature 
 * n: number of nodes on the interval 
 * nOld: number of nodes on (dim-1)-dimensional domain
 * x: nodes for the interval
 * xOld: nodes for the (dim-1)-dimensional domain
 * w: weights for the interval 
 * wOld: weights for the (dim-1)-dimensional domain
 * xNew: nodes of the new (dim)-dimensional quadrature 
 * wNew: weights of the new (dim)-dimensional quadrature
 * 
 * Computes tensor product of arbitrary quadrature over (dim-1)-dimensional 
 * domain D and 1-dimensional interval [0,1]. Quadrature of dimension dim 
 * and degree over [0,1] x D is generated. Stores interval quadrature nodes 
 * as the first coordinate.  
 */
void AddLine(int dim, int n, int nOld,double *x,double *xOld, double *w, double *wOld, double *xNew, double *wNew)  {
  int i,j,d;
	//compute tensor product for nodes
  for(i=0; i<n; i++) {
    for(j=0; j<nOld;j++) {
      xNew[i*dim*nOld+dim*j]=x[i];	
      for(d=1; d<dim; d++) {
        xNew[i*dim*nOld+dim*j+(d)]=xOld[(dim-1)*j+d-1];			
      }
    }
  }				 	 	 

	//compute tensor product for weights
  WeightsTensor(n, nOld, w,wOld, wNew); 	
} /* end AddLine */


/* AddLineLast
 * dim: dimension of the new quadrature 
 * n: number of nodes on the interval 
 * nOld: number of nodes on (dim-1)-dimensional domain
 * x: nodes for the interval
 * xOld: nodes for the (dim-1)-dimensional domain
 * w: weights for the interval 
 * wOld: weights for the (dim-1)-dimensional domain
 * xNew: nodes of the new (dim)-dimensional quadrature 
 * wNew: weights of the new (dim)-dimensional quadrature
 * 
 * Computes tensor product of arbitrary quadrature over (dim-1)-dimensional 
 * domain D and 1-dimensional interval [0,1]. Quadrature of dimension dim 
 * and degree over D x [0,1] is generated. Stores interval quadrature nodes 
 * as the last coordinate.
*/
void AddLineLast(int dim, int n, int nOld, double *x,double *xOld, double *w, double *wOld, double *xNew, double *wNew)  {
  int i,j,d;

  //compute tensor product for nodes
  for(i=0; i<n; i++) {
    for(j=0; j<nOld;j++) {
      xNew[i*dim*nOld+dim*j+dim-1]=x[i];
      for(d=0; d<dim-1; d++) {
        xNew[i*dim*nOld+dim*j+(d)]=xOld[(dim-1)*j+d];			
      }
    }
  }		 	 	 

  //compute tensor product for weights
  WeightsTensor(n,nOld, w,wOld, wNew); 	
} /* end AddLineLast */


/* AddLineLastSimplex
 * dim: dimension of the new quadrature 
 * n: number of nodes on the interval 
 * nOld: number of nodes on (dim-1)-dimensional simplex
 * x: nodes for the interval
 * xOld: nodes for the (dim-1)-dimensional simplex
 * w: weights for the interval 
 * wOld: weights for the (dim-1)-dimensional simplex
 * xNew: nodes of the new (dim)-dimensional simplex 
 * wNew: weights of the new (dim)-dimensional simplex
 * 
 * Computes tensor product of unit (dim-1)-dimensional simplex 
 * and 1-dimensional interval [0,1] and maps the product to unit 
 * simplex of dimension dim using Duffy transformation. 
 */
void AddLineSimplex(int dim, int n, int nOld, double *x,double *xOld, double *w, double *wOld, double *xNew, double *wNew) {
  int nNew=nOld*n;
  int i,j,d;
  //compute tensor product for nodes
  for(i=0; i<n; i++) {
    for(j=0; j<nOld;j++) {
      xNew[i*dim*nOld+dim*j]=x[i];	
      for(d=1; d<dim; d++) {
        xNew[i*dim*nOld+dim*j+(d)]=xOld[(dim-1)*j+d-1];
      }
    }
  }
  //compute tensor product for weights				 	 
  WeightsTensor(n, nOld, w,wOld, wNew); 

  //apply Duffy Transformation		
  for(i=0; i<nNew;i++) {
    for(d=1; d<dim; d++) {
      xNew[dim*i+d]=xNew[dim*i+d]*xNew[dim*(i)];
      wNew[i]=wNew[i]*xNew[dim*i];
    }	
  }	
	
} /*end AddLineSimplex */


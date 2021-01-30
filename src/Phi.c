#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../headers/Phi.h"
#include "../headers/LegendrePoly.h"
#include "../headers/BasisIndices.h"
#include "../headers/JacobiPoly.h"
#include "../headers/SetParams.h"
#include "../headers/Structs.h"

static int *basis;
static int dimPrev=0;
static int dimensionCount=0; 

/* PhiCube
 * x: point in (dim)-dimensional space
 * phi: array of basis functions evaluated at x
 * params: structure that contains quadrature parameters
 *
 * Generates orthogonal polynomial basis for the unit cube
*/
void PhiCube(const double *const x, double *phi, const struct quadratureParameters *params){
  int i,j,k, d;  
  int p=params->degreeOfPrecision;
  int m=params->numberOfBasisFunctions;
  int dimCube=params->dimension[0];
  int dim=params->totalDimension;   
  int order=p+1;  
  int size=BasisSize(dim,p);  //number of basis functions 
   
  //compute recursive basis indices only during the first call of the file for the same dimension
  if(dimPrev-dim!=0) { 
    if(dimensionCount>0) {
      free(basis);
    }
    basis=malloc((size*dim)*sizeof(int));
    BasisIndices(dim, p, basis);  
    dimensionCount++;
  }
  dimPrev=dim;

  double* legendre=malloc(order*dimCube*sizeof(double));
  double* dxlegendre=malloc(order*dimCube*sizeof(double));

  for(k=0; k<m; k++) {
    phi[k]=1.0;
  }
  for(d=0;d<dimCube;d++) {   
    legendrePoly(order, 2*x[d]-1, &legendre[d*order], &dxlegendre[d*order]); 
  }
  for(k=0;k<m;k++) {
    for(d=0;d<dimCube;d++) {
      phi[k]=phi[k]*legendre[basis[k*dim+d]+order*d];   
    }    
  }   
  free(legendre);
  free(dxlegendre);
} /* end PhiCube */

/* PhiPrimeCube
 * x: point in (dim)-dimensional space
 * phiPrime: array of derivatives of basis functions evaluated at x
 * params: structure that contains quadrature parameters
 *
 * Generates derivatives of orthogonal polynomial basis for the unit cube
*/
void PhiPrimeCube(const double *const x,double *phiPrime, const struct quadratureParameters *params) {
  int p=params->degreeOfPrecision;
  int m=params->numberOfBasisFunctions;
  int dim1=params->dimension[0];
  int dim2=params->dimension[1];
  int dim=dim1+dim2;   
  int i,j,k,d; 
  int order=p+1;
  int size=BasisSize(dim,p); 
  double *legendre=malloc(order*dim1*sizeof(double));
  double *dxlegendre=malloc(order*dim1*sizeof(double));  
  //compute recursive basis indices only during the first call of the file for the same dimension
  if(dimPrev-dim!=0) { 
    if(dimensionCount>0) {
      free(basis);
    }
    basis=malloc((size*dim)*sizeof(int));
    BasisIndices(dim, p, basis);  
    dimensionCount++;
  }
  dimPrev=dim;
 
  for(k=0; k<m*dim;k++) {
    phiPrime[k]=1.0;
  }
  for(d=0;d<dim1;d++) {   
    legendrePoly(order, 2*x[d]-1, &legendre[d*order], &dxlegendre[d*order]);  
  }
  for(d=0; d<dim; d++) {
    for(k=0;k<m;k++) {
  
      for(j=0;j<dim1;j++) {
        if(j!=d) {
          phiPrime[k+d*m]=phiPrime[k+d*m]*legendre[basis[k*dim+j]+order*j];  
        }
        if(j==d) {
          phiPrime[k+d*m]=2*phiPrime[k+d*m]*dxlegendre[basis[k*dim+j]+order*j];   
        }
      } 
    }
  }
  free(legendre);
  free(dxlegendre);
} /* end PhiPrimeCube */

/* PhiSimplex
 * x: point in (dim)-dimensional space
 * phi: array of basis functions evaluated at x
 * params: structure that contains quadrature parameters
 *
 * Generates orthogonal polynomial basis for the unit simplex
*/
void PhiSimplex(const double *const x, double *phi, const struct quadratureParameters *params){
  int i,j,k, d;  
  int p=params->degreeOfPrecision;
  int m=params->numberOfBasisFunctions;
  int dim1=params->dimension[0];
  int dim=params->dimension[0];
  int order=p+1; 
  int size=BasisSize(dim,p);  //number of basis functions

  //compute recursive basis indices only during the first call of the file for the same dimension
  if(dimPrev-dim!=0) { 
    if(dimensionCount>0) {
      free(basis);
    }
    basis=malloc((size*dim)*sizeof(int));
    BasisIndices(dim, p, basis);  
    dimensionCount++;
  }
  dimPrev=dim;
   
  for(k=0; k<m; k++) {
    phi[k]=1.0;
  }
 
  int dimCur=0,  xPower=0;
  double xFactor, xTemp=0, phiTemp=0, alpha=0;
  int *index;
  double *jacobiX, *polyX;
  double jacobi[order*order*(dim-1)];
  double *legendre=malloc(order*sizeof(double));
  double *dxlegendre=malloc(order*sizeof(double)); 
  if(dim==2) {
    index=&basis[0];

    for(j=0; j<order; j++) {
      JacobiPoly(order, 1.0-2*x[0], 2.0*j+1,0.0, &jacobi[j*order]);
    }
    int r,n; 

    legendrePoly(order, (2.0*x[1]-x[0])/x[0], legendre, dxlegendre);
    for(k=0;k<m;k++) { 
      r=*index;
      n=*(index+1);
      xFactor=1.0;
      for(i=0; i<r; i++) {
        xFactor=xFactor*x[0];
      }
      phi[k]=phi[k]*xFactor*jacobi[r*order+n]*legendre[r];  
      index=(index+dim);
    }
 
  } 
 
  if(dim>=3) {
    int *index;
    jacobiX=malloc((dim-1)*sizeof(double));  
    polyX=malloc((dim-1)*sizeof(double));
    legendrePoly(order, (2.0*x[dim-1]-x[dim-2])/x[dim-2], legendre, dxlegendre);
   
    for(d=0; d<dim-1;d++) {
      polyX[d]=1.0;
    } 
    for(d=0; d<dim-2; d++) {
      polyX[d]=x[dim-d-2]/x[dim-d-3];  
    }
    polyX[dim-2]=x[0];
    for(d=0; d<dim-1;d++) {
      jacobiX[d]=1.0;
    }

    for(d=1; d<dim;d++) {
      alpha=0.0;
      dimCur=d-1;
      if(d>=1&&d<dim-1) {
        jacobiX[dimCur]=jacobiX[dimCur]/x[dim-d-2];
        jacobiX[dimCur]=1.0-2.0*x[dim-d-1]*jacobiX[dimCur]; 
      }
      if(d==dim-1) {
        jacobiX[dimCur]=jacobiX[dimCur]*x[0];
        jacobiX[dimCur]=1.0-2.0*jacobiX[dimCur]; 
      }
      for(j=0; j<p+1; j++) {
        alpha=2*j+d;
        JacobiPoly(order, jacobiX[dimCur], alpha, 0.0,&jacobi[order*(d-1)*p+order*j]);
      }
    }

    for(k=0; k<m; k++) {
      index=&basis[k*dim];
      phiTemp=1.0;
      for(d=1;d<dim;d++) { 
        xPower=0;
        for(i=0; i<d; i++) {  
          xPower=xPower+*(index+i);
        } 
        xTemp=polyX[d-1];
        xFactor=1.0;
        for(i=0; i<xPower; i++) {
          xFactor=xFactor*xTemp;
        }
        phiTemp=phiTemp*jacobi[*(index+d)+order*(d-1)*p+order*xPower]*xFactor; 
      } 
      phi[k]=phi[k]*phiTemp*legendre[*index]; 
    } 
    free(jacobiX);
    free(polyX);
  }
 
  free(legendre);
  free(dxlegendre);
} /* end PhiSimplex */

/* PhiPrimeSimplex
 * x: point in (dim)-dimensional space
 * phiPrime: array of derivatives of basis functions evaluated at x
 * params: structure that contains quadrature parameters
 *
 * Generates derivatives of orthogonal polynomial basis for the unit simplex
*/
void PhiPrimeSimplex(const double *const x,double *phiPrime, const struct quadratureParameters *params) {
  int p=params->degreeOfPrecision;
  int m=params->numberOfBasisFunctions;
  int dim=params->totalDimension;   
  int i,j,k,d; 
  int order=p+1;
  int size=BasisSize(dim,p); 

  double *legendre=malloc(order*sizeof(double));
  double *dxlegendre=malloc(order*sizeof(double));
  double *jacobi, *dxjacobi;
  double factor1, factor2;
  //compute recursive basis indices only during the first call of the file for the same dimension
  if(dimPrev-dim!=0) { 
    if(dimensionCount>0) {
      free(basis);
    }
    basis=malloc((size*dim)*sizeof(int));
    BasisIndices(dim, p, basis);  
    dimensionCount++;
  }
  dimPrev=dim;

  for(k=0; k<m*dim;k++) {
    phiPrime[k]=1.0;
  }

  if(dim==2) {
    int power1, power2, index;
    jacobi=malloc(order*order*sizeof(double));
    dxjacobi=malloc(order*order*sizeof(double));
    
    for(i=0; i<order; i++) {
      JacobiPoly(order, 1.0-2*x[0], 2.0*i+1,0.0, &jacobi[i*order]);
      JacobiPolyPrime(order, 1.0-2*x[0], 2.0*i+1,0.0, &dxjacobi[i*order]);
    }
    legendrePoly(order, (2.0*x[0+1]-x[0])/x[0], legendre, dxlegendre);
  
    for(k=0;k<m;k++) {
      power1=basis[k*dim+0];
      power2=basis[k*dim+0+1]; 
      index=power2+power1*order;
      factor1=1.0;
      factor2=1.0;
      for(j=0; j<power1; j++) {
        factor1=factor1*x[0];
      }
      for(j=0; j<power1-1; j++) {
        factor2=factor2*x[0];
      }
      for(d=0; d<0; d++) {
        phiPrime[k+d*m]=phiPrime[k+d*m]*factor1*jacobi[power2+power1*order]*legendre[power1];  
      }
      phiPrime[k+0*m]=phiPrime[k+0*m]*(-2*x[0+1]/(x[0]*x[0])*dxlegendre[power1]*factor1*jacobi[index] \
      +legendre[power1]*(power1*factor2*jacobi[index]+factor1*dxjacobi[index]*(-2))  ); 
      phiPrime[k+m+0*m]=phiPrime[k+m+0*m]*factor1*jacobi[index]*dxlegendre[power1]*2.0/x[0];   
    }//end loops
    free(jacobi);
    free(dxjacobi);
  }
 
  double *phi=malloc(m*sizeof(double));
  double *jacobiPrev=malloc(order*sizeof(double));
  double *dxjacobiPrev=malloc(order*sizeof(double));
  double *jacobiX=malloc(dim*sizeof(double));
  double *polyX=malloc(dim*sizeof(double));
  double *primeFactor=malloc(dim*sizeof(double));
 
  if(dim>=3) {

    int alpha, xPower,  xPowerPrev, jacobiIndex, jacobiIndexPrev, d1, d2, cur, dimCur, next, prev, indexPrev; 
    double xFactor, xFactorPrev, v1, v2,  phiTemp, polyTemp1, polyTemp2, polyTemp3, polyTemp4,\
    polyFactor1, polyFactor2, polyFactor3, polyFactor4, polyFactor5;
    int *index;

    jacobi=malloc(order*(dim)*(p+1)*sizeof(double));
    dxjacobi=malloc(order*(dim)*(p+1)*sizeof(double));
    
    legendrePoly(order, (2.0*x[dim-1]-x[dim-2])/x[dim-2], &jacobi[0], &dxjacobi[0]); 
    for(d=0; d<dim;d++) {
      polyX[d]=1.0;
    }
    for(d=1; d<dim-1; d++) {
      polyX[d]=x[dim-d-1]/x[dim-d-2];  
    }
    polyX[dim-1]=x[0];
    polyX[0]=0;
     
    for(d=0; d<dim;d++) {
      jacobiX[d]=1.0;
    }

    for(d=1; d<dim;d++) { 
      alpha=0;
      dimCur=d;
      jacobiX[0]=(2.0*x[dim-1]-x[dim-2])/x[dim-2]; 
      if(d>=1&&d<dim-1) {
        jacobiX[dimCur]=jacobiX[dimCur]/x[dim-d-2];
        jacobiX[dimCur]=1.0-2.0*x[dim-d-1]*jacobiX[dimCur]; 
      }
  
      if(d==dim-1) {
        jacobiX[dimCur]=jacobiX[dimCur]*x[0];
        jacobiX[dimCur]=1.0-2.0*jacobiX[dimCur]; 
      }
   
      for(j=0; j<p+1; j++) { 
        alpha=2*j+d;
        if(d==0) {
          alpha=0;
        }
      JacobiPoly(order, jacobiX[dimCur], alpha, 0.0,&jacobi[order*(d)*order+order*j]);
      JacobiPolyPrime(order, jacobiX[dimCur], alpha, 0.0,&dxjacobi[order*(d)*order+order*j]);
      }

    } 
    
    for(k=0; k<m; k++) {
      index=&basis[k*dim+0];

      for(j=0;j<dim;j++) {
        xPower=0;
        if(j>0) {
          for(i=0; i<j; i++) {
            xPower=xPower+*(index+i);  
          }
        }//endif
        jacobiIndex=order*j*order+order*xPower;
        polyTemp1=polyX[j];
        xFactor=1.0;
        for(i=0; i<xPower; i++) {
          xFactor=xFactor*polyTemp1;
        }

        for(d=0; d<dim; d++) {
          primeFactor[d]=jacobi[*(index+j)+jacobiIndex]*xFactor;
        }
        if(j==0) {
          primeFactor[dim-1]=dxjacobi[*(index)+jacobiIndex]*2.0/x[dim-2];
          primeFactor[dim-2]=1.0;
        }
        if(j==1) {
          polyTemp1=x[dim-3];
          polyFactor1=x[dim-3];
          polyTemp2=x[dim-2];
          polyFactor2=1.0;
          for(i=0; i<xPower-1; i++) {
            polyFactor1=polyFactor1*polyTemp1;
            polyFactor2=polyFactor2*polyTemp2;
          }
          primeFactor[dim-2]= jacobi[*(index+1)+jacobiIndex]*xFactor*dxjacobi[*(index)+jacobiIndexPrev] \
          *(-2)*x[dim-1]/(x[dim-2]*x[dim-2])+jacobi[*(index)+jacobiIndexPrev]*(dxjacobi[*(index+1)+jacobiIndex] \
          *(-2)/x[dim-3]*xFactor +jacobi[*(index+1)+jacobiIndex]*xPower*polyFactor2/polyFactor1); 
        }
        if(dim-j-2>=0) {
          primeFactor[dim-j-2]=1.0;
        }
        if(j==dim-1) {
          d2=basis[k*dim+dim-2];
          d1=basis[k*dim+dim-1];
          polyTemp1=x[0+1]/x[0];
          polyTemp2=x[0];
          polyTemp3=x[0+1]; 
          polyFactor1=1.0;
          polyFactor3=1.0;
          polyFactor5=1.0;
          polyFactor2=1.0;
          polyFactor4=1.0;

          for(i=0; i<xPowerPrev; i++) {
            polyFactor1=polyFactor1*polyTemp1;
            polyFactor3=polyFactor3*polyTemp3;
          }
          for(i=0; i<xPower; i++) {
            polyFactor2=polyFactor2*polyTemp2;
          }
          for(i=1; i<xPower; i++) {
            polyFactor4=polyFactor4*polyTemp2;
          }
          for(i=0; i<xPowerPrev+1; i++) {
            polyFactor5=polyFactor5*polyTemp2;
          }

          v1= polyFactor1*jacobi[d2+jacobiIndexPrev]*(dxjacobi[d1+jacobiIndex]*(-2)*polyFactor2 \
          +jacobi[d1+jacobiIndex]*(xPower)*polyFactor4);
          v2= jacobi[d1+jacobiIndex]*polyFactor2*(polyFactor1*dxjacobi[d2+jacobiIndexPrev]*2.0*x[0+1]/(x[0]*x[0]) \
          +jacobi[d2+jacobiIndexPrev]*polyFactor3*(-xPowerPrev)/polyFactor5) ;
          primeFactor[0]= (v1+v2);   
        }
      
        if(j==dim-2) {
          primeFactor[0]=1.0;
        }
   
        if(dim>3&&j!=1&&j!=dim-1&&j!=0) {
          d2=basis[k*dim+0+j-1];
          d1=basis[k*dim+0+j];
          prev=dim-j-2;
          cur=dim-j-1;
          next=dim-j;

          polyTemp1=x[cur]/x[prev];
          polyTemp2=x[next]/x[cur];
          polyTemp3=x[next];
          polyTemp4=x[cur];
          polyFactor1=1.0;
          polyFactor2=1.0;
          polyFactor3=1.0;
          polyFactor4=1.0;
          polyFactor5=1.0;
          for(i=0; i<xPower; i++) {
            polyFactor1=polyFactor1*polyTemp1;
          }
          for(i=0; i<xPower-1; i++) {
            polyFactor5=polyFactor5*polyTemp1;
          }
          for(i=0; i<xPowerPrev+1; i++) {
            polyFactor4=polyFactor4*polyTemp4;
          }
          for(i=0; i<xPowerPrev; i++) {
            polyFactor2=polyFactor2*polyTemp2;
            polyFactor3=polyFactor3*polyTemp3;
          }
          v2= jacobi[d1+jacobiIndex]*polyFactor1*(polyFactor2*dxjacobi[d2+jacobiIndexPrev]*2.0*x[next]/(x[cur]*x[cur]) \
          +jacobi[d2+jacobiIndexPrev]*polyFactor3*(-xPowerPrev)/polyFactor4   ) ;
          v1= polyFactor2*jacobi[d2+jacobiIndexPrev]*(dxjacobi[d1+jacobiIndex]*(-2)/x[prev]*polyFactor1 \
          +jacobi[d1+jacobiIndex]*(xPower)*polyFactor5*1.0/x[prev]); 
          primeFactor[dim-1-j]=v1+v2;
        }//endif  
  
        indexPrev=j;
        xFactorPrev=xFactor; 
        xPowerPrev=xPower;
        jacobiIndexPrev=order*indexPrev*order+order*xPowerPrev;
        
        for(d=0; d<dim; d++) {
          phiPrime[(d+0)*m+k]=phiPrime[(d+0)*m+k]*primeFactor[d];
        } 

      }//end j loop 
    }   

    free(jacobi);
    free(dxjacobi); 
  }//end if dim>=3

  free(primeFactor);
  free(polyX);
  free(jacobiX);
  free(dxjacobiPrev);
  free(jacobiPrev);
  free(phi);
  free(dxlegendre);
  free(legendre);
} /* end PhiPrimeSimplex */

/* PhiCubeSimplex
 * x: point in (dim)-dimensional space
 * phi: array of basis functions evaluated at x
 * params: structure that contains quadrature parameters
 *
 * Generates orthogonal polynomial basis for the unit CUBESIMPLEX
*/
void PhiCubeSimplex(const double *const x, double *phi, const struct quadratureParameters *params){
  int i,j,k, d;  
  int p=params->degreeOfPrecision;
  int m=params->numberOfBasisFunctions;
  int dim1=params->dimension[0];
  int dim2=params->dimension[1];
  int dim=dim1+dim2;   
  int order=p+1;  
  int size=BasisSize(dim,p);  //number of basis functions in dim-dimensions
   
  //compute recursive basis indices only during the first call of the file for the same dimension
  if(dimPrev-dim!=0) { 
    if(dimensionCount>0) {
      free(basis);
    }
    basis=malloc((size*dim)*sizeof(int));
    BasisIndices(dim, p, basis);  
    dimensionCount++;
  }
  dimPrev=dim;

  double* legendre=malloc(order*dim1*sizeof(double));
  double* dxlegendre=malloc(order*dim1*sizeof(double));
  for(k=0; k<m; k++) {
    phi[k]=1.0;
  }
  for(d=0;d<dim1;d++) {   
    legendrePoly(order, 2*x[d]-1, &legendre[d*order], &dxlegendre[d*order]); 
  }
  for(k=0;k<m;k++) {
    for(d=0;d<dim1;d++) {
      phi[k]=phi[k]*legendre[basis[k*dim+d]+order*d];   
    }    
  }   
  free(legendre);
  free(dxlegendre); 

  double phiSimplex[m];
  PhiSimplexPolyhedralTwo(x, phiSimplex, params);
  for(k=0; k<m; k++) {
    phi[k]=phi[k]*phiSimplex[k];
  }    
} /*end PhiCubeSimplex */      


/* PhiPrimeCubeSimplex
 * x: point in (dim)-dimensional space
 * phiPrime: array of derivatives of basis functions evaluated at x
 * params: structure that contains quadrature parameters
 *
 * Generates derivatives of orthogonal polynomial basis for the unit CUBESIMPLEX
*/
void PhiPrimeCubeSimplex(const double *const x,double *phiPrime, const struct quadratureParameters *params) {
  int p=params->degreeOfPrecision;
  int m=params->numberOfBasisFunctions;
  int dim1=params->dimension[0];
  int dim2=params->dimension[1];
  int dim=dim1+dim2;   
  int i,j,k,d; 
  int order=p+1;
  int size=BasisSize(dim,p); 
  double legendre[order*dim];
  double dxlegendre[order*dim]; 
  double phiPrimeSimplex[m*dim];
  //compute recursive basis indices only during the first call of the file for the same dimension
  if(dimPrev-dim!=0) { 
    if(dimensionCount>0) {
      free(basis);
    }
    basis=malloc((size*dim)*sizeof(int));
    BasisIndices(dim, p, basis);  
    dimensionCount++;
  }
  dimPrev=dim;

  for(k=0; k<m*dim;k++) {
    phiPrime[k]=1.0;
  }
  for(d=0;d<dim1;d++) {   
    legendrePoly(order, 2*x[d]-1, &legendre[d*order], &dxlegendre[d*order]);  
  }

  for(d=0; d<dim; d++) {
    for(k=0;k<m;k++) {
  
      for(j=0;j<dim1;j++) {
        if(j!=d) {
          phiPrime[k+d*m]=phiPrime[k+d*m]*legendre[basis[k*dim+j]+order*j];  
        }
        if(j==d) {
          phiPrime[k+d*m]=2*phiPrime[k+d*m]*dxlegendre[basis[k*dim+j]+order*j];   
        }
      } 
    }
  }
  PhiPrimeSimplexPolyhedralTwo(x, phiPrimeSimplex, params);
  for(k=0; k<m*dim; k++) {
    phiPrime[k]=phiPrime[k]*phiPrimeSimplex[k];
  }  

}/* end PhiPrimeCubeSimplex */

/* PhiSimplexSimplex
 * x: point in (dim)-dimensional space
 * phi: array of basis functions evaluated at x
 * params: structure that contains quadrature parameters
 *
 * Generates orthogonal polynomial basis for the unit SIMPLEXSIMPLEX
*/
void PhiSimplexSimplex(const double *const x, double *phi, const struct quadratureParameters *params){
  int k;  
  int p=params->degreeOfPrecision;
  int m=params->numberOfBasisFunctions;
  int dim=params->totalDimension;  
  int size=BasisSize(dim,p);   //number of basis functions in dim-dimensions
  double phiSimplex[m];
  //compute recursive basis indices only during the first call of the file for the same dimension
  if(dimPrev-dim!=0) { 
    if(dimensionCount>0) {
      free(basis);
    }
    basis=malloc((size*dim)*sizeof(int));
    BasisIndices(dim, p, basis);  
    dimensionCount++;
  }
  dimPrev=dim;

  for(k=0; k<m; k++) {
    phi[k]=1.0;
  }

  PhiSimplexPolyhedralOne(x, phiSimplex, params);
  for(k=0; k<m; k++) {
    phi[k]=phi[k]*phiSimplex[k];
  }    
  PhiSimplexPolyhedralTwo(x, phiSimplex, params);
  for(k=0; k<m; k++) {
    phi[k]=phi[k]*phiSimplex[k];
  }    

} /* end PhiSimplexSimplex */

/* PhiPrimeSimplexSimplex
 * x: point in (dim)-dimensional space
 * phiPrime: array of derivatives of basis functions evaluated at x
 * params: structure that contains quadrature parameters
 *
 * Generates derivatives of orthogonal polynomial basis for the unit SIMPLEXSIMPLEX
*/
void PhiPrimeSimplexSimplex(const double *const x, double *phiPrime, const struct quadratureParameters *params){
  int k; 
  int p=params->degreeOfPrecision;
  int m=params->numberOfBasisFunctions;
  int dim=params->totalDimension;
  int size=BasisSize(dim,p); 
  double phiPrimeSimplex[m*dim];

  //compute recursive basis indices only during the first call of the file for the same dimension
  if(dimPrev-dim!=0) { 
    if(dimensionCount>0) {
      free(basis);
    }
    basis=malloc((size*dim)*sizeof(int));
    BasisIndices(dim, p, basis);  
    dimensionCount++;
  }
  dimPrev=dim;

  for(k=0; k<m*dim;k++) {
    phiPrime[k]=1.0;
  }
  PhiPrimeSimplexPolyhedralOne(x, phiPrimeSimplex, params);
  for(k=0; k<m*dim; k++) {
    phiPrime[k]=phiPrime[k]*phiPrimeSimplex[k];
  }  
  PhiPrimeSimplexPolyhedralTwo(x, phiPrimeSimplex, params);
  for(k=0; k<m*dim; k++) {
    phiPrime[k]=phiPrime[k]*phiPrimeSimplex[k];
  } 

} /* end PhiPrimeSimplexSimplex */

/* PhiCubeSimplexSimplex
 * x: point in (dim)-dimensional space
 * phi: array of basis functions evaluated at x
 * params: structure that contains quadrature parameters
 *
 * Generates orthogonal polynomial basis for the unit CUBESIMPLEXSIMPLEX
*/
void PhiCubeSimplexSimplex(const double *const x, double *phi, const struct quadratureParameters *params) {  

  int k,d;  
  int p=params->degreeOfPrecision;
  int order=p+1; 
  int m=params->numberOfBasisFunctions;
  int dimSimplex1=params->dimension[0];
  int dimSimplex2=params->dimension[1];
  int dimCube=params->dimension[2];
  int twoDims=dimSimplex1+dimSimplex2;
  int dim=params->totalDimension;  
  struct quadratureParameters *newParams =(struct quadratureParameters *)malloc(sizeof(struct quadratureParameters));
  newParams->dimension=malloc(3*sizeof(int));
  SetParams(dimSimplex1, dimSimplex2 , dimCube, p,  newParams);
  int size=BasisSize(dim,p);  //number of basis functions in dim-dimensions
  double phiSimplex[m];

   //compute recursive basis indices only during the first call of the file for the same dimension
  if(dimPrev-dim!=0) { 
    if(dimensionCount>0) {
      free(basis);
    }
    basis=malloc((size*dim)*sizeof(int));
    BasisIndices(dim, p, basis);  
    dimensionCount++;
  }
  dimPrev=dim;

  for(k=0; k<m; k++) {
    phi[k]=1.0;
  }
  PhiSimplexPolyhedralOne(x, phiSimplex, newParams);
  for(k=0; k<m; k++) {
    phi[k]=phi[k]*phiSimplex[k];
  }    
  PhiSimplexPolyhedralTwo(x, phiSimplex, newParams);
  for(k=0; k<m; k++) {
    phi[k]=phi[k]*phiSimplex[k];
  }

  double* legendre=malloc(order*dimCube*sizeof(double));
  double* dxlegendre=malloc(order*dimCube*sizeof(double));
  for(d=0;d<dimCube;d++) {   
    legendrePoly(order, 2*x[twoDims+d]-1, &legendre[(d)*order], &dxlegendre[d*order]); 
  }
  for(k=0;k<m;k++) {
    for(d=0;d<dimCube;d++) {
      phi[k]=phi[k]*legendre[basis[k*dim+twoDims+d]+order*d];   
    }    
  }   
  free(legendre);
  free(dxlegendre); 
  free(newParams->dimension);
  free(newParams); 
} /* end PhiCubeSimplexSimplex */


/* PhiPrimeCubeSimplexSimplex
 * x: point in (dim)-dimensional space
 * phiPrime: array of derivatives of basis functions evaluated at x
 * params: structure that contains quadrature parameters
 *
 * Generates derivatives of orthogonal polynomial basis for the unit CUBESIMPLEXSIMPLEX
*/
void PhiPrimeCubeSimplexSimplex(const double *const x, double *phiPrime, const struct quadratureParameters *params) {

  int j,k,d;  
  int p=params->degreeOfPrecision;
  int order=p+1;  
  int m=params->numberOfBasisFunctions;
  int dimSimplex1=params->dimension[0];
  int dimSimplex2=params->dimension[1];
  int dimCube=params->dimension[2];
  int twoDims=dimSimplex1+dimSimplex2;
  int dim=params->totalDimension;  
  int size=BasisSize(dim,p);
  double phiPrimeSimplex[m*dim];

  //compute recursive basis indices only during the first call of the file for the same dimension
  if(dimPrev-dim!=0) { 
    if(dimensionCount>0) {
      free(basis);
    }
    basis=malloc((size*dim)*sizeof(int));
    BasisIndices(dim, p, basis);  
    dimensionCount++;
  }
  dimPrev=dim;

  for(k=0; k<m*dim;k++) {
    phiPrime[k]=1.0;
  }
  PhiPrimeSimplexPolyhedralOne(x, phiPrimeSimplex, params);
  for(k=0; k<m*dim; k++) {
    phiPrime[k]=phiPrime[k]*phiPrimeSimplex[k];
  }  
  PhiPrimeSimplexPolyhedralTwo(x, phiPrimeSimplex, params);
  for(k=0; k<m*dim; k++) {
    phiPrime[k]=phiPrime[k]*phiPrimeSimplex[k];
  }

  double* legendre=malloc(order*dimCube*sizeof(double));
  double* dxlegendre=malloc(order*dimCube*sizeof(double));
  for(d=0;d<dimCube;d++) {   
    legendrePoly(order, 2*x[twoDims+d]-1, &legendre[d*order], &dxlegendre[d*order]);  
  }
  for(d=0; d<dim; d++) {
    for(k=0;k<m;k++) {
  
      for(j=twoDims;j<dim;j++) {
        if(j!=d) {
          phiPrime[k+d*m]=phiPrime[k+d*m]*legendre[basis[k*dim+j]+order*(j-twoDims)];  
        }
        if(j==d) {
          phiPrime[k+d*m]=2*phiPrime[k+d*m]*dxlegendre[basis[k*dim+j]+order*(j-twoDims)];   
        }
      } 
    }
  }

  free(legendre);
  free(dxlegendre);
} /* PhiPrimeCubeSimplexSimplex */

/* PhiSimplexPolyhedralOne
 * x: point in (dim)-dimensional space
 * phi: array of derivatives of basis functions evaluated at x
 * params: structure that contains quadrature parameters
 *
 * Generates orthogonal polynomial basis for the first [0, dim1-1] coordinates 
 * over simplex in dimTotal-dimensional space. 
*/
void PhiSimplexPolyhedralOne(const double *const x, double *phi, const struct quadratureParameters *params){
  int i,j,k, d;  
  int p=params->degreeOfPrecision;
  int m=params->numberOfBasisFunctions;
  int dim1=params->dimension[0];
  int dim2=params->dimension[1];
  int dimTwo=dim1+dim2;
  int dim=params->totalDimension;   
  int order=p+1;  
  int size=BasisSize(dim,p);  
   
  for(k=0; k<m; k++) {
    phi[k]=1.0;
  }

  int dimCur=0,  xPower=0;
  double xFactor, xTemp=0, phiTemp=0, alpha=0;
  int *index;
  double *jacobiX, *polyX;
  double jacobi[order*order*(dim1-1)];
  double *legendre=malloc(order*sizeof(double));
  double *dxlegendre=malloc(order*sizeof(double)); 

  if(dim1==2) {
    index=&basis[0];
    for(j=0; j<order; j++) {
      JacobiPoly(order, 1.0-2*x[0], 2.0*j+1,0.0, &jacobi[j*order]);
    }
    int r,n; 
    legendrePoly(order, (2.0*x[1]-x[0])/x[0], legendre, dxlegendre);
     
    for(k=0;k<m;k++) { 
      r=*index;
      n=*(index+1);
      xFactor=1.0;
      for(i=0; i<r; i++) {
        xFactor=xFactor*x[0];
      }
      phi[k]=phi[k]*xFactor*jacobi[r*order+n]*legendre[r];  
      index=(index+dim);
    }
  } 
 
  if(dim1>=3) {
    int *index;
    jacobiX=malloc((dim1-1)*sizeof(double));  
    polyX=malloc((dim1-1)*sizeof(double));
    legendrePoly(order, (2.0*x[dim1-1]-x[dim1-2])/x[dim1-2], legendre, dxlegendre);
   
    for(d=0; d<dim1-1;d++) {
      polyX[d]=1.0;
    } 
    for(d=0; d<dim1-2; d++) {
      polyX[d]=x[dim1-d-2]/x[dim1-d-3];  
    }
    polyX[dim1-2]=x[0];
    for(d=0; d<dim1-1;d++) {
      jacobiX[d]=1.0;
    }

    for(d=1; d<dim1;d++) {
      alpha=0.0;
      dimCur=d-1;
      if(d>=1&&d<dim1-1) {
        jacobiX[dimCur]=jacobiX[dimCur]/x[dim1-d-2];
        jacobiX[dimCur]=1.0-2.0*x[dim1-d-1]*jacobiX[dimCur]; 
      }
      if(d==dim1-1) {
        jacobiX[dimCur]=jacobiX[dimCur]*x[0];
        jacobiX[dimCur]=1.0-2.0*jacobiX[dimCur]; 
      }
      for(j=0; j<p+1; j++) {
        alpha=2*j+d;
        JacobiPoly(order, jacobiX[dimCur], alpha, 0.0,&jacobi[order*(d-1)*p+order*j]);
      }
    }

    for(k=0; k<m; k++) {
      index=&basis[k*dim];
      phiTemp=1.0;
      for(d=1;d<dim1;d++) { 
        xPower=0;
        for(i=0; i<d; i++) {  
          xPower=xPower+*(index+i);
        }  
        xTemp=polyX[d-1];
        xFactor=1.0;
        for(i=0; i<xPower; i++) {
          xFactor=xFactor*xTemp;
        }
        phiTemp=phiTemp*jacobi[*(index+d)+order*(d-1)*p+order*xPower]*xFactor; 
      } 
      phi[k]=phi[k]*phiTemp*legendre[*index]; 
    } 

    free(jacobiX);
    free(polyX);
  }

  free(legendre);
  free(dxlegendre);
} /* end PhiSimplexPolyhedralOne */

/* PhiSimplexPolyhedralTwo
 * x: point in (dim)-dimensional space
 * phi: array of of basis functions evaluated at x
 * params: structure that contains quadrature parameters
 *
 * Generates orthogonal polynomial basis for coordinates between [dim1, dim2-1]
 * over simplex in dimTotal-dimensional space. 
*/
void PhiSimplexPolyhedralTwo(const double *const x, double *phi, const struct quadratureParameters *params){
  int i,j,k, d;  
  int p=params->degreeOfPrecision;
  int m=params->numberOfBasisFunctions;
  int dim1=params->dimension[0];
  int dim2=params->dimension[1];
  int dimTwo=dim1+dim2;
  int dim=params->totalDimension;   
  int order=p+1;  
  int size=BasisSize(dim,p);  
   
  for(k=0; k<m; k++) {
    phi[k]=1.0;
  }
 
  int dimCur=0,  xPower=0;
  double xFactor, xTemp=0, phiTemp=0, alpha=0;
  int *index;
  double *jacobiX, *polyX;
  double jacobi[order*order*(dim2-1)];
  double *legendre=malloc(order*sizeof(double));
  double *dxlegendre=malloc(order*sizeof(double)); 

  if(dim2==2) {
    index=&basis[dim1];
    for(j=0; j<order; j++) {
      JacobiPoly(order, 1.0-2*x[dim1], 2.0*j+1,0.0, &jacobi[j*order]);
    }
    int r,n; 
    legendrePoly(order, (2.0*x[dim1+1]-x[dim1])/x[dim1], legendre, dxlegendre);
     
    for(k=0;k<m;k++) { 
      r=*index;
      n=*(index+1);
      xFactor=1.0;
      for(i=0; i<r; i++) {
        xFactor=xFactor*x[dim1];
      }
      phi[k]=phi[k]*xFactor*jacobi[r*order+n]*legendre[r];  
      index=(index+dim);
    }
 
  }  

  if(dim2>=3) {
    int *index;
    jacobiX=malloc((dim2-1)*sizeof(double));  
    polyX=malloc((dim2-1)*sizeof(double));
    legendrePoly(order, (2.0*x[dimTwo-1]-x[dimTwo-2])/x[dimTwo-2], legendre, dxlegendre);
    
    for(d=0; d<dim2-1;d++) {
      polyX[d]=1.0;
    } 
    for(d=0; d<dim2-2; d++) {
      polyX[d]=x[dimTwo-d-2]/x[dimTwo-d-3];  
    }
    polyX[dim2-2]=x[dim1];
    for(d=0; d<dim2-1;d++) {
      jacobiX[d]=1.0;
    }

    for(d=1; d<dim2;d++) {
      alpha=0.0;
      dimCur=d-1;
      if(d>=1&&d<dim2-1) {
        jacobiX[dimCur]=jacobiX[dimCur]/x[dimTwo-d-2];
        jacobiX[dimCur]=1.0-2.0*x[dimTwo-d-1]*jacobiX[dimCur]; 
      }
      if(d==dim2-1) {
        jacobiX[dimCur]=jacobiX[dimCur]*x[dim1];
        jacobiX[dimCur]=1.0-2.0*jacobiX[dimCur]; 
      }
      for(j=0; j<p+1; j++) {
        alpha=2*j+d;
        JacobiPoly(order, jacobiX[dimCur], alpha, 0.0,&jacobi[order*(d-1)*p+order*j]);
      }
    }
 
    for(k=0; k<m; k++) {
      index=&basis[k*dim+dim1];
      phiTemp=1.0;
      for(d=1;d<dim2;d++) { 
        xPower=0;
        for(i=0; i<d; i++) {  
          xPower=xPower+*(index+i);
        } 
        xTemp=polyX[d-1];
        xFactor=1.0;
        for(i=0; i<xPower; i++) {
          xFactor=xFactor*xTemp;
        }
        phiTemp=phiTemp*jacobi[*(index+d)+order*(d-1)*p+order*xPower]*xFactor; 
      } 
      phi[k]=phi[k]*phiTemp*legendre[*index]; 
    } 

    free(jacobiX);
    free(polyX);
  }

  free(legendre);
  free(dxlegendre);
} /* end PhiSimplexPolyhedralTwo */

/* PhiPrimeSimplexPolyhedralOne
 * x: point in (dim)-dimensional space
 * phiPrime: array of derivatives of basis functions evaluated at x
 * params: structure that contains quadrature parameters
 *
 * Generates derivatives of orthogonal polynomial basis for the first [0, dim1-1] coordinates 
 * over simplex in dimTotal-dimensional space. 
*/
void PhiPrimeSimplexPolyhedralOne(const double *const x, double *phiPrime, const struct quadratureParameters *params){
  int p=params->degreeOfPrecision;
  int m=params->numberOfBasisFunctions;
  int dim1=params->dimension[0];
  int dim2=params->dimension[1];
  int dimTwo=dim1+dim2;
  int dim=params->totalDimension;   
  int i,j,k,d; 
  int order=p+1;
  int size=BasisSize(dim,p); 

  double factor1, factor2;
  double legendre[order];
  double dxlegendre[order];
  double *jacobi, *dxjacobi;

  for(k=0; k<m*dim;k++) {
    phiPrime[k]=1.0;
  }

  if(dim1==2) {
    int power1, power2, index;
    jacobi=malloc(order*order*sizeof(double));
    dxjacobi=malloc(order*order*sizeof(double));
  
    for(i=0; i<order; i++) {
      JacobiPoly(order, 1.0-2*x[0], 2.0*i+1,0.0, &jacobi[i*order]);
      JacobiPolyPrime(order, 1.0-2*x[0], 2.0*i+1,0.0, &dxjacobi[i*order]);
    }
    legendrePoly(order, (2.0*x[1]-x[0])/x[0], legendre, dxlegendre);
  
    for(k=0;k<m;k++) {
      power1=basis[k*dim];
      power2=basis[k*dim+1]; 
      index=power2+power1*order;
      factor1=1.0;
      factor2=1.0;
      for(j=0; j<power1; j++) {
        factor1=factor1*x[0];
      }
      for(j=0; j<power1-1; j++) {
        factor2=factor2*x[0];
      }
      for(d=2; d<dimTwo; d++) {
        phiPrime[k+d*m]=phiPrime[k+d*m]*factor1*jacobi[power2+power1*order]*legendre[power1];  
      }
      phiPrime[k]=phiPrime[k]*(-2*x[0+1]/(x[0]*x[0])*dxlegendre[power1]*factor1*jacobi[index] \
      +legendre[power1]*(power1*factor2*jacobi[index]+factor1*dxjacobi[index]*(-2))  ); 
      phiPrime[k+m]=phiPrime[k+m]*factor1*jacobi[index]*dxlegendre[power1]*2.0/x[0]; 
   
      if(dim>dimTwo){
        for(d=dimTwo; d<dim; d++) {
          phiPrime[k+d*m]=phiPrime[k+d*m]*factor1*jacobi[power2+power1*order]*legendre[power1];  
        }
      }
    }//end loops
    free(jacobi);
    free(dxjacobi);
  }

  double *phi=malloc(m*sizeof(double));
  double *jacobiPrev=malloc(order*sizeof(double));
  double *dxjacobiPrev=malloc(order*sizeof(double));
  double *jacobiX=malloc(dim1*sizeof(double));
  double *polyX=malloc(dim1*sizeof(double));
  double *primeFactor=malloc(dim1*sizeof(double));
  if(dim1>=3) {
    int alpha, xPower,  xPowerPrev, jacobiIndex, jacobiIndexPrev, d1, d2, cur, dimCur, next, prev, indexPrev; 
    double xFactor, xFactorPrev, v1, v2,  phiTemp, polyTemp1, polyTemp2, polyTemp3, polyTemp4,\
    polyFactor1, polyFactor2, polyFactor3, polyFactor4, polyFactor5;
    int *index;
    PhiSimplexPolyhedralOne(x,phi, params); 
    jacobi=malloc(order*(dim1)*(p+1)*sizeof(double));
    dxjacobi=malloc(order*(dim1)*(p+1)*sizeof(double));
  
    legendrePoly(order, (2.0*x[dim1-1]-x[dim1-2])/x[dim1-2], &jacobi[0], &dxjacobi[0]); 
    for(d=0; d<dim1;d++) {
      polyX[d]=1.0;
    }
    for(d=1; d<dim1-1; d++) {
      polyX[d]=x[dim1-d-1]/x[dim1-d-2];  
    }
    polyX[dim1-1]=x[0];
    polyX[0]=0;
   
    for(d=0; d<dim1;d++) {
      jacobiX[d]=1.0;
    }

    for(d=1; d<dim1;d++) { 
      alpha=0;
      dimCur=d;
      jacobiX[0]=(2.0*x[dim1-1]-x[dim1-2])/x[dim1-2]; 
      if(d>=1&&d<dim1-1) {
        jacobiX[dimCur]=jacobiX[dimCur]/x[dim1-d-2];
        jacobiX[dimCur]=1.0-2.0*x[dim1-d-1]*jacobiX[dimCur]; 
      }
    
      if(d==dim1-1) {
        jacobiX[dimCur]=jacobiX[dimCur]*x[0];
        jacobiX[dimCur]=1.0-2.0*jacobiX[dimCur]; 
      }
   
      for(j=0; j<p+1; j++) { 
        alpha=2*j+d;
        if(d==0) {
          alpha=0;
        }
      JacobiPoly(order, jacobiX[dimCur], alpha, 0.0,&jacobi[order*(d)*order+order*j]);
      JacobiPolyPrime(order, jacobiX[dimCur], alpha, 0.0,&dxjacobi[order*(d)*order+order*j]);
      }
    } 
    
    for(k=0; k<m; k++) {
      index=&basis[k*dim];

      for(j=0;j<dim1;j++) {
        xPower=0;
        if(j>0) {
          for(i=0; i<j; i++) {
            xPower=xPower+*(index+i);  
          }
        }
      jacobiIndex=order*j*order+order*xPower;
      polyTemp1=polyX[j];
      xFactor=1.0;
      for(i=0; i<xPower; i++) {
        xFactor=xFactor*polyTemp1;
      }

      for(d=0; d<dim1; d++) {
        primeFactor[d]=jacobi[*(index+j)+jacobiIndex]*xFactor;
      }
      if(j==0) {
        primeFactor[dim1-1]=dxjacobi[*(index)+jacobiIndex]*2.0/x[dim1-2];
       primeFactor[dim1-2]=1.0;
      }
      if(j==1) {
        polyTemp1=x[dim1-3];
        polyFactor1=x[dim1-3];
        polyTemp2=x[dim1-2];
        polyFactor2=1.0;
        for(i=0; i<xPower-1; i++) {
         polyFactor1=polyFactor1*polyTemp1;
          polyFactor2=polyFactor2*polyTemp2;
        }
        primeFactor[dim1-2]= jacobi[*(index+1)+jacobiIndex]*xFactor*dxjacobi[*(index)+jacobiIndexPrev] \
        *(-2)*x[dim1-1]/(x[dim1-2]*x[dim1-2])+jacobi[*(index)+jacobiIndexPrev]*(dxjacobi[*(index+1)+jacobiIndex] \
        *(-2)/x[dim1-3]*xFactor +jacobi[*(index+1)+jacobiIndex]*xPower*polyFactor2/polyFactor1); 
      }
      if(dim1-j-2>=0) {
        primeFactor[dim1-j-2]=1.0;
      }
      if(j==dim1-1) {
        d2=basis[k*dim+dim1-2];
        d1=basis[k*dim+dim1-1];
        polyTemp1=x[0+1]/x[0];
        polyTemp2=x[0];
        polyTemp3=x[0+1]; 
        polyFactor1=1.0;
        polyFactor3=1.0;
        polyFactor5=1.0;
        polyFactor2=1.0;
        polyFactor4=1.0;
      
       for(i=0; i<xPowerPrev; i++) {
         polyFactor1=polyFactor1*polyTemp1;
         polyFactor3=polyFactor3*polyTemp3;
        }
       for(i=0; i<xPower; i++) {
         polyFactor2=polyFactor2*polyTemp2;
       }
       for(i=1; i<xPower; i++) {
         polyFactor4=polyFactor4*polyTemp2;
       }
       for(i=0; i<xPowerPrev+1; i++) {
         polyFactor5=polyFactor5*polyTemp2;
       }

       v1= polyFactor1*jacobi[d2+jacobiIndexPrev]*(dxjacobi[d1+jacobiIndex]*(-2)*polyFactor2 \
       +jacobi[d1+jacobiIndex]*(xPower)*polyFactor4);
       v2= jacobi[d1+jacobiIndex]*polyFactor2*(polyFactor1*dxjacobi[d2+jacobiIndexPrev]*2.0*x[0+1]/(x[0]*x[0]) \
       +jacobi[d2+jacobiIndexPrev]*polyFactor3*(-xPowerPrev)/polyFactor5) ;
       primeFactor[0]= (v1+v2);   
      }
      
      if(j==dim1-2) {
        primeFactor[0]=1.0;
      }
      
      if(dim1>3&&j!=1&&j!=dim1-1&&j!=0) {
        d2=basis[k*dim+j-1];
        d1=basis[k*dim+j];
        prev=dim1-j-2;
        cur=dim1-j-1;
        next=dim1-j;

        polyTemp1=x[cur]/x[prev];
        polyTemp2=x[next]/x[cur];
        polyTemp3=x[next];
        polyTemp4=x[cur];
        polyFactor1=1.0;
        polyFactor2=1.0;
        polyFactor3=1.0;
        polyFactor4=1.0;
        polyFactor5=1.0;
        for(i=0; i<xPower; i++) {
         polyFactor1=polyFactor1*polyTemp1;
       }
       for(i=0; i<xPower-1; i++) {
         polyFactor5=polyFactor5*polyTemp1;
       }
       for(i=0; i<xPowerPrev+1; i++) {
         polyFactor4=polyFactor4*polyTemp4;
       }
       for(i=0; i<xPowerPrev; i++) {
         polyFactor2=polyFactor2*polyTemp2;
         polyFactor3=polyFactor3*polyTemp3;
       }
      
        v2= jacobi[d1+jacobiIndex]*polyFactor1*(polyFactor2*dxjacobi[d2+jacobiIndexPrev]*2.0*x[next]/(x[cur]*x[cur])+ \
       jacobi[d2+jacobiIndexPrev]*polyFactor3*(-xPowerPrev)/polyFactor4   ) ;
       v1= polyFactor2*jacobi[d2+jacobiIndexPrev]*(dxjacobi[d1+jacobiIndex]*(-2)/x[prev]*polyFactor1+ \
       jacobi[d1+jacobiIndex]*(xPower)*polyFactor5*1.0/x[prev]); 
       primeFactor[dim1-1-j]=v1+v2;
      }//endif  
    
      indexPrev=j;
      xFactorPrev=xFactor; 
      xPowerPrev=xPower;
      jacobiIndexPrev=order*indexPrev*order+order*xPowerPrev;
    
      for(d=0; d<dim1; d++) {
        phiPrime[d*m+k]=phiPrime[d*m+k]*primeFactor[d];
      }  
    }//end j loop 
  }   
 
  for(d=dim1; d<dim; d++) {
    for(k=0;k<m;k++) {  
      phiPrime[k+d*m]=phiPrime[k+d*m]*phi[k];  
    }
  }  
  free(dxjacobi);
  free(jacobi);
  }//end if dim>=3
  free(primeFactor);
  free(polyX);
  free(jacobiX);
  free(dxjacobiPrev);
  free(jacobiPrev);
  free(phi);
}/* end PhiPrimeSimplexPolyhedralOne */

/* PhiPrimeSimplexPolyhedralTwo
 * x: point in (dim)-dimensional space
 * phiPrime: array of derivatives of basis functions evaluated at x
 * params: structure that contains quadrature parameters
 *
 * Generates derivatives of orthogonal polynomial basis for coordinates between [dim1, dim2-1]
 * over simplex in dimTotal-dimensional space. 
*/
void PhiPrimeSimplexPolyhedralTwo(const double *const x, double *phiPrime, const struct quadratureParameters *params){
  int p=params->degreeOfPrecision;
  int m=params->numberOfBasisFunctions;
  int dim1=params->dimension[0];
  int dim2=params->dimension[1];
  int dimTwo=dim1+dim2;
  int dim=params->totalDimension;   
  int i,j,k,d; 
  int order=p+1;
  int size=BasisSize(dim,p); 

  double factor1, factor2;
  double legendre[order];
  double dxlegendre[order];
  double *jacobi, *dxjacobi;

  for(k=0; k<m*dim;k++) {
    phiPrime[k]=1.0;
  }

  if(dim2==2) {
    int power1, power2, index;
    jacobi=malloc(order*order*sizeof(double));
    dxjacobi=malloc(order*order*sizeof(double));
   
    for(i=0; i<order; i++) {
      JacobiPoly(order, 1.0-2*x[dim1], 2.0*i+1,0.0, &jacobi[i*order]);
      JacobiPolyPrime(order, 1.0-2*x[dim1], 2.0*i+1,0.0, &dxjacobi[i*order]);
    }
    legendrePoly(order, (2.0*x[dim1+1]-x[dim1])/x[dim1], legendre, dxlegendre);
   
    for(k=0;k<m;k++) {
      power1=basis[k*dim+dim1];
      power2=basis[k*dim+dim1+1]; 
      index=power2+power1*order;
      factor1=1.0;
      factor2=1.0;
      for(j=0; j<power1; j++) {
        factor1=factor1*x[dim1];
      }
      for(j=0; j<power1-1; j++) {
        factor2=factor2*x[dim1];
      }
      for(d=0; d<dim1; d++) {
        phiPrime[k+d*m]=phiPrime[k+d*m]*factor1*jacobi[power2+power1*order]*legendre[power1];  
      }
      phiPrime[k+dim1*m]=phiPrime[k+dim1*m]*(-2*x[dim1+1]/(x[dim1]*x[dim1])*dxlegendre[power1]*factor1*jacobi[index] \
      +legendre[power1]*(power1*factor2*jacobi[index]+factor1*dxjacobi[index]*(-2))  ); 
      phiPrime[k+m+dim1*m]=phiPrime[k+m+dim1*m]*factor1*jacobi[index]*dxlegendre[power1]*2.0/x[dim1];   
      
      if(dim>dimTwo){
        for(d=dimTwo; d<dim; d++) {
         phiPrime[k+d*m]=phiPrime[k+d*m]*factor1*jacobi[power2+power1*order]*legendre[power1];  
        }
      }
    }//end loops
    free(jacobi);
    free(dxjacobi);
  }
  
  double *phi=malloc(m*sizeof(double));
  double *jacobiPrev=malloc(order*sizeof(double));
  double *dxjacobiPrev=malloc(order*sizeof(double));
  double *jacobiX=malloc(dim2*sizeof(double));
  double *polyX=malloc(dim2*sizeof(double));
  double *primeFactor=malloc(dim2*sizeof(double));  
  if(dim2>=3) {

    int alpha, xPower,  xPowerPrev, jacobiIndex, jacobiIndexPrev, d1, d2, cur, dimCur, next, prev, indexPrev; 
    double xFactor, xFactorPrev, v1, v2,  phiTemp, polyTemp1, polyTemp2, polyTemp3, polyTemp4,\
    polyFactor1, polyFactor2, polyFactor3, polyFactor4, polyFactor5;
    int *index;
    PhiSimplexPolyhedralTwo(x,phi, params); 
    jacobi=malloc(order*(dim2)*(p+1)*sizeof(double));
    dxjacobi=malloc(order*(dim2)*(p+1)*sizeof(double));
    
    legendrePoly(order, (2.0*x[dimTwo-1]-x[dimTwo-2])/x[dimTwo-2], &jacobi[0], &dxjacobi[0]); 
    for(d=0; d<dim2;d++) {
      polyX[d]=1.0;
    }
    for(d=1; d<dim2-1; d++) {
      polyX[d]=x[dimTwo-d-1]/x[dimTwo-d-2];  
    }
    polyX[dim2-1]=x[dim1];
    polyX[0]=0;
     
    for(d=0; d<dim2;d++) {
      jacobiX[d]=1.0;
    }

    for(d=1; d<dim2;d++) { 
      alpha=0;
      dimCur=d;
      jacobiX[0]=(2.0*x[dimTwo-1]-x[dimTwo-2])/x[dimTwo-2]; 
      if(d>=1&&d<dim2-1) {
        jacobiX[dimCur]=jacobiX[dimCur]/x[dimTwo-d-2];
        jacobiX[dimCur]=1.0-2.0*x[dimTwo-d-1]*jacobiX[dimCur]; 
      }
    
      if(d==dim2-1) {
        jacobiX[dimCur]=jacobiX[dimCur]*x[dim1];
        jacobiX[dimCur]=1.0-2.0*jacobiX[dimCur]; 
      }
     
      for(j=0; j<p+1; j++) { 
        alpha=2*j+d;
        if(d==0) {
          alpha=0;
        }
        JacobiPoly(order, jacobiX[dimCur], alpha, 0.0,&jacobi[order*(d)*order+order*j]);
        JacobiPolyPrime(order, jacobiX[dimCur], alpha, 0.0,&dxjacobi[order*(d)*order+order*j]);
      }
    } 
    
    for(k=0; k<m; k++) {
      index=&basis[k*dim+dim1];

      for(j=0;j<dim2;j++) {
        xPower=0;
        if(j>0) {
          for(i=0; i<j; i++) {
            xPower=xPower+*(index+i);  
          }
        }
        jacobiIndex=order*j*order+order*xPower;
        polyTemp1=polyX[j];
        xFactor=1.0;
        for(i=0; i<xPower; i++) {
          xFactor=xFactor*polyTemp1;
        }

        for(d=0; d<dim2; d++) {
          primeFactor[d]=jacobi[*(index+j)+jacobiIndex]*xFactor;
        }
        if(j==0) {
          primeFactor[dim2-1]=dxjacobi[*(index)+jacobiIndex]*2.0/x[dimTwo-2];
          primeFactor[dim2-2]=1.0;
        }
        if(j==1) {
          polyTemp1=x[dimTwo-3];
          polyFactor1=x[dimTwo-3];
          polyTemp2=x[dimTwo-2];
          polyFactor2=1.0;
          for(i=0; i<xPower-1; i++) {
            polyFactor1=polyFactor1*polyTemp1;
            polyFactor2=polyFactor2*polyTemp2;
          }
          primeFactor[dim2-2]= jacobi[*(index+1)+jacobiIndex]*xFactor*dxjacobi[*(index)+jacobiIndexPrev] \
          *(-2)*x[dimTwo-1]/(x[dimTwo-2]*x[dimTwo-2])+jacobi[*(index)+jacobiIndexPrev]*(dxjacobi[*(index+1)+jacobiIndex] \
          *(-2)/x[dimTwo-3]*xFactor +jacobi[*(index+1)+jacobiIndex]*xPower*polyFactor2/polyFactor1); 
        }  
        if(dim2-j-2>=0) {
          primeFactor[dim2-j-2]=1.0;
        }
        if(j==dim2-1) {
          d2=basis[k*dim+dimTwo-2];
          d1=basis[k*dim+dimTwo-1];
          polyTemp1=x[dim1+1]/x[dim1];
          polyTemp2=x[dim1];
          polyTemp3=x[dim1+1]; 
          polyFactor1=1.0;
          polyFactor3=1.0;
          polyFactor5=1.0;
          polyFactor2=1.0;
          polyFactor4=1.0;

          for(i=0; i<xPowerPrev; i++) {
            polyFactor1=polyFactor1*polyTemp1;
            polyFactor3=polyFactor3*polyTemp3;
          }
          for(i=0; i<xPower; i++) {
            polyFactor2=polyFactor2*polyTemp2;
          }
          for(i=1; i<xPower; i++) {
            polyFactor4=polyFactor4*polyTemp2;
          }
          for(i=0; i<xPowerPrev+1; i++) {
            polyFactor5=polyFactor5*polyTemp2;
          }

          v1= polyFactor1*jacobi[d2+jacobiIndexPrev]*(dxjacobi[d1+jacobiIndex]*(-2)*polyFactor2 \
          +jacobi[d1+jacobiIndex]*(xPower)*polyFactor4);
          v2= jacobi[d1+jacobiIndex]*polyFactor2*(polyFactor1*dxjacobi[d2+jacobiIndexPrev]*2.0*x[dim1+1]/(x[dim1]*x[dim1]) 
          +jacobi[d2+jacobiIndexPrev]*polyFactor3*(-xPowerPrev)/polyFactor5) ;
          primeFactor[0]= (v1+v2);   
        }
  
        if(j==dim2-2) {
          primeFactor[0]=1.0;
        }

        if(dim2>3&&j!=1&&j!=dim2-1&&j!=0) {
          d2=basis[k*dim+dim1+j-1];
          d1=basis[k*dim+dim1+j];
          prev=dimTwo-j-2;
          cur=dimTwo-j-1;
          next=dimTwo-j;

          polyTemp1=x[cur]/x[prev];
          polyTemp2=x[next]/x[cur];
          polyTemp3=x[next];
          polyTemp4=x[cur];
          polyFactor1=1.0;
          polyFactor2=1.0;
          polyFactor3=1.0;
          polyFactor4=1.0;
          polyFactor5=1.0;
          for(i=0; i<xPower; i++) {
            polyFactor1=polyFactor1*polyTemp1;
          }
          for(i=0; i<xPower-1; i++) {
            polyFactor5=polyFactor5*polyTemp1;
          }
          for(i=0; i<xPowerPrev+1; i++) {
            polyFactor4=polyFactor4*polyTemp4;
          }
          for(i=0; i<xPowerPrev; i++) {
            polyFactor2=polyFactor2*polyTemp2;
            polyFactor3=polyFactor3*polyTemp3;
          }
    
          v2= jacobi[d1+jacobiIndex]*polyFactor1*(polyFactor2*dxjacobi[d2+jacobiIndexPrev]*2.0*x[next]/(x[cur]*x[cur])+ \
          jacobi[d2+jacobiIndexPrev]*polyFactor3*(-xPowerPrev)/polyFactor4   ) ;
          v1= polyFactor2*jacobi[d2+jacobiIndexPrev]*(dxjacobi[d1+jacobiIndex]*(-2)/x[prev]*polyFactor1+ \
          jacobi[d1+jacobiIndex]*(xPower)*polyFactor5*1.0/x[prev]); 
          primeFactor[dim2-1-j]=v1+v2;
        }//endif  
     
        indexPrev=j;
        xFactorPrev=xFactor; 
        xPowerPrev=xPower;
        jacobiIndexPrev=order*indexPrev*order+order*xPowerPrev;
        
        for(d=0; d<dim2; d++) {
          phiPrime[(d+dim1)*m+k]=phiPrime[(d+dim1)*m+k]*primeFactor[d];
        }  
      }//end j loop 
    }   
 
    for(d=0; d<dim1; d++) {
      for(k=0;k<m;k++) {  
        phiPrime[k+d*m]=phiPrime[k+d*m]*phi[k];  
      }
    } 
    if(dim>dimTwo){
      for(d=dimTwo; d<dim; d++) {
        for(k=0;k<m;k++) {  
          phiPrime[k+d*m]=phiPrime[k+d*m]*phi[k];  
        }
      }
    }   
    free(jacobi);
    free(dxjacobi);
  }//end if dim>=3
  free(primeFactor);
  free(polyX);
  free(jacobiX);
  free(dxjacobiPrev);
  free(jacobiPrev);
  free(phi);
} /* end PhiPrimeSimplexPolyhedralTwo */
 

#include "../headers/JacobiPoly.h"

/*
 * Evaluates the Jacobi polynomials and their first derivatives
 * by the recursion formula
 * Parameters
 *    order    computes polynomials of order 0..order-1
 *    x        point where the polynomials are evaluated
 *    p        function values of the Legendre polynomials
 *    dp       derivatives of Legendre polynomials
 */
void JacobiPoly(int order, double x, double alpha, double beta, double *p){
  int k,n;
  double fac1, fac2, fac3;

  p[0] = 1.0;
  if(order>1) {
    p[1] = (x-1.0)/2.0*(alpha+beta+2.0)+alpha+1.0;
  }
  for ( k=1; k<order-1; k++ ) {
    n=k+1;
    fac1 =  (2.0*n+alpha+beta-1.0)*((2.0*n+alpha+beta)*(2.0*n+alpha+beta-2.0)*x+alpha*alpha-beta*beta);
    fac2 = -2.0*(n+alpha-1.0)*(n+beta-1.0)*(2.0*n+alpha+beta);
    fac3=1.0/(2.0*n*(n+alpha+beta)*(2.0*n+alpha+beta-2.0));
    p[k+1] = fac3*(fac1*p[k] + fac2*p[k-1]);
  }    
} /* end JacobiPoly */

void JacobiPolyPrime(int order, double x, double alpha, double beta,  double *dp) {
  int k;
  double p[order];
  JacobiPoly(order, x, alpha+1, beta+1, p);
  dp[0]= 0.0;
  dp[1]=(alpha+beta+2.0)/2.0;
  for(k=2; k<order; k++) {
    dp[k]=0.5*(1.0+alpha+beta+k)*p[k-1];   
  }
} /*end JacobiPolyPrime */

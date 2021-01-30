#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "../headers/LegendrePoly.h"


/*
 * Evaluate the Legendre polynomials and their first derivatives
 * by the recursion formula
 * Parameters
 *    order    computes polynomials of order 0..order-1
 *    x        point where the polynomials are evaluated
 *    p        function values of the Legendre polynomials
 *    dp       derivatives of Legendre polynomials
 */
void legendrePoly(int order, double x,  double *p,double *dp){
  int k;
  double fac1, fac2;

  p[0] = 1;
  dp[0] = 0;
  p[1] = x;
  dp[1] = 1;

  for ( k=1; k<order-1; k++ ) {
    fac1 = (2.0*k+1.0)/(k+1.0);
    fac2 = k/(k+1.0);
    p[k+1] = fac1*x*p[k] - fac2*p[k-1];
    dp[k+1] = fac1*(p[k] + x*dp[k]) - fac2*dp[k-1];
  }
} /* legendrePoly */


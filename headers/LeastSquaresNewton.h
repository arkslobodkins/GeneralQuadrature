/* Arkadijs Slobodkins
 * SMU Mathematics
 * February 2021
 */

#ifndef LEAST_SQUARES_NEWTON_H
#define LEAST_SQUARES_NEWTON_H
struct functionHandles;
struct quadratureParameters;
void  dgels_( char *, int * ,int *, int *, double *, int *, 
double *, int *, double *, int *, int *); 
void LeastSquaresNewton(double *, double *,int, int *, int *, struct functionHandles *, struct quadratureParameters *);
double ScaledTwoNorm(int , double *);
double TwoNorm(int, double *);
double InfNorm(int, double *);
#endif



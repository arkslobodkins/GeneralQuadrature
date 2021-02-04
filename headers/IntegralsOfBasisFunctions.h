/* Arkadijs Slobodkins
 * SMU Mathematics
 * February 2021
 */

#ifndef INTEGRALS_OF_BASIS_FUNCTIONS_H
#define INTEGRALS_OF_BASIS_FUNCTIONS_H
struct quadratureParameters;
void IntegralsCube(double *, struct quadratureParameters *params);
void IntegralsSimplex(double *, struct quadratureParameters *params);
void IntegralsCubeSimplex(double *, struct quadratureParameters *params);
void IntegralsSimplexSimplex(double *, struct quadratureParameters *params);
void IntegralsCubeSimplexSimplex(double *, struct quadratureParameters *params);
#endif


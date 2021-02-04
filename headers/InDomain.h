/* Arkadijs Slobodkins
 * SMU Mathematics
 * February 2021
 */

#ifndef IN_DOMAIN_H
#define IN_DOMAIN_H
typedef int boolean;
boolean InCube(const int, const double  *, const int *);
boolean InSimplex(const int, const double  *, const int *);
boolean InCubeSimplex(const int, const double  *, const int *);
boolean InSimplexSimplex(const int, const double  *, const int *);
boolean InCubeSimplexSimplex(const int, const double  *, const int *);
#endif

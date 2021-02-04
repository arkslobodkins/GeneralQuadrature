/* Arkadijs Slobodkins
 * SMU Mathematics
 * February 2021
 */

#include <stdlib.h>
#include "../headers/ImplementDomain.h"
#include "../headers/SetDomain.h"
#include "../headers/SetParams.h"
#include "../headers/ComputeDomain.h"
#include "../headers/Structs.h"


/* ImplementCube
 * Initializes structure parameters
 * and passes them to ComputeCube where 
 * quadrature rules are computed and printed to a file. 
*/
void ImplementCube(int degree,  int dim) {

  struct functionHandles *domainFunctions=(struct functionHandles *)malloc(sizeof(struct functionHandles));
  struct quadratureParameters *params =(struct quadratureParameters *)malloc(sizeof(struct quadratureParameters));
  params->dimension=malloc(3*sizeof(int));
  SetCube(domainFunctions);
  SetParams(dim, 0, 0, degree, params);
  ComputeCube(domainFunctions, params);
  free(params->dimension);
  free(params);
  free(domainFunctions); 
} /* end ImplementCube */

/* ImplementSimplex
 * Initializes structure parameters
 * and passes them to ComputeSimplex where 
 * quadrature rules are computed and printed to a file. 
*/
void ImplementSimplex(int degree, int dim) {

  struct functionHandles *domainFunctions=(struct functionHandles *)malloc(sizeof(struct functionHandles));
  struct quadratureParameters *params =(struct quadratureParameters *)malloc(sizeof(struct quadratureParameters));
  params->dimension=malloc(3*sizeof(int));
  SetSimplex(domainFunctions);
  SetParams(dim, 0, 0,  degree, params);
  ComputeSimplex(domainFunctions, params);
  free(params->dimension);
  free(params);
  free(domainFunctions); 
} /* end ImplementSimplex */

/* ImplementCubeSimplex
 * Initializes structure parameters
 * and passes them to ComputeCubeSimplex where 
 * quadrature rules are computed and printed to a file. 
*/
void ImplementCubeSimplex(int degree, int dim1, int dim2) {

  struct functionHandles *domainFunctions=(struct functionHandles *)malloc(sizeof(struct functionHandles));
  struct quadratureParameters *params =(struct quadratureParameters *)malloc(sizeof(struct quadratureParameters));
  params->dimension=malloc(3*sizeof(int));
  SetCubeSimplex(domainFunctions);
  SetParams(dim1, dim2, 0, degree, params);
  ComputeCubeSimplex(domainFunctions, params);
  free(params->dimension);
  free(params);
  free(domainFunctions); 
} /* end ImplementCubeSimplex */

/* ImplementSimplexSimplex
 * Initializes structure parameters
 * and passes them to ComputeSimplexSimplex where 
 * quadrature rules are computed and printed to a file. 
*/
void ImplementSimplexSimplex(int degree, int dim1, int dim2) {

  struct functionHandles *domainFunctions=(struct functionHandles *)malloc(sizeof(struct functionHandles));
  struct quadratureParameters *params =(struct quadratureParameters *)malloc(sizeof(struct quadratureParameters));
  params->dimension=malloc(3*sizeof(int));
  SetSimplexSimplex(domainFunctions);
  SetParams(dim1, dim2, 0, degree, params);
  ComputeSimplexSimplex(domainFunctions, params);
  free(params->dimension);
  free(params);
  free(domainFunctions); 
} /* end ImplementSimplexSimplex */


/* ImplementCubeSimplexSimplex
 * Initializes structure parameters
 * and passes them to ComputeCubeSimplexSimplex where 
 * quadrature rules are computed and printed to a file. 
*/
void ImplementCubeSimplexSimplex(int degree, int dim1, int dim2, int dim3) {

  struct functionHandles *domainFunctions=(struct functionHandles *)malloc(sizeof(struct functionHandles));
  struct quadratureParameters *params =(struct quadratureParameters *)malloc(sizeof(struct quadratureParameters));
  params->dimension=malloc(3*sizeof(int));
  SetCubeSimplexSimplex(domainFunctions);
  SetParams(dim1, dim2, dim3, degree, params);
  ComputeCubeSimplexSimplex(domainFunctions, params);
  free(params->dimension);
  free(params);
  free(domainFunctions); 
} /* end ImplementCubeSimplexSimplex */



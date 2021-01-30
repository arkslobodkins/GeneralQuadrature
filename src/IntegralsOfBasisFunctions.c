#include <string.h>
#include <math.h>
#include "../headers/IntegralsOfBasisFunctions.h"
#include "../headers/Structs.h"

/* Routines that compute analytical integrals of all orthogonal basis functions 
 * for respective domain and parameters. 
*/

void IntegralsCube(double *integrals, struct quadratureParameters *params){
  int i;
  int m=params->numberOfBasisFunctions;
  memset(integrals, 0, m*sizeof(integrals[0]));
  integrals[0]=1;
} /*end IntegralsCube */

void IntegralsSimplex(double *integrals,struct quadratureParameters *params){
  int i;
  int m=params->numberOfBasisFunctions;
  int dim1=params->dimension[0];

  memset(integrals, 0, m*sizeof(integrals[0]));
  integrals[0]=1;
  for(i=1; i<=dim1;i++) {
    integrals[0]=integrals[0]/i;  
  }
} /* end IntergalsSimplex */

void IntegralsCubeSimplex(double *integrals,struct quadratureParameters *params){
  int i;
  int m=params->numberOfBasisFunctions;
  int dim2=params->dimension[1];
  
  memset(integrals, 0, m*sizeof(integrals[0]));
  integrals[0]=1;
  for(i=1; i<=dim2;i++) {
    integrals[0]=integrals[0]/i;  
  }
} /* end IntegralsCubeSimplex */

void IntegralsSimplexSimplex(double *integrals,struct quadratureParameters *params){
  int i;
  int m=params->numberOfBasisFunctions;
  int dim1=params->dimension[0];
  int dim2=params->dimension[1];

  memset(integrals, 0, m*sizeof(integrals[0]));
  integrals[0]=1;
  for(i=1; i<=dim1;i++) {
    integrals[0]=integrals[0]/i;  
  }

  for(i=1; i<=dim2;i++) {
    integrals[0]=integrals[0]/i;  
  }
} /* end IntegralsSimplexSimplex */


void IntegralsCubeSimplexSimplex(double *integrals,struct quadratureParameters *params){
  int i;
  int m=params->numberOfBasisFunctions;
  int dim2=params->dimension[0];
  int dim3=params->dimension[1];

  memset(integrals, 0, m*sizeof(integrals[0]));
  integrals[0]=1;
  for(i=1; i<=dim2;i++) {
    integrals[0]=integrals[0]/i;  
  }

  for(i=1; i<=dim3;i++) {
    integrals[0]=integrals[0]/i;  
  }
} /* end IntegralsCubeSimplexSimplex */

#include <stdlib.h>
#include "../headers/InDomain.h"

/* InCube
 * Tests whether nodes are inside unit cube. 
*/
boolean InCube(const int numberOfNodes, const double *nodes, const int *dimension) {
  int i,j;
  int dim=dimension[0];
  boolean domain=1;

  for(i=0; i<numberOfNodes;i++) { 
    for(j=0; j<dim;j++) {

      if(nodes[dim*i+j]<0 || nodes[dim*i+j]>1) {
        domain=0;
        goto returnValue;
      }
    } 
  }

  returnValue: return domain;
} /* end InCube */

/* InSimplex
 * Tests whether nodes are inside unit simplex. 
*/
boolean InSimplex(const int numberOfNodes, const double *nodes, const int *dimension) {
  int i,j;
  int dim=dimension[0];
  boolean domain=1;

  for(i=0; i<numberOfNodes;i++) { 
    for(j=0; j<dim;j++) {
      if(nodes[dim*i+j]<0||nodes[dim*i+j]>1) {
        domain=0;
        goto returnValue;
      }
    } 

    for(j=1; j<dim;j++) { 
      if(nodes[dim*i+j]>nodes[dim*i+j-1]) {
        domain=0;
        goto returnValue;
      }
    }
  }

  returnValue: return domain;
} /* end InSimplex */

/* InCubeSimplex
 * Tests whether nodes are inside (unit cube of dimension[0]) X (unit simplex of dimension[1]). 
*/
boolean InCubeSimplex(const int numberOfNodes, const double *nodes, const int *dimension) {
  int i,j;
  int dimCube=dimension[0];
  int dimSimplex=dimension[1];
  int dimTotal=dimCube+dimSimplex;
  double nodesCube[numberOfNodes*dimCube];
  double nodesSimplex[numberOfNodes*dimSimplex];
  boolean domain;
  for(i=0; i<numberOfNodes;i++) { 
    for(j=0; j<dimCube;j++) {
      nodesCube[i*dimCube+j]=nodes[i*dimTotal+j];
    } 
    for(j=0; j<dimSimplex;j++) {  
      nodesSimplex[i*dimSimplex+j]=nodes[i*dimTotal+dimCube+j];
    }
  }
  domain=InCube(numberOfNodes, nodesCube, &dimCube)&InSimplex(numberOfNodes, nodesSimplex, &dimSimplex);

  return domain;
} /* end InCubeSimplex */

/* InSimplexSimplex
 * Tests whether nodes are inside (unit simplex of dimension[0]) X (unit simplex of dimension[1]). 
*/
boolean InSimplexSimplex(const int numberOfNodes, const double *nodes, const int *dimension) {
  int i,j;
  int dimSimplex1=dimension[0];
  int dimSimplex2=dimension[1];
  int dimTotal=dimSimplex1+dimSimplex2;
  double *nodesSimplex1=malloc(numberOfNodes*dimSimplex1*sizeof(double));
  double *nodesSimplex2=malloc(numberOfNodes*dimSimplex2*sizeof(double));
  boolean domain;
  for(i=0; i<numberOfNodes;i++) { 
    for(j=0; j<dimSimplex1;j++) {
      nodesSimplex1[i*dimSimplex1+j]=nodes[i*dimTotal+j];
    } 
    for(j=0; j<dimSimplex2;j++) { 
      nodesSimplex2[i*dimSimplex2+j]=nodes[i*dimTotal+dimSimplex1+j];
    }
  }

  domain=InSimplex(numberOfNodes, nodesSimplex1, &dimSimplex1) & InSimplex(numberOfNodes, nodesSimplex2, &dimSimplex2);
  free(nodesSimplex1);
  free(nodesSimplex2);
  return domain;
} /* end InSimplexSimplex */


/* InCubeSimplexSimplex
 * Tests whether nodes are inside (unit simplex of dimension[0]) X (unit simplex of dimension[1]) X (unit cube of dimension[2]). 
*/
boolean InCubeSimplexSimplex(const int numberOfNodes, const double *nodes, const int *dimension) {
  int i,j;
  int dimCube=dimension[2];
  int dimSimplex1=dimension[0];
  int dimSimplex2=dimension[1];
  int dimTotal=dimCube+dimSimplex1+dimSimplex2;
  
  double *nodesCube=malloc(numberOfNodes*dimCube*sizeof(double));
  double *nodesSimplex1=malloc(numberOfNodes*dimSimplex1*sizeof(double));
  double *nodesSimplex2=malloc(numberOfNodes*dimSimplex2*sizeof(double));
  boolean domain;

  for(i=0; i<numberOfNodes;i++) { 
    for(j=0; j<dimSimplex1;j++) {
      nodesSimplex1[i*dimSimplex1+j]=nodes[i*dimTotal+j];
    } 
    for(j=0; j<dimSimplex2;j++) { 
      nodesSimplex2[i*dimSimplex2+j]=nodes[i*dimTotal+dimSimplex1+j];
    }
    for(j=0; j<dimCube;j++) {
      nodesCube[i*dimCube+j]=nodes[i*dimTotal+dimSimplex1+dimSimplex2+j];
    } 
  }

  domain=InCube(numberOfNodes, nodesCube, &dimCube) & InSimplex(numberOfNodes, nodesSimplex1, &dimSimplex1) \
  & InSimplex(numberOfNodes, nodesSimplex2, &dimSimplex2);
  free(nodesCube);
  free(nodesSimplex1);
  free(nodesSimplex2);

  return domain;
} /* InCubeSimplexSimplex */

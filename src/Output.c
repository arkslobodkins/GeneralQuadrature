/* Arkadijs Slobodkins
 * SMU Mathematics
 * February 2021
 */

#include <stdio.h>
#include "../headers/Output.h"
#include "../headers/Structs.h"
#include "../headers/TestIntegral.h"

/* Output 
 * xInititial: initial nodes
 * wInitial: initial weights
 * nInitial: initial number of nodes and weights
 * xFinal: final nodes
 * wFinal: final weights
 * nFinal: final number of nodes and weights
 * functions: structure that contains function handles for the respective domain
 * params: structure that contains quadrature parameters
 * eliminationHistory: structure that contains information about every successful iteration of the Node Elimination
 * 
 * Prints history and quadrature parameters to a file. 
*/
void Output(double* xInitial, double* wInitial, int nInitial, double *xFinal, double *wFinal, int nFinal, struct functionHandles *functions, struct quadratureParameters *params, struct eliminationHistory *history, char *shape) {
  char str[50];
  int i,d;
  int eliminations=history->totalEliminations;
  int dim=params->totalDimension;
  int deg=params->degreeOfPrecision;
  int m=params->numberOfBasisFunctions;
  double res;
  FILE* FID;
  sprintf(str, "../results/history_%s_dim%i_deg%i.txt", shape, dim, deg);
  FID = fopen(str,"w");
    
  fprintf(FID, "degree of precision=%i \n", deg);
  fprintf(FID, "number of basis functions=%i \n", m);
  fprintf(FID, "dimension=%i \n", dim);
  fprintf(FID, "initial number of nodes=%i \n", nInitial);
  fprintf(FID, "final number of nodes=%i \n \n", nFinal);
   
  for(i=0; i<eliminations; i++) {
    fprintf(FID,"total number of nodes[%i]=%i \n", i, history->nodesTotal[i]);
    fprintf(FID,"successNodexIndex[%i]=%i \n", i, history->successNodeIndex[i]);
    fprintf(FID,"Converged in %i iterations \n",history->successNewtonIterations[i]);
    fprintf(FID," \n");
  }
  res=TestIntegral(nFinal, xFinal, wFinal, functions, params);   
  fprintf(FID,"final residual=%0.16f \n", res);

  fclose(FID);
} /* end output */

/* DumpCubatureRule
 * xFinal: final nodes
 * wFinal: final weights
 * nFinal: final number of nodes and weights
 * dim: dimension of the quadrature rule
 *  
 * prints final quadrature to a file 
*/  
void DumpCubatureRule(double *xFinal, double *wFinal, int nFinal,  struct quadratureParameters *params, char  *shape) {
  int i, j;
  FILE* FID;
  char str[50];
  int dim=params->totalDimension;
  int deg=params->degreeOfPrecision;
  sprintf(str, "../results/quadrature_%s_dim%i_deg%i.dat", shape, dim, deg);
  FID = fopen(str,"w");
  fprintf(FID, "%i %i \n", dim, nFinal);
  //print new quadrature nodes and weights    
  for(i=0;i<nFinal;i++) { 
    for(j=0;j<dim;j++) {
      fprintf(FID, "%.16le ", xFinal[(dim)*i+j]);         
    }
    fprintf(FID, "%.16le\n",wFinal[i]);
  }
  fclose(FID);
 
} /* end dumpCubatureRule */

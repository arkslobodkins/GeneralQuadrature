/* Arkadijs Slobodkins
 * SMU Mathematics
 * February 2021
 */

#ifndef STRUCTS_H
#define STRUCTS_H
#include "InDomain.h"
struct quadratureParameters{ 
  int *dimension;
	int  totalDimension;
  int degreeOfPrecision;
	int numberOfBasisFunctions;
}P;

struct functionHandles {
	void(*EvalBasis)(const double *const, double *, const struct quadratureParameters *);
	void(*EvalBasisDerivatives)(const double *const, double *, const struct quadratureParameters *);
	void(*IntegralsOfBasisFunctions)(double *, struct quadratureParameters *);
	boolean(*InDomain)(const int, const double *, const int *);
} F;

struct eliminationHistory {
	int *nodesTotal;
	int *successNodeIndex;
	int *successNewtonIterations;
	int totalEliminations;
} E;
#endif

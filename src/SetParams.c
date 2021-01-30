#include <stdlib.h>
#include "../headers/SetParams.h"
#include "../headers/BasisIndices.h"
#include "../headers/Structs.h"

void SetParams(int dim1, int dim2 , int dim3,  int degree,  struct quadratureParameters *params) {	
  params->dimension[0]=dim1;
  params->dimension[1]=dim2;
  params->dimension[2]=dim3;
  params->totalDimension=dim1+dim2+dim3;
  params->degreeOfPrecision=degree;
  params->numberOfBasisFunctions= BasisSize((params)->totalDimension, degree);
} /* end SetParams */

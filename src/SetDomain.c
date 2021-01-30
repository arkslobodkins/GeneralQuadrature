#include <stdlib.h>
#include "../headers/SetDomain.h"
#include "../headers/SetParams.h"
#include "../headers/Phi.h"
#include "../headers/InDomain.h"
#include "../headers/IntegralsOfBasisFunctions.h"
#include "../headers/Structs.h"

void SetCube(struct functionHandles *domainFunctions) {
  domainFunctions->EvalBasis=&PhiCube;
  domainFunctions->EvalBasisDerivatives=&PhiPrimeCube;
  domainFunctions->InDomain=&InCube;
  domainFunctions->IntegralsOfBasisFunctions=&IntegralsCube;
} /* end SetCube */

void SetSimplex(struct functionHandles *domainFunctions) {
  domainFunctions->EvalBasis=&PhiSimplex;
  domainFunctions->EvalBasisDerivatives=&PhiPrimeSimplex;
  domainFunctions->InDomain=&InSimplex;
  domainFunctions->IntegralsOfBasisFunctions=&IntegralsSimplex;
} /* end SetSimplex */

void SetCubeSimplex(struct functionHandles *domainFunctions) {
  domainFunctions->EvalBasis=&PhiCubeSimplex;
  domainFunctions->EvalBasisDerivatives=&PhiPrimeCubeSimplex;
  domainFunctions->InDomain=&InCubeSimplex;
  domainFunctions->IntegralsOfBasisFunctions=&IntegralsCubeSimplex;
} /* end SetCubeSimplex */

void SetSimplexSimplex(struct functionHandles *domainFunctions) {
  domainFunctions->EvalBasis=&PhiSimplexSimplex;
  domainFunctions->EvalBasisDerivatives=&PhiPrimeSimplexSimplex;
  domainFunctions->InDomain=&InSimplexSimplex;
  domainFunctions->IntegralsOfBasisFunctions=&IntegralsSimplexSimplex;
} /* end SetSimplexSimplex */

void SetCubeSimplexSimplex(struct functionHandles *domainFunctions) {
  domainFunctions->EvalBasis=&PhiCubeSimplexSimplex;
  domainFunctions->EvalBasisDerivatives=&PhiPrimeCubeSimplexSimplex;
  domainFunctions->InDomain=&InCubeSimplexSimplex;
  domainFunctions->IntegralsOfBasisFunctions=&IntegralsCubeSimplexSimplex;
} /* end SetCubeSimplexSimplex */


/* Arkadijs Slobodkins
 * SMU Mathematics
 * February 2021
 */

#ifndef COMPUTE_DOMAIN_H
#define COMPUTE_DOMAIN_H
struct functionHandles;
struct quadratureParameters;        
void ComputeCube(struct functionHandles *, struct quadratureParameters *);
void ComputeSimplex(struct functionHandles *, struct quadratureParameters *);
void ComputeCubeSimplex(struct functionHandles *, struct quadratureParameters *);
void ComputeSimplexSimplex(struct functionHandles *, struct quadratureParameters *);
void ComputeCubeSimplexSimplex(struct functionHandles *, struct quadratureParameters *);
#endif

#ifndef GET_JACOBIAN_H
#define GET_JACOBIAN_H
struct functionHandles;
struct quadratureParameters;
void GetJacobian(int ,double *, double *, double *,  struct functionHandles *, struct quadratureParameters *);
#endif

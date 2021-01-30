#ifndef PHI_H
#define PHI_H
struct quadratureParameters;
void PhiCube(const double *const, double *, const struct quadratureParameters *);
void PhiPrimeCube(const double *const, double *, const struct quadratureParameters *);
void PhiSimplex(const double *const, double *, const struct quadratureParameters *);
void PhiPrimeSimplex(const double *const, double *, const struct quadratureParameters *);
void PhiCubeSimplex(const double *const, double *, const struct quadratureParameters *);
void PhiPrimeCubeSimplex(const double *const, double *, const struct quadratureParameters *);
void PhiSimplexSimplex(const double *const, double *, const struct quadratureParameters *);
void PhiPrimeSimplexSimplex(const double *const, double *, const struct quadratureParameters *);
void PhiCubeSimplexSimplex(const double *const, double *, const struct quadratureParameters *);
void PhiPrimeCubeSimplexSimplex(const double *const, double *, const struct quadratureParameters *);
void PhiSimplexPolyhedralOne(const double *const, double *, const struct quadratureParameters *);
void PhiSimplexPolyhedralTwo(const double *const, double *, const struct quadratureParameters *);
void PhiPrimeSimplexPolyhedralOne(const double *const, double *, const struct quadratureParameters *);
void PhiPrimeSimplexPolyhedralTwo(const double *const, double *, const struct quadratureParameters *);
#endif


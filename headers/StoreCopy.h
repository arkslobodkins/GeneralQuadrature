/* Arkadijs Slobodkins
 * SMU Mathematics
 * February 2021
 */

#ifndef STORE_COPY_H
#define STORE_COPY_H
void CopyNodes( int, int, double *, double *);
void CopyWeights( int, double *, double *);
void CopyNodesAndWeights(int, int, double *, double *);
void UnrollNodesAndWeights(int, int, double *, double *, double *);
void RollNodesAndWeights(int, int, double *, double *, double *);
#endif

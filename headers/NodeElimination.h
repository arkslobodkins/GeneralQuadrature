/* Arkadijs Slobodkins
 * SMU Mathematics
 * February 2021
 */

#ifndef NODE_ELIMINATION_H
#define NODE_ELIMINATION_H
struct functionHandles;
struct quadratureParameters;
struct eliminationHistory; 
void dgeqr2_(int *,int *, double *, int *, double *, double *, int *); 
void dormqr_(char *, char *, int *, int *, int *, double *, int *, double *, double *, int *, double *, int  *, int *); 	
void NodeElimination(int, double *,double *,   int *, double *, double *, struct functionHandles *, struct quadratureParameters *,struct eliminationHistory *);
#endif

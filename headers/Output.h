/* Arkadijs Slobodkins
 * SMU Mathematics
 * February 2021
 */

#ifndef OUTPUT_H
#define OUTPUT_H
struct functionHandles;
struct quadratureParameters;
struct eliminationHistory history; 
void Output(double*, double* , int , double *, double *, int , struct functionHandles *, struct quadratureParameters *, struct eliminationHistory *);
void DumpCubatureRule(double *, double *, int , int);
#endif

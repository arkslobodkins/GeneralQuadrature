#ifndef ADD_DIMENSION_H
#define ADD_DIMENSION_H
void AddLine(int d,int n, int Nold, double *x,double *Xnew, double *w, double *Wnew, double *Xinitial, double *Winitial);
void AddLineLast(int d,int n, int Nold, double *x,double *Xnew, double *w, double *Wnew, double *Xinitial, double *Winitial);
void AddLineSimplex(int d,int n, int Nold, double *x,double *Xnew, double *w, double *Wnew, double *Xinitial, double *Winitial);
#endif 

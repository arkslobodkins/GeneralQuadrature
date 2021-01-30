#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "../headers/ImplementDomain.h"

/* main
 * Receives domain type, dimension sizes, and degree
 * of precision as command line arguments. Performs tests
 * to test appropriateness of the input arguments. If inputs are 
 * valid, it proceeds to recursive initial guess procedure and Node Elimination 
 * algorithm.
*/

int main(int argc, char *argv[]) {
 
  char *SHAPE, *deg, *d1, *d2, *d3;
  char flag[0];
  int i, p;
  char *d[3];
  int dim[3];

  if(argv[1]!=NULL) {
    SHAPE=argv[1];
  }
  i=2;
  int counter=0;
  if(argv[2]==NULL) {
    printf("Not enough input parameters for %s, exiting program. \n", SHAPE);
    exit(-1);
  }
  while(argv[i+1]!=NULL&&counter<3) {
    d[counter]=argv[i];
    dim[counter]=atoi(argv[i]);
    i++;
    counter++;
  }
  deg=argv[i];
  p=atoi(deg);
  if(p<1) {
    printf("Degree of accuracy should be >=1, exiting program. \n");
    exit(-1);
  }
  int domainType=-1;
  int test=1;
  if(strcmp("CUBE", SHAPE)==0){
    
          
  if(argv[3]==NULL) {
    printf("Not enough input parameters for %s, exiting program. \n", SHAPE);
    exit(-1);
  }

    
  if(dim[0]<1){
    printf("Dimension for %s should be >=1. \n ", SHAPE);
    test=0;
  }
  if (test==0) {
    printf("Exiting program. \n");
    exit(-1);
  }
                

  if(argv[4] !=NULL) {
    printf("Warning: Received redundant parameters for %s. Enter character 'Y' if you would like to continue. Exiting otherwise. \n", SHAPE);
    scanf("%c", flag);
    if(flag[0]!='Y') {
      printf("Exiting program. \n");
      exit(-1);
    }
  }
                
  domainType=0;

  } else if(strcmp("SIMPLEX",SHAPE)==0){

    if(argv[3]==NULL) {
      printf("Not enough input parameters for %s, exiting program \n", SHAPE);
      exit(-1);
    }
 
    if(dim[0]<2){
      printf("Dimension for %s should be >=2. ", SHAPE);
      test=0;
    }
    if (test==0) {
      printf("Exiting program. \n");
      exit(-1);
    }
                
    if(argv[4] !=NULL) {
      printf("Warning: Received redundant parameters for %s. Enter character 'Y' if you would like to continue. Exiting otherwise. \n", SHAPE);
      scanf("%c", flag);
      if(flag[0]!='Y') {
        printf("Exiting program. \n");
        exit(-1);
      }
    }
    domainType=1;
  } 
  else if(strcmp("CUBESIMPLEX", SHAPE)==0){

    if(argv[4]==NULL || argv[3]==NULL) {
      printf("Not enough input parameters for %s, exiting program. \n", SHAPE);
      exit(-1);
    }

    if(dim[0]<1){
      printf("First dimension for %s should be >=1.\n", SHAPE);
      test=0;
    }
    if(dim[1]<2) {
      printf("Second dimension for %s should be >=2.\n", SHAPE);
      test=0;
    }
    if (test==0) {
      printf("Exiting program. \n");
      exit(-1);
    }
                
    if(argv[5] !=NULL) {
      printf("Warning: Received redundant parameters for %s. Enter character 'Y' if you would like to continue. Exiting otherwise. \n", SHAPE);
      scanf("%c", flag);
      if(flag[0]!='Y') {
        printf("Exiting program. \n");
        exit(-1);
      }
    }

    domainType=2;
  }
  else if(strcmp("SIMPLEXSIMPLEX", SHAPE)==0){

    if(argv[4]==NULL || argv[3] ==NULL) {
      printf("Not enough input parameters for %s, exiting program. \n", SHAPE);
      exit(-1);
    }
    if(dim[0]<2||dim[1]<2){
      printf("Both dimensions for %s should be >=2, exiting program. \n", SHAPE);
      exit(-1);
    }
    domainType=3;
    if(argv[5] !=NULL) {
      printf("Warning: Received redundant parameters for %s. Enter character 'Y' if you would like to continue. Exiting otherwise. \n", SHAPE);
      scanf("%c", flag);
      if(flag[0]!='Y') {
        printf("Exiting program \n");
        exit(-1);
      }
    }
  }
 
  else if(strcmp("CUBESIMPLEXSIMPLEX", SHAPE)==0){

    if(argv[5]==NULL || argv[4]==NULL ||argv[3]==NULL) {
      printf("Not enough input parameters for %s, exiting program \n", SHAPE);
      exit(-1);
    }

    if(dim[0]<1){
      printf("First dimension for %s should be >=1. \n", SHAPE);
      test=0;
    }
    if(dim[1]<2) {
      printf("Second dimension for %s should be >=2. \n", SHAPE);
      test=0;
    }
    if(dim[2]<2) {
      printf("Third dimension for %s should be >=2. \n", SHAPE);
      test=0;
    }

    if (test==0) {
      printf("Exiting program. \n");
      exit(-1);
    }

    if(argv[6] !=NULL) {
      printf("Warning: Received redundant parameters for %s. Enter character 'Y' if you would like to continue. Exiting otherwise. \n", SHAPE);
      scanf("%c", flag);
      if(flag[0]!='Y') {
        printf("Exiting program \n");
        exit(-1);
      }
    }
    domainType=4;
  }

  switch(domainType) {
    case 0:
      ImplementCube(p, dim[0]);
      break;
    case 1:
      ImplementSimplex(p, dim[0]);
      break;
    case 2:
      ImplementCubeSimplex(p, dim[0], dim[1]);
      break;
    case 3:
      ImplementSimplexSimplex(p, dim[0], dim[1]);
      break;
    case 4:
      ImplementCubeSimplexSimplex(p, dim[0], dim[1], dim[2]);
      break;
    default: 
      printf("Domain is not specified properly. Acceptable parameters:  \n \
      'CUBE' 'SIMPLEX' 'CUBESIMPLEX 'SIMPLEXSIMPLEX', 'CUBESIMPLEXSIMPLEX ");   
  }

return 0;       
} /* end main */


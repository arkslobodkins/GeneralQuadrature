/* Arkadijs Slobodkins
 * SMU Mathematics
 * February 2021
 */

#include "../headers/InsertionSort.h"

/* InsertionSort 
 * numberOfEntries: number of entries in array
 * array: array to be rearranged
 * arrayIndex: indices to be rearranged
 * Sorts current entries of array in ascending order    
 */
void InsertionSort(const int numberOfEntries, double *array,int *arrayIndex) {
  int i,j;
  double tempValue;

//Sort s in ascending order
  for(i=0;i<numberOfEntries;i++) {
    tempValue=array[i];
    j=i;
    while((j>0)&&(array[j-1]>tempValue)) {      
      array[j]=array[j-1];
      arrayIndex[j]=arrayIndex[j-1];
      j--;
    } 
    array[j]=tempValue;
    arrayIndex[j]=i;  
  }
} /*end Insertion sort */  
 


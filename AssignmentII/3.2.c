#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/time.h>
#include <stddef.h>
#include <stdint.h>
#include "omp.h"


#define sizet 10000000
double start; 
double end; 
double inputt[sizet];


void generate_random(double *input, size_t size)
{
  for (size_t i = 0; i < size; i++) {
    input[i] = rand() / (double)(RAND_MAX);
  }
}

double serial_sum(double *x, size_t size)
{
  double sum_val = 0.0;

  #pragma omp parallel for
  for (size_t i = 0; i < size; i++) {
    sum_val += x[i];
  }

  return sum_val;
}

void main(){
    
    double sum;
    generate_random(inputt, sizet);

    start = omp_get_wtime(); 

    sum=serial_sum(inputt,sizet);

    end = omp_get_wtime(); 

    printf("SERIAL  The sum of array is %f\n", sum);
    printf("Work took %f seconds\n", end - start);

}

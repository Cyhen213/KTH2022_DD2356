#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/time.h>
#include <stddef.h>
#include <stdint.h>
#include "omp.h"


#define sizet 10000000
#define MAX_THREADS 32
double start; 
double end; 
double inputt[sizet];



void generate_random(double *input, size_t size)
{
  for (size_t i = 0; i < size; i++) {
    input[i] = rand() / (double)(RAND_MAX);
  }
}

double omp_local_sum(double *x, size_t size)
{
  double sum_val[MAX_THREADS], mval;
  #pragma omp parallel shared(mval)
  {
    int id = omp_get_thread_num();
    //printf("The number of threads is %d\n", id);
    sum_val[id] = 0.0;

    #pragma omp for
    for (size_t i = 0; i < size; i++) {
      sum_val[id] += x[i];
    }
  }
  mval=sum_val[0];
  for (int i=1; i<MAX_THREADS; i++){
      mval += sum_val[i];
  }

  return mval;
}

void main(){
    
    double sum;
    generate_random(inputt, sizet);

    start = omp_get_wtime(); 

    sum=omp_local_sum(inputt,sizet);

    end = omp_get_wtime(); 

    printf("USING MAX_THREADS The sum of array is %f\n", sum);
    //printf("The number of threads is %d\n", id);
    printf("Work took %f seconds\n", end - start);

}

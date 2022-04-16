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
typedef struct { double val; char pad[128];} tvals;
tvals suminfo[MAX_THREADS];


void generate_random(double *input, size_t size)
{
  for (size_t i = 0; i < size; i++) {
    input[i] = rand() / (double)(RAND_MAX);
  }
}

double opt_local_sum(double *x, size_t size)
{
  //tvals suminfo[MAX_THREADS];//called function cannot allocate the space
  double tval;
  #pragma omp parallel shared(suminfo)
  {
    int id = omp_get_thread_num();

    suminfo[id].val = 0.0;

    #pragma omp for
    for (size_t i = 0; i < size; i++) {
      suminfo[id].val += x[i];
    }
  }
  tval=suminfo[0].val;
  for (int i=1; i<MAX_THREADS; i++){
      tval += suminfo[i].val;
  }

  return tval;
}

void main(){

    double sum;
    generate_random(inputt, sizet);

    start = omp_get_wtime(); 

    sum=opt_local_sum(inputt,sizet);

    end = omp_get_wtime(); 

    printf("USING OPT MAX_THREADS The sum of array is %f\n", sum);
    //printf("The number of threads is %d\n", id);
    printf("Work took %f seconds\n", end - start);

}

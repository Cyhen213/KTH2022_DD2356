#include <stdio.h>
#include <stdlib.h>
#include "/Users/gaogao/opt/anaconda3/pkgs/llvm-openmp-12.0.0-h0dcd299_1/lib/clang/12.0.0/include/omp.h"
#define MAX_THREADS 20
void generate_random();
int main()
{
    double start; 
    double end; 
    size_t size=10000000;
    int a;
    double *x;
    x=(double*)malloc(size*sizeof(double));
    double local_sum[MAX_THREADS];
    double lsum;
    generate_random(x,size);
    omp_set_num_threads(MAX_THREADS);
    #pragma omp parallel shared(local_sum)
    {
      int id=omp_get_thread_num();
      local_sum[id]=0;
      start = omp_get_wtime(); 
      #pragma omp for
        for (size_t i = 0; i < size; i++) 
        { 
          local_sum [id]+= x[i];
        }
        lsum=local_sum[0];
        for(int i=0;i<MAX_THREADS;i++)
        {
          lsum+=local_sum[i];
        }
      end = omp_get_wtime(); 
    }
    printf("Work took %f seconds\n", end - start);
}
void generate_random(double *input, size_t size)
{

  for (size_t i = 0; i < size; i++) {
    input[i] = rand() / (double)(RAND_MAX);
  }
}

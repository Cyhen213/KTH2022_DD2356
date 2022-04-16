#include <stdio.h>
#include <stdlib.h>
#include "/Users/gaogao/opt/anaconda3/pkgs/llvm-openmp-12.0.0-h0dcd299_1/lib/clang/12.0.0/include/omp.h"
double serial_sum();
void generate_random();
int main()
{
    double start; 
    double end; 
    size_t size=10000000;
    int a;
    double *input;
    input=(double*)malloc(size*sizeof(double));
    generate_random(input,size);
    omp_set_num_threads(1);

    start = omp_get_wtime(); 
    a=serial_sum(input,size);
    end = omp_get_wtime(); 
    printf("Work took %f seconds\n", end - start);
}
void generate_random(double *input, size_t size)
{

  for (size_t i = 0; i < size; i++) {
    input[i] = rand() / (double)(RAND_MAX);
  }
}

double serial_sum(double *x, size_t size)
{
  double sum_val = 0.0;
  #pragma omp parallel 
  {
    #pragma omp for
    for (size_t i = 0; i < size; i++) 
    {
      #pragma omp critical
      {
        sum_val += x[i];
      }
    }
  }
  return sum_val;
} 
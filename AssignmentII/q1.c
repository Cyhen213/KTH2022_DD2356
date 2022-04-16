#include <stdio.h>
#include "/Users/gaogao/opt/anaconda3/pkgs/llvm-openmp-12.0.0-h0dcd299_1/lib/clang/12.0.0/include/omp.h"
int main(int argc, char* argv[]) {
    omp_set_num_threads(4);
    #pragma omp parallel
    {
        int id=omp_get_thread_num();
        printf("Hello world from %d\n",id);
    }
    return 0;
}

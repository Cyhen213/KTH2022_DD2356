
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <mpi.h>

#define SEED     921
#define NUM_ITER 1000000000

int main(int argc, char* argv[])
{
    int count = 0;
    double x, y, z, pi;
    double start,end;

    int rank, size, i, provided, NUM_perth;
    int Count_sum;

    MPI_Init_thread(&argc, &argv, MPI_THREAD_SINGLE, &provided);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    srand(SEED*rank); // Important: Multiply SEED by "rank" when you introduce MPI!
    NUM_perth=NUM_ITER/size;
    
    
    start = MPI_Wtime();



        //Calculate PI following a Monte Carlo method
    for (int iter = 0; iter < NUM_perth; iter++)
    {
        // Generate random (X,Y) points
        x = (double)random() / (double)RAND_MAX;
        y = (double)random() / (double)RAND_MAX;
        z = sqrt((x*x) + (y*y));
        
        // Check if point is in unit circle
        if (z <= 1.0)
        {
            count++;
        }
    }

    //double start = MPI_Wtime();
    
    MPI_Reduce(&count, &Count_sum, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);  //This time there is no difference in rank0



    if(rank==0){
            // Estimate Pi and display the result
        pi = ((double)count / (double)NUM_ITER) * 4.0;
    }

        end = MPI_Wtime();
        double time = end - start; 
        
    if(rank==0){
        printf("The result is Pi = %f\n", pi);
        printf("Execution time is %f\n", time);
    }



    MPI_Finalize();
    return 0;
}


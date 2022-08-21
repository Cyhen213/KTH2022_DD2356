#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

int main(int argc, char *argv[])
{
	/* -------------------------------------------------------------------------------------------
		MPI Initialization 
	--------------------------------------------------------------------------------------------*/
	int size,rank;
        MPI_Init(&argc,&argv);
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &size);
        MPI_Status stat;
     
	if(size!= 2){
		if(rank == 0){
			printf("This program requires exactly 2 MPI ranks, but you are attempting to use %d! Exiting...\n",size);
		}
		MPI_Finalize();
		exit(0);
	}

	/* -------------------------------------------------------------------------------------------
		Loop from 8 B to 1 GB
	--------------------------------------------------------------------------------------------*/

	for(int i=0; i<=27; i++){
		long int N = 1 << i;
                double *local_mem;
                double *shared_mem;
                MPI_Alloc_mem(N * sizeof(double), MPI_INFO_NULL, &local_mem);
	        MPI_Alloc_mem(N * sizeof(double), MPI_INFO_NULL, &shared_mem);
                for(int i=0;i<N;++i){
                  local_mem[i]=1;
                  shared_mem[i]=0;
               }
               MPI_Win win;
	       MPI_Win_create(local_mem, sizeof(double)*N,sizeof(double),MPI_INFO_NULL, MPI_COMM_WORLD, &win);
           
	       int loop_count = 50;
	       double start_time, stop_time, elapsed_time;
              
               start_time= MPI_Wtime();
               for(int i=1; i<=loop_count; i++){
                  MPI_Win_fence(0, win);
                  MPI_Get(shared_mem, N, MPI_DOUBLE, 0,(rank+1)%size, N, MPI_DOUBLE, win); 
                  MPI_Win_fence(0, win);		
	       }
	       stop_time = MPI_Wtime();
	       elapsed_time = stop_time - start_time;
	       long int num_B = 8*N;
	       long int B_in_GB = 1 << 30;
	       double num_GB = (double)num_B / (double)B_in_GB;
	       double avg_time_per_transfer = elapsed_time / (2.0*(double)loop_count);
	       if(rank == 0) printf("%10li\t%15.9f\n", num_B, avg_time_per_transfer);
	       free(local_mem);
       }
       MPI_Finalize();
       return 0;
}

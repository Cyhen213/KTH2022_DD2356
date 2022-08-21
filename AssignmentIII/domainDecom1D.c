
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(int argc, char *argv[]){

    int rank, size, i, provided;
    
    // number of cells (global)
    int nxc = 128; // make sure nxc is divisible by size
    double L = 2*3.1415; // Length of the domain
    double error = 0.01;
    

    MPI_Init_thread(&argc, &argv, MPI_THREAD_SINGLE, &provided);

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // number of nodes (local to the process): 0 and nxn_loc-1 are ghost cells /////////
    int nxn_loc = nxc/size + 3; // number of nodes is number cells + 1; we add also 2 ghost cells
    double L_loc = L/((double) size);
    double dx = L / ((double) nxc);
    
    // define out function
    double *f = calloc(nxn_loc, sizeof(double)); // allocate and fill with z
    double *dfdx = calloc(nxn_loc, sizeof(double)); // allocate and fill with z

    for (i=1; i<(nxn_loc-1); i++)
      f[i] = sin(L_loc*rank + (i-1) * dx);
    
    // need to communicate and fill ghost cells f[0] and f[nxn_loc-1]
    // communicate ghost cells
    // ...
    // ...  
    int dest_prev,dest_next;
    if (rank == 0) {
        dest_prev = size-1;
        dest_next = 1;    
    }
    else if (rank == size-1) {
        dest_prev = rank - 1;
        dest_next = 0;
    }
    else {
        dest_prev = rank - 1;
        dest_next = rank + 1;
    }

    MPI_Send(&f[2], 1, MPI_DOUBLE, dest_prev, 2, MPI_COMM_WORLD);  //
    MPI_Recv(&f[nxn_loc-1], 1, MPI_DOUBLE, dest_next, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Send(&f[nxn_loc-3], 1, MPI_DOUBLE, dest_next, 1, MPI_COMM_WORLD);
    MPI_Recv(&f[0], 1, MPI_DOUBLE, dest_prev, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE); 
    // here we finish the calculations

    // calculate first order derivative using central difference
    // here we need to correct value of the ghost cells!
    for (i=1; i<(nxn_loc-1); i++)
      dfdx[i] = (f[i+1] - f[i-1])/(2*dx);

    
    // Print f values
    if (rank==0){ // print only rank 0 for convenience
        printf("My rank %d of %d\n", rank, size );
        printf("Here are my values for f including ghost cells\n");
        for (i=0; i<nxn_loc; i++)
	    printf("%f\n", f[i]);
        printf("\n");
    }


        //Testing values of the ghost cells and domain edge deviative

    if (fabs(f[0]-sin(L_loc*rank-dx))>error || fabs(f[nxn_loc-1]-sin(L_loc*(rank+1)+dx))>error) {
        printf("Check NOT OK in rank %d\n",rank);
    }
    else printf("Check OK in rank %d \n", rank);

    if(rank==0){
      if(fabs(dfdx[1]-cos(0))>error){
         printf("Left derivative check NOT OK in rank \n");
      }
      else printf("Left derivative check OK in rank \n");
    }

    if(rank==size-1){
      if(fabs(dfdx[nxn_loc-2]-cos(L))>error){
         printf("Right derivative check NOT OK in rank \n");
      }
      else printf("Right derivative check OK in rank \n");
    }
      
      

    MPI_Finalize();
}






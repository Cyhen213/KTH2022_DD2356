#include <mpi.h>  
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <time.h>

#define SEED     921
#define NUM_ITER 1000000000

int main(int argc, char* argv[])
{
    int count = 0;
    int count_buf=0;
    int Count_sum=0;
    double x, y, z, pi;

    int rank, size,  provided, NUM_perth;
    int i,j,power,k;

    /*for(i=0;i<100;i++){
        if (size==(int)pow(2,i)){
            power=i;
            printf("power is %d\n", power);
            break;
        }
    }*/

    MPI_Init_thread(&argc, &argv, MPI_THREAD_SINGLE, &provided);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    srand(SEED*rank); // Important: Multiply SEED by "rank" when you introduce MPI!
    NUM_perth=NUM_ITER/size;
    
        for(j=0;j<100;j++){
        if (size==(int)pow(2,j)){
            power=j;
            //printf("power is %d\n", power);
            break;
        }
    }
    int once_flag=1;


    double start = MPI_Wtime();

        //Calculate PI following a Monte Carlo method
    for (int iter=0; iter<NUM_perth; iter++)
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


    /*
    for(i=1;i<power+1;i++){

        if((rank%(int)pow(2,i-1)==0)&&(rank%(int)pow(2,i)!=0)){                         
            MPI_Send(&count, 1, MPI_INT, rank-(int)pow(2,i-1), 0, MPI_COMM_WORLD);
        }
        else if (rank%(int)pow(2,i)==0) {
            //int count1;
            //for(j=1;j<size/pow(2,i);j++){
                MPI_Recv(&count_buf, 1, MPI_INT, rank+(int)pow(2,i-1), 0, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
                printf("saving rank is %d\n", rank);
                count+=count_buf;
            //}
        }
    }*/


     

      
        for(i=1;i<power+1;i++){   //2,4,8,...

        if (once_flag==1){
            int difference = (int)pow(2,i);
            if(rank%(int)pow(2,i)!=0){
                MPI_Send(&count, 1, MPI_INT, rank-(int)pow(2,i-1), 0, MPI_COMM_WORLD);
                once_flag=0;
                //printf("rank is %d\n", rank);
            }
        else {

            MPI_Recv(&count_buf, 1, MPI_INT, rank+(int)(pow(2,i-1)), 0, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            count+=count_buf;
            //if(i==power) {count=Count_sum;}
            //printf("saving rank is %d\n", rank);

        }
        }  
    }
   



    if(rank==0){
            // Estimate Pi and display the result
        pi = ((double)count / (double)NUM_ITER) * 4.0;
    }

        double end = MPI_Wtime();
        
    if(rank==0){
        printf("The result is Pi = %f\n", pi);
        double time = end - start; 
        printf("Execution time is %f\n", time);
    }

    
    MPI_Finalize();
    return 0;
}
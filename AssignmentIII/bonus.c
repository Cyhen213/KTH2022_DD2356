#include <stdio.h>
#include <stdlib.h>
//#include <mpi.h>
#include "/Users/gaogao/opt/anaconda3/pkgs/mpich-3.3.2-hc856adb_0/include/mpi.h"
#include <math.h>
#include <string.h>
#include <assert.h>
void print_mat(double *m,int n)
{
  int row,col;
  for(row=0;row<n;row++)
  {
    for(col=0;col<n;col++)
    {
      printf("%f",m[row*n+col]);
      printf("  ");
    }
    printf("\n");
  }
  printf("\n");
}

void generate_matrix(double *m,int n,int seed){
    int row,col;
  srand(seed);
    for(row=0;row<n;row++){
        for(col=0;col<n;col++){
            m[row*n+col]=((double)rand()*100/(double)(RAND_MAX));
        }
    }
}

void send_block_AB(double *MatrixA,double *MatrixB,double *blockA,double *blockB,double *tempA,double *tempB,int n,int block,int num_processor,int sqrt_num_processor,MPI_Comm processGrid)
{

    int row_min,row_max,col_min,col_max;
    for(int i=0;i<num_processor;i++)
    {
        row_min=(i/sqrt_num_processor)*block;
        row_max=row_min+block;
        col_min=(i%sqrt_num_processor)*block;
        col_max=col_min+block;
        for(int j=row_min;j<row_max;j++)
        {
            for(int k=col_min;k<col_max;k++)
            {
              int idx;
                idx=(j-row_min)*block+k-col_min;
                tempA[idx]=MatrixA[n*j+k];
                tempB[idx]=MatrixB[n*j+k];
            }
        }
        if (i==0)
        {
          memcpy(blockA,tempA,block*block*sizeof(double));
          memcpy(blockB,tempB,block*block*sizeof(double));
        }
        else
        {
          MPI_Send(tempA,block*block,MPI_DOUBLE,i,1,processGrid);
          MPI_Send(tempB,block*block,MPI_DOUBLE,i,2,processGrid);
        }
    }
}
void matmul(double *a,double *b,double *c,int n){
        for(int i=0;i<n;i++){
          for(int j=0;j<n;j++){
            for(int k=0;k<n;k++){
              c[n*i+j]+=a[n*i+k]*b[n*k+j];
            }
          }
      }
}

                 
int main(int argc, char *argv[])
{
            int num_processor;
            int parRank;
            double start, end;
            MPI_Init(&argc,&argv);
            MPI_Comm_size(MPI_COMM_WORLD,&num_processor);
            MPI_Comm_rank(MPI_COMM_WORLD, &parRank);
            MPI_Request request;
            MPI_Status stat;
            MPI_Comm processGrid,processRow,processCol;

            int seed=213;
            int n = atoi(argv[1]);
            int num_dim=2;
            int dimensions[2]={0,0};
            int periods[2]={1,1};
            int reorder=1;
            int coords[2];
            int comm_rank;
            MPI_Dims_create(num_processor,num_dim,dimensions);
            MPI_Cart_create(MPI_COMM_WORLD,num_dim,dimensions,periods,reorder,&processGrid);
            MPI_Comm_rank(processGrid,&comm_rank);

            int block=dimensions[0];
            int sqrt_num_processor=sqrt(num_processor);
            int block_size=n*n/(block*block);
            double *blockA=(double *)calloc(block_size,sizeof(double));
            double *blockB=(double *)calloc(block_size,sizeof(double));
            double *blockC=(double *)calloc(block_size,sizeof(double));
          
            double *tempA=(double *)calloc(block_size,sizeof(double));
            double *tempB=(double *)calloc(block_size,sizeof(double));

            double *MatrixA=(double *)calloc(n*n,sizeof(double));
            double *MatrixB=(double *)calloc(n*n,sizeof(double));
            double *MatrixC=(double *)calloc(n*n,sizeof(double));
          //  generate_matrix(MatrixC,n,seed);
          
            if(comm_rank==0){
              generate_matrix(MatrixA,n,seed);
              generate_matrix(MatrixB,n,seed);
              print_mat(MatrixA,n);
              print_mat(MatrixB,n);
              start = MPI_Wtime();
              send_block_AB(MatrixA,MatrixB,blockA,blockB,tempA,tempB,n,n/block,num_processor,sqrt_num_processor,processGrid);
            }
            int row_rank,col_rank;
            MPI_Cart_coords(processGrid,comm_rank,num_dim,coords);
            MPI_Recv(blockA,block_size,MPI_DOUBLE,0,1,processGrid,&stat);
            MPI_Recv(blockB,block_size,MPI_DOUBLE,0,2,processGrid,&stat);
            
            MPI_Comm_split(processGrid,coords[0],coords[1],&processRow);
            MPI_Comm_rank(processRow,&row_rank);
            MPI_Comm_split(processGrid,coords[1],coords[1],&processCol);
            MPI_Comm_rank(processCol,&col_rank);

            MPI_Request req1,req2;
            int times,m;
            for(times=0;times<block;times++){
              m=(coords[0]+times)%block;
              memcpy(tempA,blockA,sizeof(double)*block_size);
              MPI_Bcast(tempA,block_size,MPI_DOUBLE,m,processRow);
              MPI_Irecv(tempB,block_size,MPI_DOUBLE,col_rank==block-1?0:col_rank+1,3,processCol,&req1);
              MPI_Isend(blockB,block_size,MPI_DOUBLE,col_rank==0?block-1:col_rank-1,3,processCol,&req2);
              matmul(blockA,blockB,blockC,n);
              MPI_Wait(&req1,&stat);
              MPI_Wait(&req2,&stat);
              memcpy(blockB,tempB,sizeof(double)*block_size);
            }
            MPI_Isend(blockC,block_size,MPI_DOUBLE,0,4,processGrid,&request);
            if(comm_rank==0){
              for(int i=0;i<block*block;i++){
                MPI_Recv(MatrixC,block_size,MPI_DOUBLE,i,4,processGrid,&stat);
              }
              end=MPI_Wtime();
            }
            print_mat(MatrixC,n);
            free(blockA);
            free(blockB);
            free(blockC);
            MPI_Wait(&request, &stat);
            free(tempA);
            free(tempB);
            free(MatrixA);
            free(MatrixB);
            free(MatrixC);
            MPI_Finalize();
  
}


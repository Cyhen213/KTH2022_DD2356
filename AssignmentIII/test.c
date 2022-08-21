#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>
#include <string.h>

void generate_matrix(double *m,int n,int seed){
    int row,col;
  srand(seed);
    for(row=0;row<n;row++){
        for(col=0;col<n;col++){
            m[row*n+col]=((double)rand()*100/(double)(RAND_MAX));
        }
    }
}
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
void matmul(double *a,double *b,double *c,int n){
        for(int i=0;i<n;i++){
          for(int j=0;j<n;j++){
            for(int k=0;k<n;k++){
              c[n*i+j]+=a[n*i+k]*b[n*k+j];
            }
          }
      }
}
void send_block_AB(double *MatrixA,double *MatrixB,double *blockA,double *blockB,double *tempA,double *tempB,int n,int block,int num_processor,int sqrt_num_processor,MPI_Comm processGrid,MPI_Request req)
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
         // MPI_Send(tempA,block*block,MPI_DOUBLE,i,1,processGrid);
          MPI_Isend(tempA,block*block,MPI_DOUBLE,i,1,processGrid,&req);

       //   MPI_Send(tempB,block*block,MPI_DOUBLE,i,2,processGrid);
          MPI_Isend(tempB,block*block,MPI_DOUBLE,i,2,processGrid,&req);

        }
    }
}

int main(int argc, char *argv[])
{
  int num_processor,firstRank,provided;
  MPI_Init_thread(&argc, &argv, MPI_THREAD_SINGLE, &provided);
  MPI_Comm_size(MPI_COMM_WORLD,&num_processor);
  MPI_Comm_rank(MPI_COMM_WORLD, &firstRank);
  
  int n=9;
  int block=3;
 
  int sqrt_num_processor;
  

  int block_size=n*n/(block*block);
  int seed=1;
  int num_dim=2;
  int dimensions[2]={0,0};
  
  int periods[2]={1,1};
  int reorder=0;
  int coords[2];
  
  
  double start, end;

  
  MPI_Request req;
  MPI_Status stat;
  int comm_rank;
  MPI_Comm processGrid;
  MPI_Dims_create(num_processor,num_dim,dimensions);
  MPI_Cart_create(MPI_COMM_WORLD,num_dim,dimensions,periods,reorder,&processGrid);
  MPI_Comm_rank(processGrid,&comm_rank);
  MPI_Cart_shift(cartesian, 0, 1, &neighbour[0], &neighbour[1]);

  MPI_Cart_coords(processGrid,comm_rank,num_dim,coords);
  MPI_Cart_sub(cartesian, remain_dims, &sub); // create sub
  MPI_Comm_rank(sub, &subrank); // the rank in sub
  double *MatrixA=(double *)calloc(n*n,sizeof(double));
  double *MatrixB=(double *)calloc(n*n,sizeof(double));
  double *MatrixC=(double *)calloc(n*n,sizeof(double));

  sqrt_num_processor=sqrt(num_processor);
  double *blockA=(double *)calloc(block_size,sizeof(double));
  double *blockB=(double *)calloc(block_size,sizeof(double));
  double *blockC=(double *)calloc(block_size,sizeof(double));

  double *tempA=(double *)calloc(block_size,sizeof(double));
  double *tempB=(double *)calloc(block_size,sizeof(double));

  MPI_Datatype newtype;
  MPI_Type_vector(lblock, lblock, N, MPI_DOUBLE, &newtype);
  MPI_Type_commit(&newtype);

  if(comm_rank==0){
    generate_matrix(MatrixA,n,seed);
    generate_matrix(MatrixB,n,seed);
    start = MPI_Wtime();
    for (int i = 1; i < n; ++i) {
        MPI_Cart_coords(processGrid, i, num_dim, coords);
        MPI_Send(&MatrixA[(N * (N * coords[0] + coords[1])) / block], 1, newtype, i, 0, processGrid);
        MPI_Send(&MatrixB[(N * (N * coords[0] + coords[1])) / block], 1, newtype, i, 1, processGrid);
        MPI_Send(&MatrixC[(N * (N * coords[0] + coords[1])) / block], 1, newtype, i, 1, processGrid);
    }
    
    MPI_Cart_coords(cartesian, 0, ndims, coords);
    for (int i = 0; i < lblock; ++i) {
        for (int j = 0; j < lblock; ++j) {
            blockA[i * block + j] = MatrixA[((coords[0] * block) + i) * n + (coords[1] * block + j)];
            blockB[i * block + j] = MatrixB[((coords[0] * block) + i) * n + (coords[1] * block + j)];
            blockC[i * lblock + j] = MatrixC[((coords[0] * block) + i) * n + (coords[1] * block + j)];
        }
  
  }
  else if(comm_rank!=0){
    MPI_Recv(blockA, block * block, MPI_DOUBLE, 0, 0, processGrid, MPI_STATUS_IGNORE);
    MPI_Recv(blockB, block * block, MPI_DOUBLE, 0, 1, processGrid, MPI_STATUS_IGNORE);
    MPI_Recv(blockC, block * block, MPI_DOUBLE, 0, 1, processGrid, MPI_STATUS_IGNORE);
  }
  //MPI_Recv(blockA,block_size,MPI_DOUBLE,0,1,processGrid,&stat);
  //MPI_Recv(blockB,block_size,MPI_DOUBLE,0,2,processGrid,&stat);
  
  int row_rank,col_rank;
  MPI_Comm processRow,processCol;
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
  MPI_Isend(blockC,block_size,MPI_DOUBLE,0,4,processGrid,&req);
  if(comm_rank==0){
    for(int i=0;i<block*block;i++){
      MPI_Recv(MatrixC,block_size,MPI_DOUBLE,i,4,processGrid,&stat);
    }
    end=MPI_Wtime();
    print_mat(MatrixC,n);

  }
  free(blockA);
  free(blockB);
  free(blockC);
  MPI_Wait(&req, &stat);
  free(tempA);
  free(tempB);
  free(MatrixA);
  free(MatrixB);
  free(MatrixC);
  MPI_Finalize();
}


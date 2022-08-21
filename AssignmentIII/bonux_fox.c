#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>
#include <string.h>
#define SIZE 9
#include <time.h>



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
›
void normal_cal(double *a,double *b,double *c,double *MatrixC,int n)
{
  clock_t start, finish;
  double duration;

  start = clock();
  for(int i=0;i<n;i++){
    for(int j=0;j<n;j++){
      for(int k=0;k<n;k++){
        c[n*i+j]+=a[n*i+k]*b[n*k+j];
        if(c[n*i+j]-MatrixC[n*i+j]>1e-6 || MatrixC[n*i+j]-c[n*i+j]>1e-6 ){
          printf("wrong");
        }
      }
    }
  }
  printf("correct\n");
 // finish= clock();
  //duration = (double)(finish- start) / CLOCKS_PER_SEC;
//  return duration;
}
int main(int argc, char *argv[])
{
  int n;
  n=SIZE;
 
  double *A=(double *)malloc(n*n*sizeof(double));
  double *B=(double *)malloc(n*n*sizeof(double));
  double *C=(double *)malloc(n*n*sizeof(double));
  
  int num_processor,rank,provided;
  MPI_Init_thread(&argc, &argv, MPI_THREAD_SINGLE, &provided);
  MPI_Comm_size(MPI_COMM_WORLD,&num_processor);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  
  
  int block;
  block=sqrt(num_processor);
  int sqrt_num_processor;
  sqrt_num_processor=sqrt(num_processor);
  int block_len,block_size;
  block_len=n/block;
  block_size=block_len*block_len;
  
  int seed=1;
  
  int dimensions[2]={0,0};
  int remain_dims[2] = {0, 1};
  int periods[2]={1,1};
  int reorder=0;
  int coords[2];
  
  double start1, end1;
  MPI_Comm cartesian, sub;

  int newrank, subrank;
  int neighbour[2];
  int cart_rank;
  MPI_Dims_create(num_processor,2,dimensions);
  MPI_Cart_create(MPI_COMM_WORLD,2,dimensions,periods,reorder,&cartesian);
  
  MPI_Comm_rank(cartesian,&cart_rank);
  MPI_Cart_shift(cartesian, 0, 1, &neighbour[0], &neighbour[1]);
  
  MPI_Cart_coords(cartesian,cart_rank,2,coords);
  MPI_Cart_sub(cartesian, remain_dims, &sub); // create sub

  MPI_Comm_rank(sub, &subrank); // the rank in sub


  double *MatrixA=(double *)malloc(n*n*sizeof(double));
  double *MatrixB=(double *)malloc(n*n*sizeof(double));
  double *MatrixC=(double *)malloc(n*n*sizeof(double));


  double *blockA=(double *)malloc(block_size*sizeof(double));
  double *blockB=(double *)malloc(block_size*sizeof(double));
  

  double *tempA=(double *)malloc(block_size*sizeof(double));
  double *tempB=(double *)malloc(block_size*sizeof(double));
  
  double *blockC=(double *)calloc(block_size,sizeof(double));
  
  MPI_Datatype newtype;
  MPI_Type_vector(block_len, block_len, n, MPI_DOUBLE, &newtype);
  MPI_Type_commit(&newtype);

  if(cart_rank==0){
    generate_matrix(MatrixA,n,23);
    generate_matrix(MatrixB,n,23);

    start1 = MPI_Wtime();
    for (int i = 1; i < n; ++i) {
        MPI_Cart_coords(cartesian, i, 2, coords);
        MPI_Send(&MatrixA[(n * (n * coords[0] + coords[1])) / block_len], 1, newtype, i, 0, cartesian);
        MPI_Send(&MatrixB[(n * (n * coords[0] + coords[1])) / block_len], 1, newtype, i, 1, cartesian);
        MPI_Send(&MatrixC[(n * (n * coords[0] + coords[1])) / block_len], 1, newtype, i, 1, cartesian);
    }
    
    MPI_Cart_coords(cartesian, 0, 2, coords);
    for (int i = 0; i < block_len; ++i) {
        for (int j = 0; j < block_len; ++j) {
            blockA[i * block_len + j] = MatrixA[((coords[0] * block_len) + i) * n + (coords[1] * block + j)];
            blockB[i * block_len + j] = MatrixB[((coords[0] * block_len) + i) * n + (coords[1] * block + j)];
            blockC[i * block_len + j] = MatrixC[((coords[0] * block_len) + i) * n + (coords[1] * block + j)];
        }
    }
  }
 else {
    MPI_Recv(blockA, block_len * block_len, MPI_DOUBLE, 0, 0, cartesian, MPI_STATUS_IGNORE);
    MPI_Recv(blockB, block_len * block_len, MPI_DOUBLE, 0, 1, cartesian, MPI_STATUS_IGNORE);
    MPI_Recv(blockC, block_len * block_len, MPI_DOUBLE, 0, 1, cartesian, MPI_STATUS_IGNORE);
  }
  //MPI_Recv(blockA,block_size,MPI_DOUBLE,0,1,processGrid,&stat);
  //MPI_Recv(blockB,block_size,MPI_DOUBLE,0,2,processGrid,&stat);
  
  for (int i = 0; i < block_len; ++i)
  {
    
      memcpy(tempA, blockA, block_size * sizeof(double));
      int root = (coords[0] + i) % block_len;
      MPI_Bcast(tempA, block_size, MPI_DOUBLE, root, sub);
      matmul(tempA,blockB,blockC,block_len);
      MPI_Sendrecv(blockB, block_size, MPI_DOUBLE, neighbour[0], 0, tempB, block_size, MPI_DOUBLE, neighbour[1], MPI_ANY_TAG, cartesian, MPI_STATUS_IGNORE);
      memcpy(blockB, tempB, block_size * sizeof(double));
  }
  
  if (cart_rank != 0)
  {
      MPI_Send(blockC, block_size, MPI_DOUBLE, 0, 9, cartesian);
  }
  else
  {
    for (int i = 0; i < block_len; ++i) {
        for (int j = 0; j < block_len; ++j) {
              MatrixC[((coords[0] * block_len) + i) * n + (coords[1] * block_len + j)] = blockC[i * block_len + j];
          }
      }
      for (int i = 1; i < num_processor; ++i) {
          MPI_Cart_coords(cartesian, i, 2, coords);
          MPI_Recv(&MatrixC[(n * (n * coords[0] + coords[1])) / block_len], 1, newtype, i, 9, cartesian, MPI_STATUS_IGNORE);
      }
      end1 = MPI_Wtime();
    normal_cal(A,B,C,MatrixC,n);
      printf("Time of fox %f s\n", end1 - start1);
      
  }


  free(MatrixA);
  free(MatrixB);
  free(MatrixC);
  free(blockA);
  free(tempA);
  free(blockB);
  free(tempB);
  free(blockC);
  MPI_Finalize();
}


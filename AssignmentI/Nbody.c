#define DIM 2
typedef double vect_t[DIM];
#include <string.h> /* For memset */
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>
#define delta_t 0.05
#define N 500
#define n_steps 100
#define G 6.67408e-11
double *mass;
vect_t *forces,*vel,*pos;
double x_diff,y_diff,dist,dist_cubed;
double ts, t;

double mysecond(){
  struct timeval tp;
  struct timezone tzp;
  int i;

  i = gettimeofday(&tp,&tzp);
  return ( (double) tp.tv_sec + (double) tp.tv_usec * 1.e-6 );
}


int main(){
    forces = (vect_t *)malloc(N*sizeof(vect_t));
    mass=(double *)malloc(N*sizeof(double));
    vel=(vect_t *)malloc(N*sizeof(vect_t));
    pos=(vect_t *)malloc(N*sizeof(vect_t));
    for(int q=0;q<N;q++){
        pos[q][0]=(rand() / (double)(RAND_MAX)) * 2 - 1;
        pos[q][1] = (rand() / (double)(RAND_MAX)) * 2 - 1;
        
        vel[q][0] = (rand() / (double)(RAND_MAX)) * 2 - 1;
        vel[q][1] = (rand() / (double)(RAND_MAX)) * 2 - 1;
        
        mass[q] = fabs((rand() / (double)(RAND_MAX)) * 2 - 1);
    }
    for(int step=1;step<=n_steps;step++){
        forces=memset(forces,0,N*sizeof(vect_t));
        for (int q=0;q<N-1;q++){
            for(int k=0;k<N;k++){
                if(k!=q){
                    x_diff = pos[q][0] - pos[k][0];
                    y_diff = pos[q][1] - pos[k][1];
                    dist = sqrt(x_diff*x_diff + y_diff*y_diff);
                    dist_cubed= dist*dist*dist;
                    forces[q][0]-=G*mass[q]*mass[k]/dist_cubed * x_diff;
                    forces[q][1]-=forces[q][1]-G*mass[q]*mass[k]/dist_cubed *y_diff;
                }
            }
        }
    }
    for(int q=0;q<N;q++){
        pos[q][0] += delta_t*vel[q][0];
        pos[q][1] += delta_t*vel[q][1];
        vel[q][0] += delta_t/mass[q]*forces[q][0];
        vel[q][1] += delta_t/mass[q]*forces[q][1];
    }
    ts = mysecond();
    for(int step=1;step<=n_steps;step++){
        forces=memset(forces,0,N*sizeof(vect_t));
        for (int q=0;q<N-1;q++){
            for(int k=0;k<N;k++){
                if(k!=q){
                    x_diff = pos[q][0] - pos[k][0];
                    y_diff = pos[q][1] - pos[k][1];
                    dist = sqrt(x_diff*x_diff + y_diff*y_diff);
                    dist_cubed= dist*dist*dist;
                    forces[q][0]-=G*mass[q]*mass[k]/dist_cubed * x_diff;
                    forces[q][1]-=forces[q][1]-G*mass[q]*mass[k]/dist_cubed *y_diff;
                    }
                }
            }
        }
    for(int q=0;q<N-1;q++){
        pos[q][0] += delta_t*vel[q][0];
        pos[q][1] += delta_t*vel[q][1];
        vel[q][0] += delta_t/mass[q]*forces[q][0];
        vel[q][1] += delta_t/mass[q]*forces[q][1];
    }
    
    t = mysecond() - ts;
    printf("Time for Nbody problem when N=,%d, time : T = %f /s\n", N, t);
}


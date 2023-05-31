#include <cstdio>
#include <string>
#include <cstdlib>
#include <math.h>
#include <omp.h>
#include <sys/time.h>

#define NSPEEDS 9
#define FN1 "./png/iter"
#define FN2 ".dat"
#define FN3 ".png"
#define TARGETDAT " ./png/iter.dat"
#define TARGETPNG " ./png/iter.png"
#define NUMOFTHREADS 48

/* struct to hold the parameter values */
typedef struct{
  int    nx;            /* no. of cells in x-direction */
  int    ny;            /* no. of cells in y-direction */
  int    maxIters;      /* no. of iterations */
  int    reynolds_dim;  /* dimension for Reynolds number */
  float density;       /* density per link */
  float accel;         /* density redistribution */
  float omega;         /* relaxation parameter */
  int framerate;       /* num of frame per second*/
} t_param;

/* struct to hold the 'speed' values */
typedef struct{
  float speeds[NSPEEDS];
} t_speed;

void initialise(const char* paramfile, const char* obstaclefile, t_param* params, t_speed** cells_ptr, t_speed** tmp_cells_ptr, int** obstacles_ptr);
void timestep(const t_param params, t_speed* cells, t_speed* tmp_cells, int* obstacles, int tt);
void accelerate_flow(const t_param params, t_speed* cells, int* obstacles);
void propagate(const t_param params, t_speed* cells, t_speed* tmp_cells);
void rebound(const t_param params, t_speed* cells, t_speed* tmp_cells, int* obstacles);
void collision(const t_param params, t_speed* cells, t_speed* tmp_cells, int* obstacles);
void calc_values(const t_param params, t_speed* cells, int* obstacles, float* u, float* vorticity,float* local_density,float* u_x,float* u_y);
void write_values(const t_param params, float* u, float* vorticity, int tt);
void finalise(const t_param* params, t_speed* cells, t_speed* tmp_cells,
             int* obstacles, float* u, float* vorticity, float* local_density, float* u_x, float* u_y);
float av_velocity(const t_param params, t_speed* cells, int* obstacles);
float calc_reynolds(const t_param params, t_speed* cells, int* obstacles);
void usage(const char* exe);
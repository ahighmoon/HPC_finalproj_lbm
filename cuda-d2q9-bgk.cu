/*
** Code to implement a d2q9-bgk lattice boltzmann scheme.
** 'd2' inidates a 2-dimensional grid, and
** 'q9' indicates 9 velocities per grid cell.
** 'bgk' refers to the Bhatnagar-Gross-Krook collision step.
**
** The 'speeds' in each cell are numbered as follows:
**
** 6 2 5
**  \|/
** 3-0-1
**  /|\
** 7 4 8
**
** A 2D grid:
**
**           cols
**       --- --- ---
**      | D | E | F |
** rows  --- --- ---
**      | A | B | C |
**       --- --- ---
**
** 'unwrapped' in row major order to give a 1D array:
**
**  --- --- --- --- --- ---
** | A | B | C | D | E | F |
**  --- --- --- --- --- ---
**
** Grid indicies are:
**
**          ny
**          ^       cols(ii)
**          |  ----- ----- -----
**          | | ... | ... | etc |
**          |  ----- ----- -----
** rows(jj) | | 1,0 | 1,1 | 1,2 |
**          |  ----- ----- -----
**          | | 0,0 | 0,1 | 0,2 |
**          |  ----- ----- -----
**          ----------------------> nx
**
** Note the names of the input parameter and obstacle files
** are passed on the command line, e.g.:
**
**   ./d2q9-bgk input.params obstacles.dat
**
** Be sure to adjust the grid dimensions in the parameter file
** if you choose a different obstacle file.
*/

#include <cstdio>
#include <string>
#include <cstdlib>
#include <math.h>
#include <sys/time.h>

#define NSPEEDS         9
#define FN1 "./png/iter"
#define FN2 ".dat"
#define FN3 ".png"
#define TARGETDAT " ./png/iter.dat"
#define TARGETPNG " ./png/iter.png"
#define BLOCK_DIM 8

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

/*
** function prototypes
*/

/* load params, allocate memory, load obstacles & initialise fluid particle densities */
int initialise(const char* paramfile, const char* obstaclefile,
               t_param* params, t_speed** cells_ptr, t_speed** tmp_cells_ptr,
               int** obstacles_ptr);

/*
** The main calculation methods.
** timestep calls, in order, the functions:
** accelerate_flow(), propagate(), rebound() & collision()
*/
int accelerate_flow(const t_param params, t_speed* cells, int* obstacles);
__global__ void propagate(const t_param params, t_speed* cells, t_speed* tmp_cells);
__global__ void rebound(const t_param params, t_speed* cells, t_speed* tmp_cells, int* obstacles);
__global__ void collision(const t_param params, t_speed* cells, t_speed* tmp_cells, int* obstacles);
__global__ void calc_values(const t_param params, t_speed* cells, int* obstacles, float* u, float* vorticity,float* local_density,float* u_x,float* u_y);
int timestep(const t_param params, t_speed* cells, t_speed* tmp_cells, int* obstacles, int tt, dim3 gridDim, dim3 blockDim);
int write_values(const t_param params, float* u, float* vorticity, int tt);

/* compute average velocity */
float av_velocity(const t_param params, t_speed* cells, int* obstacles);

/* calculate Reynolds number */
float calc_reynolds(const t_param params, t_speed* cells, int* obstacles);

/* utility functions */
void usage(const char* exe);

int main(int argc, char* argv[])
{
  char*    paramfile = NULL;    /* name of the input parameter file */
  char*    obstaclefile = NULL; /* name of a the input obstacle file */
  t_param  params;              /* struct to hold parameter values */
  t_speed* cells     = NULL;    /* grid containing fluid densities */
  t_speed* tmp_cells = NULL;    /* scratch space */
  int*     obstacles = NULL;    /* grid indicating which cells are blocked */
  struct timeval timstr;                                                             /* structure to hold elapsed time */
  double compute_tic, compute_toc, write_time, write_tic, write_toc, total_tic, total_toc; /* elapsed time */
  write_time = 0;

  /* parse the command line */
  if (argc != 3){
    usage(argv[0]);
  }
  else{
    paramfile = argv[1];
    obstaclefile = argv[2];
  }
 
  /* initialise our data structures and load values from file */
  t_param*  paramsgpu;              /* struct to hold parameter values */
  t_speed* cellsgpu;
  t_speed* tmp_cellsgpu;
  int*     obstaclesgpu;
  float*  ugpu;
  float*  vorticitygpu;
  float*  local_densitygpu;
  float*  u_xgpu;
  float*  u_ygpu;

  initialise(paramfile, obstaclefile, &params, &cells, &tmp_cells, &obstacles);
  long memsz = params.ny * params.nx;
  float* u = (float*) malloc(sizeof(float) * memsz);
  float* vorticity = (float*) malloc(sizeof(float) * memsz);
  
  cudaMalloc(&paramsgpu, sizeof(t_param));
  cudaMalloc(&cellsgpu, sizeof(t_speed) * memsz);
  cudaMalloc(&tmp_cellsgpu, sizeof(t_speed) * memsz);
  cudaMalloc(&obstaclesgpu, sizeof(int) * memsz);
  cudaMalloc(&ugpu, sizeof(float) * memsz);
  cudaMalloc(&vorticitygpu, sizeof(float) * memsz);
  cudaMalloc(&local_densitygpu, sizeof(float) * memsz);
  cudaMalloc(&u_xgpu, sizeof(float) * memsz);
  cudaMalloc(&u_ygpu, sizeof(float) * memsz);
  
  cudaMemcpy(paramsgpu, &params, sizeof(t_param), cudaMemcpyHostToDevice);
  cudaMemcpy(cellsgpu, cells, sizeof(t_speed)*memsz, cudaMemcpyHostToDevice);
  cudaMemcpy(tmp_cellsgpu, tmp_cells, sizeof(t_speed)*memsz, cudaMemcpyHostToDevice);
  cudaMemcpy(obstaclesgpu, obstacles, sizeof(int)*memsz, cudaMemcpyHostToDevice);

  dim3 blockDim(BLOCK_DIM, BLOCK_DIM);
  dim3 gridDim(params.nx/BLOCK_DIM, params.ny/BLOCK_DIM);
  cudaDeviceSynchronize();

  /* Total/init time starts here: */
  gettimeofday(&timstr, NULL);
  total_tic = timstr.tv_sec + (timstr.tv_usec / 1000000.0);

  for (int tt = 0; tt < params.maxIters; tt++){
    timestep(params, cellsgpu, tmp_cellsgpu, obstaclesgpu, tt, gridDim, blockDim);
    /*
    if ((tt + 1) % params.framerate == 0){
      calc_values<<<gridDim,blockDim>>>(params, cellsgpu, obstaclesgpu, ugpu, vorticitygpu, local_densitygpu, u_xgpu, u_ygpu);
      gettimeofday(&timstr, NULL);
      write_tic = timstr.tv_sec + (timstr.tv_usec / 1000000.0);
      cudaMemcpy(u, ugpu, sizeof(float) * memsz, cudaMemcpyDeviceToHost);
      cudaMemcpy(vorticity, vorticitygpu, sizeof(float) * memsz, cudaMemcpyDeviceToHost);
      write_values(params, u, vorticity, 1000+(tt+1)/params.framerate);
      gettimeofday(&timstr, NULL);
      write_toc = timstr.tv_sec + (timstr.tv_usec / 1000000.0);
      write_time += write_toc - write_tic;
    }
    */
  }
  cudaDeviceSynchronize();
  gettimeofday(&timstr, NULL);
  total_toc = timstr.tv_sec + (timstr.tv_usec / 1000000.0);
  
  cudaMemcpy(cells, cellsgpu, sizeof(t_speed)*memsz, cudaMemcpyDeviceToHost);
  cudaMemcpy(obstacles, obstaclesgpu, sizeof(int)*memsz, cudaMemcpyDeviceToHost);

  printf("==done==\n");
  printf("Reynolds number:\t\t%.12E\n", calc_reynolds(params, cells, obstacles));
  printf("Total elapsed time:\t\t\t%.6lf (s)\n", total_toc - total_tic);
  printf("Total write time:\t\t\t%.6lf (s)\n", write_time);
  printf("Total compute time:\t\t\t%.6lf (s)\n", total_toc - total_tic - write_time);
  printf("Producing animation ...\n");

  /* make two graphs*/
  /*
  std::string FN;
  std::string command;
  std::string mv = "mv";
  for (int graph = 0; graph < 2; graph ++){
    for (int iter = 1001; iter <= 1000 + params.maxIters / params.framerate; iter++){
      FN = FN1 + std::to_string(iter) + FN2;
      command = mv + " " + FN + TARGETDAT;
      system(command.c_str());
      switch (graph){
      case 0:
        system("gnuplot velocity.plt");
        break;
      case 1:
        system("gnuplot vortex.plt");
        break;
      }
      command = mv + TARGETPNG + " " + FN1 + std::to_string(iter) + FN3;
      system(command.c_str());
      command = mv + TARGETDAT + " " + FN;
      system(command.c_str());
    }
    if (graph == 1) {
      system("cd png && convert -delay 4 -loop 1 *.png vortex.gif");
      system("cd png && rm ./*.png *.dat");
    }
    else{
      system("cd png && convert -delay 4 -loop 1 *.png velocity.gif");
    }
  }
  */
  cudaFree(paramsgpu);
  cudaFree(cellsgpu);
  cudaFree(tmp_cellsgpu);
  cudaFree(obstaclesgpu);
  cudaFree(ugpu);
  cudaFree(vorticitygpu);
  cudaFree(local_densitygpu);
  cudaFree(u_xgpu);
  cudaFree(u_ygpu);
  free(u);
  free(vorticity);
  free(cells);
  free(tmp_cells);
  free(obstacles);
  return EXIT_SUCCESS;
}

int timestep(const t_param paramsgpu, t_speed* cellsgpu, t_speed* tmp_cellsgpu, int* obstaclesgpu, int tt, dim3 gridDim, dim3 blockDim){
  propagate<<<gridDim,blockDim>>>(paramsgpu, cellsgpu, tmp_cellsgpu);
  //accelerate_flow(params, cells, obstacles);
  rebound<<<gridDim,blockDim>>>(paramsgpu, cellsgpu, tmp_cellsgpu, obstaclesgpu);
  collision<<<gridDim,blockDim>>>(paramsgpu, cellsgpu, tmp_cellsgpu, obstaclesgpu);
  return EXIT_SUCCESS;
}

__global__
void propagate(const t_param params, t_speed* cells, t_speed* tmp_cells){
  /* loop over _all_ cells */
  float w1 = params.density  / 9.f;
  float w2 = params.density / 36.f;

  long ii = threadIdx.x + blockDim.x * blockIdx.x;
  long jj = threadIdx.y + blockDim.y * blockIdx.y;
  //printf("From t %d,%d of b %d,%d. My ii and jj are %ld and %ld\n", threadIdx.x, threadIdx.y, blockIdx.x, blockIdx.y, ii, jj);

  if (ii < params.nx && jj < params.ny){
    int index = ii + jj*params.nx;
    /* determine indices of axis-direction neighbours
    ** respecting periodic boundary conditions (wrap around) */
    int y_n = (jj + 1) % params.ny;
    int x_e = (ii + 1) % params.nx;
    int y_s = (jj == 0) ? (params.ny - 1) : (jj - 1);
    int x_w = (ii == 0) ? (params.nx - 1) : (ii - 1);

    /* propagate densities from neighbouring cells, following
    ** appropriate directions of travel and writing into
    ** scratch space grid */
    tmp_cells[index].speeds[0] = cells[index].speeds[0]; /* central cell, no movement */
    tmp_cells[index].speeds[1] = cells[x_w + jj*params.nx].speeds[1]; /* east */
    tmp_cells[index].speeds[2] = cells[ii + y_s*params.nx].speeds[2]; /* north */
    tmp_cells[index].speeds[3] = cells[x_e + jj*params.nx].speeds[3]; /* west */
    tmp_cells[index].speeds[4] = cells[ii + y_n*params.nx].speeds[4]; /* south */
    tmp_cells[index].speeds[5] = cells[x_w + y_s*params.nx].speeds[5]; /* north-east */
    tmp_cells[index].speeds[6] = cells[x_e + y_s*params.nx].speeds[6]; /* north-west */
    tmp_cells[index].speeds[7] = cells[x_e + y_n*params.nx].speeds[7]; /* south-west */
    tmp_cells[index].speeds[8] = cells[x_w + y_n*params.nx].speeds[8]; /* south-east */
    
    if (ii==0){
      tmp_cells[index].speeds[1] = w1*3; /* east */
      tmp_cells[index].speeds[5] = w2; /* north-east */
      tmp_cells[index].speeds[8] = w2; /* south-east */
    }
    else if(ii==params.nx-1){
      tmp_cells[index].speeds[3] = w1; /* east */
      tmp_cells[index].speeds[6] = w2; /* north-east */
      tmp_cells[index].speeds[7] = w2; /* south-east */
    }
  }
}

__global__ 
void rebound(const t_param params, t_speed* cells, t_speed* tmp_cells, int* obstacles){
  /* loop over the cells in the grid */
  //#pragma omp parallel for num_threads(NUMOFTHREADS)
  
  long ii = threadIdx.x + blockDim.x * blockIdx.x;
  long jj = threadIdx.y + blockDim.y * blockIdx.y;
  int index = ii + jj*params.nx;
  /* if the cell contains an obstacle */
  if (obstacles[index]){
    /* called after propagate, so taking values from scratch space
    ** mirroring, and writing into main grid */
    cells[index].speeds[1] = tmp_cells[index].speeds[3];
    cells[index].speeds[2] = tmp_cells[index].speeds[4];
    cells[index].speeds[3] = tmp_cells[index].speeds[1];
    cells[index].speeds[4] = tmp_cells[index].speeds[2];
    cells[index].speeds[5] = tmp_cells[index].speeds[7];
    cells[index].speeds[6] = tmp_cells[index].speeds[8];
    cells[index].speeds[7] = tmp_cells[index].speeds[5];
    cells[index].speeds[8] = tmp_cells[index].speeds[6];
  }
}

__global__ 
void collision(const t_param params, t_speed* cells, t_speed* tmp_cells, int* obstacles){
  const float c_sq = 1.f / 3.f; /* square of speed of sound */
  const float w0 = 4.f / 9.f;  /* weighting factor */
  const float w1 = 1.f / 9.f;  /* weighting factor */
  const float w2 = 1.f / 36.f; /* weighting factor */

  /* loop over the cells in the grid
  ** NB the collision step is called after
  ** the propagate step and so values of interest
  ** are in the scratch-space grid */  
  long ii = threadIdx.x + blockDim.x * blockIdx.x;
  long jj = threadIdx.y + blockDim.y * blockIdx.y;
  int index = ii + jj*params.nx;
  /* don't consider occupied cells */
  if (!obstacles[index]){
    /* compute local density total */
    float local_density = 0.f;
    t_speed tmp_index = tmp_cells[index];

    for (int kk = 0; kk < NSPEEDS; kk++){
      local_density += tmp_index.speeds[kk];
    }
    /* compute x velocity component */
    float u_x = (tmp_index.speeds[1]
                  + tmp_index.speeds[5]
                  + tmp_index.speeds[8]
                  - (tmp_index.speeds[3]
                      + tmp_index.speeds[6]
                      + tmp_index.speeds[7]))
                  / local_density;
    /* compute y velocity component */
    float u_y = (tmp_index.speeds[2]
                  + tmp_index.speeds[5]
                  + tmp_index.speeds[6]
                  - (tmp_index.speeds[4]
                      + tmp_index.speeds[7]
                      + tmp_index.speeds[8]))
                  / local_density;

    /* velocity squared */
    float u_sq = u_x * u_x + u_y * u_y;

    /* directional velocity components */
    float u[NSPEEDS];
    u[1] =   u_x;        /* east */
    u[2] =         u_y;  /* north */
    u[3] = - u_x;        /* west */
    u[4] =       - u_y;  /* south */
    u[5] =   u_x + u_y;  /* north-east */
    u[6] = - u_x + u_y;  /* north-west */
    u[7] = - u_x - u_y;  /* south-west */
    u[8] =   u_x - u_y;  /* south-east */

    /* equilibrium densities */
    float d_equ[NSPEEDS];
    /* zero velocity density: weight w0 */
    d_equ[0] = w0 * local_density
                * (1.f - u_sq / (2.f * c_sq));
    /* axis speeds: weight w1 */
    float temp;
    for (int d = 1; d < NSPEEDS; d++){
      temp = local_density * (1.f + u[d] / c_sq
                                      + (u[d] * u[d]) / (2.f * c_sq * c_sq)
                                      - u_sq / (2.f * c_sq));
      if (d < 5) d_equ[d] = w1 * temp;
      else d_equ[d] = w2 * temp;
    }

    /* relaxation step */
    for (int kk = 0; kk < NSPEEDS; kk++){
      cells[index].speeds[kk] = tmp_index.speeds[kk]
                                              + params.omega
                                              * (d_equ[kk] - tmp_index.speeds[kk]);
    }
  }
}

float av_velocity(const t_param params, t_speed* cells, int* obstacles){
  int    tot_cells = 0;  /* no. of cells used in calculation */
  float tot_u;          /* accumulated magnitudes of velocity for each cell */

  /* initialise */
  tot_u = 0.f;

  /* loop over all non-blocked cells */
  for (int jj = 0; jj < params.ny; jj++){
    for (int ii = 0; ii < params.nx; ii++){
      int index = ii + jj*params.nx;
      /* ignore occupied cells */
      if (!obstacles[index]){
        /* local density total */
        float local_density = 0.f;

        for (int kk = 0; kk < NSPEEDS; kk++){
          local_density += cells[index].speeds[kk];
        }

        /* x-component of velocity */
        float u_x = (cells[index].speeds[1]
                      + cells[index].speeds[5]
                      + cells[index].speeds[8]
                      - (cells[index].speeds[3]
                         + cells[index].speeds[6]
                         + cells[index].speeds[7]))
                     / local_density;
        /* compute y velocity component */
        float u_y = (cells[index].speeds[2]
                      + cells[index].speeds[5]
                      + cells[index].speeds[6]
                      - (cells[index].speeds[4]
                         + cells[index].speeds[7]
                         + cells[index].speeds[8]))
                     / local_density;
        /* accumulate the norm of x- and y- velocity components */
        //printf("index = %d, local_density = %f, u_x = %f\n", index, local_density, u_x);
        tot_u += sqrtf((u_x * u_x) + (u_y * u_y));
        /* increase counter of inspected cells */
        ++tot_cells;
      }
    }
  }
 
  return tot_u / (float)tot_cells;
}

int initialise(const char* paramfile, const char* obstaclefile,
               t_param* params, t_speed** cells_ptr, t_speed** tmp_cells_ptr,
               int** obstacles_ptr){
  FILE*   fp;            /* file pointer */
  int    xx, yy;         /* generic array indices */
  int    blocked;        /* indicates whether a cell is blocked by an obstacle */

  /* open the parameter file */
  fp = fopen(paramfile, "r");

  fscanf(fp, "%d\n", &(params->nx));
  fscanf(fp, "%d\n", &(params->ny));
  fscanf(fp, "%d\n", &(params->maxIters));
  fscanf(fp, "%d\n", &(params->reynolds_dim));
  fscanf(fp, "%f\n", &(params->density));
  fscanf(fp, "%f\n", &(params->accel));
  fscanf(fp, "%f\n", &(params->omega));
  fscanf(fp, "%d\n", &(params->framerate));
  fclose(fp);

  /* main grid */
  *cells_ptr = (t_speed*)malloc(sizeof(t_speed) * (params->ny * params->nx));
  /* 'helper' grid, used as scratch space */
  *tmp_cells_ptr = (t_speed*)malloc(sizeof(t_speed) * (params->ny * params->nx));
  /* the map of obstacles */
  *obstacles_ptr = (int*)malloc(sizeof(int) * (params->ny * params->nx));
  /* initialise densities */
  float w0 = params->density * 4.f / 9.f;
  float w1 = params->density      / 9.f;
  float w2 = params->density      / 36.f;

  //#pragma omp parallel for num_threads(NUMOFTHREADS)
  for (int jj = 0; jj < params->ny; jj++){
    for (int ii = 0; ii < params->nx; ii++){
      /* center */
      (*cells_ptr)[ii + jj*params->nx].speeds[0] = w0;
      /* axis directions */
      (*cells_ptr)[ii + jj*params->nx].speeds[1] = w1;
      (*cells_ptr)[ii + jj*params->nx].speeds[2] = w1;
      (*cells_ptr)[ii + jj*params->nx].speeds[3] = w1;
      (*cells_ptr)[ii + jj*params->nx].speeds[4] = w1;
      /* diagonals */
      (*cells_ptr)[ii + jj*params->nx].speeds[5] = w2;
      (*cells_ptr)[ii + jj*params->nx].speeds[6] = w2;
      (*cells_ptr)[ii + jj*params->nx].speeds[7] = w2;
      (*cells_ptr)[ii + jj*params->nx].speeds[8] = w2;
    }
  }

  /* first set all cells in obstacle array to zero */
  for (int jj = 0; jj < params->ny; jj++){
    for (int ii = 0; ii < params->nx; ii++){
      (*obstacles_ptr)[ii + jj*params->nx] = 0;
    }
  }

  /* open the obstacle data file */
  fp = fopen(obstaclefile, "r");

  /* read-in the blocked cells list */
  while (fscanf(fp, "%d %d %d\n", &xx, &yy, &blocked) != EOF){
    /* assign to array */
    (*obstacles_ptr)[xx + yy*params->nx] = blocked;
  }
  fclose(fp);
  return EXIT_SUCCESS;
}

float calc_reynolds(const t_param params, t_speed* cells, int* obstacles){
  const float viscosity = 1.f / 6.f * (2.f / params.omega - 1.f);
  float avvel= av_velocity(params, cells, obstacles);
  //printf("avvel = %f\n", avvel);
  return avvel * params.reynolds_dim / viscosity;
}


__global__
void calc_values(const t_param params, t_speed* cells, int* obstacles, float* u, float* vorticity,float* local_density,float* u_x,float* u_y){
  long ii = threadIdx.x + blockDim.x * blockIdx.x;
  long jj = threadIdx.y + blockDim.y * blockIdx.y;
  
  if (ii < params.nx && jj < params.ny){
    int index = ii + jj*params.nx;

    if (obstacles[index]){
      u_x[index] = u_y[index] = 0.f;
    }
    else{
      local_density[index] = 0.f;
      for (int kk = 0; kk < NSPEEDS; kk++){
        local_density[index] += cells[index].speeds[kk];
      }
      u_x[index] = (cells[index].speeds[1]+ cells[index].speeds[5]+ cells[index].speeds[8]- (cells[index].speeds[3]
                  + cells[index].speeds[6]+ cells[index].speeds[7]))/ local_density[index];
      u_y[index] = (cells[index].speeds[2]+ cells[index].speeds[5]+ cells[index].speeds[6] - (cells[index].speeds[4]
              + cells[index].speeds[7]+ cells[index].speeds[8]))/ local_density[index];
      u[index] = sqrtf(u_x[index]*u_x[index]+u_y[index]*u_y[index]);
    }
  }

  //calculate the vorticity
  if (ii < params.nx && jj < params.ny){
    int index = ii + jj*params.nx;
    if (ii==0 || jj==0 || ii==params.nx-1 || jj==params.ny-1){
      vorticity[index] = 0;
    }
    else{
      double duy_dx = (u_y[ii+(jj+1)*params.nx] - u_y[ii+(jj-1)*params.nx]) / 2;
      double dux_dy = (u_x[ii+1+jj*params.nx] - u_x[ii-1+jj*params.nx]) / 2;
      vorticity[index] = duy_dx - dux_dy;
    }
  }
}

int write_values(const t_param params, float* u, float* vorticity, int tt){
  FILE* fp;
  std::string FN = FN1 + std::to_string(tt) + FN2;
  fp = fopen(FN.c_str(), "w");
  for (int jj = 0; jj < params.ny; ++jj) {
    for (int ii = 0; ii < params.nx; ++ii) {
      int index = ii + jj * params.nx;
      fprintf(fp, "%d %d %.12E %.12E\n", ii, jj, u[index], vorticity[index]);
    }
  }
  fclose(fp);
  return EXIT_SUCCESS;
}

void usage(const char* exe){
  fprintf(stderr, "Usage: %s <paramfile> <obstaclefile>\n", exe);
  exit(EXIT_FAILURE);
}

/******************************************************************************
* FILE: omp_bug4.c
* DESCRIPTION:
*   This very simple program causes a segmentation fault.
* AUTHOR: Blaise Barney  01/09/04
* LAST REVISED: 04/06/05
******************************************************************************/
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "utils.h"
#define N 1048

double** malloc_matrix(long n_row, long n_col){
  double** p = (double**)aligned_malloc(n_row * sizeof(double));
  for(long i = 0; i < n_row; i++)
    p[i] = (double*)aligned_malloc(n_col * sizeof(double));
  return p;
}

int main (int argc, char *argv[]) 
{
int nthreads, tid, i, j;

/* Fork a team of threads with explicit variable scoping */
#pragma omp parallel shared(nthreads) private(i,j,tid)
  {

  /* Obtain/print thread info */
  tid = omp_get_thread_num();
  if (tid == 0) 
    {
    nthreads = omp_get_num_threads();
    printf("Number of threads = %d\n", nthreads);
    }
  printf("Thread %d starting...\n", tid);

  double **a = malloc_matrix(N, N);
  /* Each thread works on its own private copy of the array */
  for (i=0; i<N; i++)
    for (j=0; j<N; j++)
      a[i][j] = tid + i + j;

  /* For confirmation */
  printf("Thread %d done. Last element= %f\n",tid,a[N-1][N-1]);

  aligned_free(a);

  }  /* All threads join master thread and disband */

}


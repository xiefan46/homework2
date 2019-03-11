// g++ -fopenmp -O3 -march=native MMult1.cpp && ./a.out

#include <stdio.h>
#include <math.h>
#include <omp.h>
#include "utils.h"

#define BLOCK_SIZE 12


// Note: matrices are stored in column major order; i.e. the array elements in
// the (m x n) matrix C are stored in the sequence: {C_00, C_10, ..., C_m0,
// C_01, C_11, ..., C_m1, C_02, ..., C_0n, C_1n, ..., C_mn}
void MMult0(long M, long N, long K, double **a, double **b, double **c) {
    /*for (long j = 0; j < n; j++) {
      for (long p = 0; p < k; p++) {
        for (long i = 0; i < m; i++) {
          double A_ip = a[i+p*m];
          double B_pj = b[p+j*k];
          double C_ij = c[i+j*m];
          C_ij = C_ij + A_ip * B_pj;
          c[i+j*m] = C_ij;
        }
      }
    }*/
    for (long i = 0; i < M; i++) {
        for (long j = 0; j < N; j++) {
            for (long k = 0; k < K; k++) {
                c[i][j] += a[i][k] * b[k][j];
            }
        }
    }
}

void MMult1(long M, long N, long K, double **a, double **b, double **c) {
    long m_block = M / BLOCK_SIZE;
    long n_block = N / BLOCK_SIZE;
    long k_block = K / BLOCK_SIZE;
    for (long ib = 0; ib < m_block; ib++) {
        for (long jb = 0; jb < n_block; jb++) {
            for (long kb = 0; kb < k_block; kb++) {
                //multi small block
                for (long i = 0; i < BLOCK_SIZE; i++) {
                    for (long j = 0; j < BLOCK_SIZE; j++) {
                        long ii = ib * BLOCK_SIZE + i;
                        long jj = jb * BLOCK_SIZE + j;
                        for (long k = 0; k < BLOCK_SIZE; k++) {
                            long kk = kb * BLOCK_SIZE + k;
                            c[ii][jj] += a[ii][kk] * b[kk][jj];
                        }
                    }
                }
            }
        }
    }
}

void MMult1_with_openmp(long M, long N, long K, double **a, double **b, double **c) {
    long m_block = M / BLOCK_SIZE;
    long n_block = N / BLOCK_SIZE;
    long k_block = K / BLOCK_SIZE;

    #pragma omp parallel for
    for (long ib = 0; ib < m_block; ib++) {
        for (long jb = 0; jb < n_block; jb++) {
            for (long kb = 0; kb < k_block; kb++) {
                //multi small block
                for (long i = 0; i < BLOCK_SIZE; i++) {
                    for (long j = 0; j < BLOCK_SIZE; j++) {
                        long ii = ib * BLOCK_SIZE + i;
                        long jj = jb * BLOCK_SIZE + j;
                        for (long k = 0; k < BLOCK_SIZE; k++) {
                            long kk = kb * BLOCK_SIZE + k;
                            c[ii][jj] += a[ii][kk] * b[kk][jj];
                        }
                    }
                }
            }
        }
    }
}

double** malloc_matrix(long n_row, long n_col){
    double** p = (double**)aligned_malloc(n_row * sizeof(double));
    for(long i = 0; i < n_row; i++)
        p[i] = (double*)aligned_malloc(n_col * sizeof(double));
    return p;
}

void init_matrix_rand(double **m, long n_row, long n_col){
    for(long i = 0; i < n_row; i++)
        for(long j = 0; j < n_col; j++)
            m[i][j] = drand48();
}

void init_matrix_zero(double **m, long n_row, long n_col){
    for(long i = 0; i < n_row; i++)
        for(long j = 0; j < n_col; j++)
            m[i][j] = 0;
}


int main(int argc, char **argv) {
    const long PFIRST = BLOCK_SIZE;
    const long PLAST = 2000;
    const long PINC = std::max(50 / BLOCK_SIZE, 1) * BLOCK_SIZE; // multiple of BLOCK_SIZE

    printf(" Dimension       Time    Gflop/s           GB/s             Time in the ref solution        Error\n");
    for (long p = PFIRST; p < PLAST; p += PINC) {
        long m = p, n = p, k = p;
        long NREPEATS = 1e9 / (m * n * k) + 1;
        /*double *a = (double *) aligned_malloc(m * k * sizeof(double)); // m x k
        double *b = (double *) aligned_malloc(k * n * sizeof(double)); // k x n
        double *c = (double *) aligned_malloc(m * n * sizeof(double)); // m x n
        double *c_ref = (double *) aligned_malloc(m * n * sizeof(double)); // m x n*/

        double **a = malloc_matrix(m, k);
        double **b = malloc_matrix(k, n);
        double **c = malloc_matrix(m, n);
        double **c_ref = malloc_matrix(m, n);

        // Initialize matrices
        init_matrix_rand(a, m,k);
        init_matrix_rand(b, k, n);
        init_matrix_zero(c, m, n);
        init_matrix_zero(c_ref, m, n);

        Timer t2;
        t2.tic();
        for (long rep = 0; rep < NREPEATS; rep++) { // Compute reference solution
            MMult0(m, n, k, a, b, c_ref);
        }
        double time_ref = t2.toc();

        Timer t;
        t.tic();
        for (long rep = 0; rep < NREPEATS; rep++) {
            MMult1_with_openmp(m, n, k, a, b, c);
        }
        double time = t.toc();
        double flops = 10 * m * n * k * NREPEATS / (time * 1e9); // TODO: calculate from m, n, k, NREPEATS, time
        double bandwidth = 4 * 8 * m * n * k * NREPEATS / (time * 1e9); // TODO: calculate from m, n, k, NREPEATS, time
        printf("%10d %10f %10f %10f        %10f", p, time, flops, bandwidth, time_ref);

        double max_err = 0;
        for (long i = 0; i < m ; i++){
            for(long j = 0; j < n; j++){
                max_err = std::max(max_err, fabs(c[i][j] - c_ref[i][j]));
            }
        }

        printf("                 %10e\n", max_err);

        aligned_free(a);
        aligned_free(b);
        aligned_free(c);
        aligned_free(c_ref);
    }

    return 0;
}

// * Using MMult0 as a reference, implement MMult1 and try to rearrange loops to
// maximize performance. Measure performance for different loop arrangements and
// try to reason why you get the best performance for a particular order?
//
//
// * You will notice that the performance degrades for larger matrix sizes that
// do not fit in the cache. To improve the performance for larger matrices,
// implement a one level blocking scheme by using BLOCK_SIZE macro as the block
// size. By partitioning big matrices into smaller blocks that fit in the cache
// and multiplying these blocks together at a time, we can reduce the number of
// accesses to main memory. This resolves the main memory bandwidth bottleneck
// for large matrices and improves performance.
//
// NOTE: You can assume that the matrix dimensions are multiples of BLOCK_SIZE.
//
//
// * Experiment with different values for BLOCK_SIZE (use multiples of 4) and
// measure performance.  What is the optimal value for BLOCK_SIZE?
//
//
// * Now parallelize your matrix-matrix multiplication code using OpenMP.
//
//
// * What percentage of the peak FLOP-rate do you achieve with your code?
//
//
// NOTE: Compile your code using the flag -march=native. This tells the compiler
// to generate the best output using the instruction set supported by your CPU
// architecture. Also, try using either of -O2 or -O3 optimization level flags.
// Be aware that -O2 can sometimes generate better output than using -O3 for
// programmer optimized code.

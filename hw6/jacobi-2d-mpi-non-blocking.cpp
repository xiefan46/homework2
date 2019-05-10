/* MPI-parallel Jacobi smoothing to solve -u''=f
 * Global vector has N unknowns, each processor works with its
 * part, which has lN = N/p unknowns.
 * Author: Georg Stadler
 */
#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include <string.h>
#include <stdlib.h>

int two2one(int i, int j, int lN) {
    return i * lN + j;
}

/* compuate global residual, assuming ghost values are updated */
double compute_residual(double *lu, int lN, double invhsq) {
    double tmp, gres = 0.0, lres = 0.0;

    /*for (i = 1; i <= lN; i++){
      tmp = ((2.0*lu[i] - lu[i-1] - lu[i+1]) * invhsq - 1);
      lres += tmp * tmp;*/

    for (int i = 1; i <= lN; i++) {
        for (int j = 1; j <= lN; j++) {
            int mid = two2one(i, j, lN);
            int top = two2one(i + 1, j, lN);
            int down = two2one(i - 1, j, lN);
            int left = two2one(i, j - 1, lN);
            int right = two2one(i, j + 1, lN);
            tmp = ((4.0 * lu[mid] - lu[top] - lu[down] - lu[left] -
                    lu[right]) * invhsq - 1);
            lres += tmp * tmp;
        }
    }

    /* use allreduce for convenience; a reduce would also be sufficient */
    MPI_Allreduce(&lres, &gres, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return sqrt(gres);
}


int main(int argc, char *argv[]) {
    //printf("start\n");
    int mpirank, i, p, N, lN, iter, max_iters;
    MPI_Status status, status1;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpirank);
    MPI_Comm_size(MPI_COMM_WORLD, &p);

    /* get name of host running MPI process */
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int name_len;
    MPI_Get_processor_name(processor_name, &name_len);
    printf("Rank %d/%d running on %s.\n", mpirank, p, processor_name);

    sscanf(argv[1], "%d", &N);
    sscanf(argv[2], "%d", &max_iters);

    MPI_Request request_out1, request_in1;
    MPI_Request request_out2, request_in2;
    MPI_Request request_out3, request_in3;
    MPI_Request request_out4, request_in4;

    /* compute number of unknowns handled by each process */
    lN = (int) sqrt(N * N / p);

    if ((N * N % p != 0) && mpirank == 0) {
        printf("N: %d, local N: %d\n", N, lN);
        printf("Exiting. N ^ 2 must be a multiple of p\n");
        MPI_Abort(MPI_COMM_WORLD, 0);
    }
    //printf("read parameters\n");
    /* timing */
    MPI_Barrier(MPI_COMM_WORLD);
    double tt = MPI_Wtime();

    /* Allocation of vectors, including left/upper and right/lower ghost points */
    double *lu = (double *) calloc(sizeof(double), (lN + 2) * (lN + 2));
    double *lunew = (double *) calloc(sizeof(double), (lN + 2) * (lN + 2));
    double *left_edge_send = (double *) calloc(sizeof(double), lN);
    double *left_edge_recv = (double *) calloc(sizeof(double), lN + 2);
    double *right_edge_send = (double *) calloc(sizeof(double), lN + 2);
    double *right_edge_recv = (double *) calloc(sizeof(double), lN + 2);

    double *lutemp;

    double h = 1.0 / (N + 1);
    double hsq = h * h;
    double invhsq = 1. / hsq;
    double gres, gres0, tol = 1e-5;

    //printf("calculate the location\n");
    //calculate the location of this mpi process
    int num_process_per_edge = (int) sqrt(p);
    int mpi_i_index = mpirank / num_process_per_edge;
    int mpi_j_index = mpirank % num_process_per_edge;

    /* initial residual */
    gres0 = compute_residual(lu, lN, invhsq);
    gres = gres0;

    for (iter = 0; iter < max_iters && gres / gres0 > tol; iter++) {

        /* Jacobi step for local points */
        /*for (i = 1; i <= lN; i++){
          lunew[i]  = 0.5 * (hsq + lu[i - 1] + lu[i + 1]);
        }*/
        for (int i = 1; i <= lN; i++) {
            for (int j = 1; j <= lN; j++) {
                int mid = two2one(i, j, lN);
                int top = two2one(i + 1, j, lN);
                int down = two2one(i - 1, j, lN);
                int left = two2one(i, j - 1, lN);
                int right = two2one(i, j + 1, lN);
                lunew[mid] = 0.25 * (hsq + lu[top] + lu[down] + lu[left] + lu[right]);
            }
        }

        /* communicate ghost values */


        /*if (mpirank < p - 1) {
          //If not the last process, send/recv bdry values to the right
          MPI_Send(&(lunew[lN]), 1, MPI_DOUBLE, mpirank+1, 124, MPI_COMM_WORLD);
          MPI_Recv(&(lunew[lN+1]), 1, MPI_DOUBLE, mpirank+1, 123, MPI_COMM_WORLD, &status);
        }
        if (mpirank > 0) {
          //If not the first process, send/recv bdry values to the left
          MPI_Send(&(lunew[1]), 1, MPI_DOUBLE, mpirank-1, 123, MPI_COMM_WORLD);
          MPI_Recv(&(lunew[0]), 1, MPI_DOUBLE, mpirank-1, 124, MPI_COMM_WORLD, &status1);
        }*/
        int topLeft = two2one(1, 1, lN);
        int topRight = two2one(1, lN, lN);
        int downLeft = two2one(lN, 1, lN);
        int downRight = two2one(lN, lN, lN);
        //printf("calculate edge\n");
        //has top process
        //MPI_Irecv(&(lunew[lN+1]), 1, MPI_DOUBLE, mpirank+1, 123, MPI_COMM_WORLD, &request_in1);
        //MPI_Isend(&(lunew[lN]), 1, MPI_DOUBLE, mpirank+1, 124, MPI_COMM_WORLD, &request_out1);
        if (mpi_i_index - 1 >= 0) {
            int top_process_rank = mpirank - num_process_per_edge;
            MPI_Isend(&(lunew[topLeft]), lN, MPI_DOUBLE, top_process_rank, 123, MPI_COMM_WORLD, &request_in1);
            MPI_Irecv(&(lunew[topLeft]), lN, MPI_DOUBLE, top_process_rank, 124, MPI_COMM_WORLD, &request_out1);
        }
        //printf("finish first communication\n");
        //has down process
        if (mpi_i_index + 1 < num_process_per_edge) {
            int down_process_rank = mpirank + num_process_per_edge;
            MPI_Isend(&(lunew[downLeft]), lN, MPI_DOUBLE, down_process_rank, 124, MPI_COMM_WORLD, &request_in2);
            MPI_Irecv(&(lunew[downLeft]), lN, MPI_DOUBLE, down_process_rank, 123, MPI_COMM_WORLD, &request_out2);
        }
        //printf("finish second communication\n");
        //has left process
        if (mpi_j_index - 1 >= 0) {

            for (i = 1; i <= lN; i++)
                left_edge_send[i] = lunew[two2one(i, 1, lN)];
            MPI_Isend(left_edge_send, lN, MPI_DOUBLE, mpirank - 1, 125, MPI_COMM_WORLD, &request_in3);
            MPI_Irecv(left_edge_recv, lN, MPI_DOUBLE, mpirank - 1, 126, MPI_COMM_WORLD, &request_out3);
            for (i = 1; i <= lN; i++)
                lunew[two2one(i, 1, lN)] = left_edge_recv[i];
        }
        //printf("finish third communication\n");
        //has right process
        if (mpi_j_index + 1 < num_process_per_edge) {
            for (i = 1; i <= lN; i++)
                right_edge_send[i] = lunew[two2one(i, lN, lN)];
            MPI_Isend(right_edge_send, lN, MPI_DOUBLE, mpirank + 1, 126, MPI_COMM_WORLD, &request_in4);
            MPI_Irecv(right_edge_recv, lN, MPI_DOUBLE, mpirank + 1, 125, MPI_COMM_WORLD, &request_out4);
            for (i = 1; i <= lN; i++)
                lunew[two2one(i, lN, lN)] = right_edge_recv[i];
        }
        //printf("finish forth communication\n");


        /* copy newu to u using pointer flipping */
        lutemp = lu;
        lu = lunew;
        lunew = lutemp;
        if (0 == (iter % 10)) {
            gres = compute_residual(lu, lN, invhsq);
            if (0 == mpirank) {
                printf("Iter %d: Residual: %g\n", iter, gres);
            }
        }
    }

    /* Clean up */
    free(lu);
    free(lunew);
    free(left_edge_recv);
    free(left_edge_send);
    free(right_edge_recv);
    free(right_edge_send);

    /* timing */
    MPI_Barrier(MPI_COMM_WORLD);
    double elapsed = MPI_Wtime() - tt;
    if (0 == mpirank) {
        printf("Time elapsed is %f seconds.\n", elapsed);
    }
    MPI_Finalize();
    return 0;
}

#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>

int getExpectedNumber(int initial_number, int num_proc, int N) {
    return initial_number + (num_proc - 1) * num_proc / 2 * N;
}

/*
 * Send an integer in ring order
 */
void ring_integer(int argc, char *argv[]) {
    const int TAG = 0;
    double total_time = 0.0;
    int rank, size;
    int N = atoi(argv[1]);
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    //printf( "I am %d of %d\n", rank, size );
    int number;
    int initial_number = 0;
    total_time -= MPI_Wtime();
    for (int i = 0; i < N; i++) {
        if (rank == 0) {
            if (i == 0)
                number = initial_number;
            MPI_Send(&number, 1, MPI_INT, (rank + 1) % size, TAG, MPI_COMM_WORLD);
            MPI_Recv(&number, 1, MPI_INT, size - 1, TAG, MPI_COMM_WORLD,
                     MPI_STATUS_IGNORE);
        } else {
            MPI_Recv(&number, 1, MPI_INT, rank - 1, TAG, MPI_COMM_WORLD,
                     MPI_STATUS_IGNORE);
            number += rank;
            MPI_Send(&number, 1, MPI_INT, (rank + 1) % size, TAG, MPI_COMM_WORLD);
        }

    }
    total_time += MPI_Wtime();
    if (rank == 0) {
        printf("Expected integer is %d\n", getExpectedNumber(initial_number, size, N));
        printf("Final integer is %d. Total time is : %lf\n", number, total_time);
    }
    MPI_Finalize();
}

/*
 * Communicate a large array of about 2MByte in a ring
 */
void communicate_array(int argc, char *argv[]) {
    const int TAG = 0;
    double total_time = 0.0;
    int rank, size;
    int N = atoi(argv[1]);
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    const int array_size = 1024;

    int* array2 = NULL;
    int* array = NULL;
    //int* array2 = (int*)malloc(sizeof(int) * array_size);

    total_time -= MPI_Wtime();
    for (int i = 0; i < N; i++) {

        if (rank == 0) {
            if(i == 0){
                //init array
                array = (int*)malloc(sizeof(int) * array_size);
                for(int i = 0; i < array_size; i++)
                    array[i] = i;

                MPI_Send(&array, array_size, MPI_INT, 1, TAG, MPI_COMM_WORLD);

            }
            else
                MPI_Send(&array2, array_size, MPI_INT, 1, TAG, MPI_COMM_WORLD);
            MPI_Recv(&array2, array_size, MPI_INT, size - 1, TAG, MPI_COMM_WORLD,
                     MPI_STATUS_IGNORE);
        } else {

            MPI_Recv(&array2, array_size, MPI_INT, rank - 1, TAG, MPI_COMM_WORLD,
                     MPI_STATUS_IGNORE);

            MPI_Send(&array2, array_size, MPI_INT, (rank + 1) % size, TAG, MPI_COMM_WORLD);

        }

    }
    total_time += MPI_Wtime();
    if (rank == 0) {
        printf(" Total time is : %lf\n", total_time);
        free(array);
    }
    //free(array2);
    MPI_Finalize();
}


int main(int argc, char *argv[]) {
    if (argc < 2) {
        printf("Expected usage : int_ring ${N rounds}\n");
        exit(1);
    }
    ring_integer(argc, argv);
    return 0;
}

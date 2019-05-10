// Parallel sample sort
#include <stdio.h>
#include <unistd.h>
#include <mpi.h>
#include <stdlib.h>
#include <algorithm>

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);

    int rank, p;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &p);

    // Number of random numbers per processor (this should be increased
    // for actual tests or could be passed in through the command line
    int N = 100;

    int *vec = (int *) malloc(N * sizeof(int));
    // seed random number generator differently on every core
    srand((unsigned int) (rank + 393919));

    // fill vector with random integers
    for (int i = 0; i < N; ++i) {
        vec[i] = rand();
    }
    printf("rank: %d, first entry: %d\n", rank, vec[0]);

    // sort locally
    std::sort(vec, vec + N);
    printf("local sort success\n");
    // sample p-1 entries from vector as the local splitters, i.e.,
    // every N/P-th entry of the sorted vector
    int *local_sample = (int *) malloc(sizeof(int) * (p - 1));
    int step = N / p;
    for (int i = 0; i < p - 1; i++) {
        local_sample[i] = vec[i * step];
    }
    printf("local sampling success\n");
    // every process communicates the selected entries to the root
    // process; use for instance an MPI_Gather
    int *gather_samples;
    int n_samples = p * (p - 1);
    if (rank == 0)
        gather_samples = (int *) malloc(sizeof(int) * n_samples);
    MPI_Gather(local_sample, p - 1, MPI_INT, gather_samples, p - 1, MPI_INT, 0,
               MPI_COMM_WORLD);

    // root process does a sort and picks (p-1) splitters (from the
    // p(p-1) received elements)
    int *split = (int *) malloc(sizeof(int) * (p - 1));
    if (rank == 0) {
        printf("gather all samples successfully\n");
        std::sort(gather_samples, gather_samples + n_samples);
        for (int i = 0; i < p - 1; i++)
            split[i] = gather_samples[i * p];
    }

    // root process broadcasts splitters to all other processes
    MPI_Bcast(split, p - 1, MPI_INT, 0, MPI_COMM_WORLD);
    printf("bcast all samples successfully\n");

    // every process uses the obtained splitters to decide which
    // integers need to be sent to which other process (local bins).
    // Note that the vector is already locally sorted and so are the
    // splitters; therefore, we can use std::lower_bound function to
    // determine the bins efficiently.
    //
    // Hint: the MPI_Alltoallv exchange in the next step requires
    // send-counts and send-displacements to each process. Determining the
    // bins for an already sorted array just means to determine these
    // counts and displacements. For a splitter s[i], the corresponding
    // send-displacement for the message to process (i+1) is then given by,
    // sdispls[i+1] = std::lower_bound(vec, vec+N, s[i]) - vec;
    int *send_counts = (int *) malloc(sizeof(int) * p);
    int *recv_counts = (int *) malloc(sizeof(int) * p);
    int *sdispls = (int *) malloc(sizeof(int) * p);
    int *rdispls = (int *) malloc(sizeof(int) * p);

    sdispls[0] = 0;
    for (int i = 0; i < p - 1; i++)
        sdispls[i + 1] = std::lower_bound(vec, vec + N, split[i]) - vec;
    // send and receive: first use an MPI_Alltoall to share with every
    // process how many integers it should expect, and then use
    // MPI_Alltoallv to exchange the data
    for (int i = 0; i < p - 1; i++)
        send_counts[i] = sdispls[i + 1] - sdispls[i];
    send_counts[p - 1] = N - sdispls[p - 2];
    MPI_Alltoall(send_counts, 1, MPI_INT, recv_counts, 1, MPI_INT, MPI_COMM_WORLD);
    printf("all to all\n");
    rdispls[0] = 0;
    for (int i = 0; i < p - 1; i++)
        rdispls[i + 1] = recv_counts[i] + rdispls[i];

    int ssize = 0;
    int rsize = 0;
    for (int i = 0; i < p; i++) {
        ssize += send_counts[i];
        rsize += recv_counts[i];
    }
    int *recv_vec = (int*)malloc(sizeof(int) * rsize);
    MPI_Alltoallv(vec, send_counts, sdispls, MPI_INT,
                  recv_vec, recv_counts, rdispls, MPI_INT,
                  MPI_COMM_WORLD);
    printf("all to all v\n");
    // do a local sort of the received data
    std::sort(recv_vec, recv_vec + recv_counts[rank]);

    // every process writes its result to a file
    FILE *fd = NULL;
    char filename[256];
    snprintf(filename, 256, "output%02d.txt", rank);
    fd = fopen(filename, "w+");

    if (NULL == fd) {
        printf("Error opening file \n");
        return 1;
    }

    for (int i = 0; i < rsize; i++)
        fprintf(fd, "  %d\n", recv_vec[i]);

    fclose(fd);
    printf("process %d write result to output file successfully\n", rank);
    free(vec);
    free(local_sample);
    if (rank == 0)
        free(gather_samples);
    free(split);
    free(sdispls);
    MPI_Finalize();
    return 0;
}

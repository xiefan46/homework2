#include <iostream>
#include <cstdlib>
#include "mpi.h"
#include <limits.h>

using namespace std;

const int MAX_ITER_ROUNDS = 10000;

struct Point {
    int x, y;
};

Point *createRandomPoints(int n) {
    Point *points = (Point *) malloc(sizeof(Point) * n);
    for (int i = 0; i < n; i++) {
        points[i].x = rand();
        points[i].y = rand();
    }
    return points;
}

Point *createRandomCenters(Point *points, int k, int n) {
    int step = n / k;
    Point *centers = (Point *) malloc(sizeof(Point) * k);
    for (int i = 0; i < k; i++) {
        int index = i * step + rand() % step;
        centers[i].x = points[index].x;
        centers[i].y = points[index].y;
    }
    return centers;
}

void toBuffer(Point *points, int *buffer, int n) {
    for (int i = 0; i < n; i++) {
        buffer[2 * i] = points[i].x;
        buffer[2 * i + 1] = points[i].y;
    }
}

void toPoints(int *buffer, Point *points, int n) {
    for (int i = 0; i < n; i++) {
        points[i].x = buffer[2 * i];
        points[i].y = buffer[2 * i + 1];
    }
}

int getDist(Point point1, Point point2) {
    return (point1.x - point2.x) * (point1.x - point2.x)
           + (point1.y - point2.y) * (point1.y - point2.y);
}

void assignParents(int *parents_local, Point *points_local, Point *centers, int n_local, int k) {

    for (int i = 0; i < n_local; i++) {
        int parent = -1;
        int min_dist = INT_MAX;
        for (int j = 0; j < k; j++) {
            int dist = getDist(points_local[i], centers[j]);
            if (dist < min_dist) {
                min_dist = dist;
                parent = j;
            }
        }
        parents_local[i] = parent;
    }
}

void reset(long *array, int k) {
    for (int i = 0; i < k; i++)
        array[i] = 0;
}

void findNewCenters(Point *points, Point *new_centers, int *parents, int n, int k,
                    long *x_count, long *y_count, long *count) {
    reset(x_count, k);
    reset(y_count, k);
    reset(count, k);
    for (int i = 0; i < n; i++) {
        int index = parents[i];
        x_count[index] += points[i].x;
        y_count[index] += points[i].y;
        count[index]++;
    }
    for (int i = 0; i < k; i++) {
        new_centers[i].x = x_count[i] / count[i];
        new_centers[i].y = y_count[i] / count[i];
    }
}

bool isConverge(Point *new_centers, Point *centers, int k) {
    for (int i = 0; i < k; i++) {
        if (new_centers[i].x != centers[i].x || new_centers[i].y != centers[i].y)
            return false;
    }
    return true;
}

void copy(Point *new_centers, Point *centers, int k) {
    for (int i = 0; i < k; i++) {
        centers[i].x = new_centers[i].x;
        centers[i].y = new_centers[i].y;
    }
}


int main(int argc, char *argv[]) {
    if (argc < 3) {
        cout << "Error. Usage : ./kmeans ${number of points} ${paramter k}" << endl;
        exit(1);
    }

    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    int n = atoi(argv[1]);
    int k = atoi(argv[2]);
    //cout << "Read Parameters Succesfully" << endl;
    cout << "Process " << rank << " start!" << endl;
    double start = MPI_Wtime();
    Point *points;
    Point *centers;
    Point *new_centers;
    int *parents;
    int *buffer;
    long *x_count;
    long *y_count;
    long *count;

    //step1 : initialization
    int *centerBuffer = (int *) malloc(sizeof(Point) * k);
    int n_local = n / size;
    int *buffer_local = (int *) malloc(sizeof(Point) * n_local);
    int *parents_local = (int *) malloc(sizeof(int) * n_local);
    Point *points_local = (Point *) malloc(sizeof(Point) * n_local);
    if (rank == 0) {
        //create random points
        buffer = (int *) malloc(sizeof(Point) * n);
        points = createRandomPoints(n);
        centers = createRandomCenters(points, k, n);
        new_centers = (Point *) malloc(sizeof(Point) * k);
        toBuffer(points, buffer, n);
        toBuffer(centers, centerBuffer, k);
        parents = (int *) malloc(sizeof(int) * n);
        x_count = (long *) malloc(sizeof(long) * k);
        y_count = (long *) malloc(sizeof(long) * k);
        count = (long *) malloc(sizeof(long) * k);
        //cout << "Master node init successfully" << endl;
    }
    MPI_Bcast(centerBuffer, k * 2, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Scatter(buffer, n_local * 2, MPI_INT, buffer_local, n_local * 2, MPI_INT, 0, MPI_COMM_WORLD);
    toPoints(buffer_local, points_local, n_local);
    //cout << "Broadcast information successfully" << endl;
    if (rank != 0) {
        centers = (Point *) malloc(sizeof(Point) * k);
        toPoints(centerBuffer, centers, k);
    }

    //cout << "Start iteration step" << endl;

    //start k-means
    for (int round = 0; round < MAX_ITER_ROUNDS; round++) {
        //step 2 : assignment step (assign each point to a center)
        assignParents(parents_local, points_local, centers, n_local, k);
        //cout<<"Assign parents successfully"<<endl;
        MPI_Gather(parents_local, n_local, MPI_INT, parents, n_local, MPI_INT, 0, MPI_COMM_WORLD);
        //cout<<"Assignment step"<<endl;
        //step 3: update step
        if (rank == 0) {
            findNewCenters(points, new_centers, parents, n, k, x_count, y_count, count);
            if (isConverge(new_centers, centers, k))
                break;
            copy(new_centers, centers, k);
            toBuffer(centers, centerBuffer, k);
        }
        MPI_Bcast(centerBuffer, k * 2, MPI_INT, 0, MPI_COMM_WORLD);
        if (rank != 0)
            toPoints(centerBuffer, centers, k);
    }


    if (rank == 0) {
        free(buffer);
        free(points);
        free(new_centers);
        free(count);
        free(x_count);
        free(y_count);
    }

    free(centers);
    free(centerBuffer);
    free(points_local);
    MPI_Finalize();
    if (rank == 0)
        cout << "K-means Algorithm done. Time cost : " << (MPI_Wtime() - start) / 1000 << "s" << endl;
    return 0;
}
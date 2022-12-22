#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>


void Matrix_fill(double *A, int N, int M) 
{
    int j = 0, i = 0;
    for (int k = 0; i < M * N; ++i) {
       A[k] = k;        
    }
}

double **Matrix_alloc(int N, int M) 
{
    double ** res = malloc(sizeof(double *) * M);
    for (int i = 0; i < M; ++i) {
        res[i] = malloc(sizeof(double) * N);
    }
    return res;
}

void
transpose(int N, double *src, double dst[][N], int rank, int nproc)
{
    size_t i = 0, j = 0, k = 0;
    double row[N];
    MPI_Request request;
    while (rank + i < N) {
        for (j = 0; j < N; ++j)
        {
            k = j * N + rank + i;
            if (j * N + rank + i < N * N) {
                if (!rank) {
                    dst[i][j] = row[j];
                } else {
                    row[j] = src[k];
                }
            }
        }
        if (rank) {
            MPI_Isend(row, N, MPI_DOUBLE, 0, i, MPI_COMM_WORLD, &request);
        } 
        i += nproc;
    }
}


void Matrix_print(double **A, int N, int M)
{
    for (int i = 0; i < M; ++i) {
        for (int j = 0; j < N; ++j) {
            printf("%f ", A[i][j]);
        }
        printf("\n");
    }
    printf("---------------\n");
}

void Matrix_free(double **A, int N)
{
    for (int i = 0; i < N; ++i) {
        free(A[i]);
    }
    free(A);
}

int main(int argc, char *argv[])
{
    int err;
    if ((err = MPI_Init(&argc, &argv))) {
        MPI_Abort(MPI_COMM_WORLD, err);
    }

    int nproc, rank;

    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int dim[11] = {100, 500, 600, 700, 800, 900, 6000, 8000, 10000, 12000, 15000};

    double curr, avg;

    MPI_Request request;
    for (int i = 0; i < 11; ++i) {
        int N = dim[i], M = dim[i];
        double src[N * M];
        double dst[M][N];
        avg = 0;
        for (int k = 0; k < 10; k++) {
            if (!rank) {
                Matrix_fill(src, N, M);
            }
            
            MPI_Barrier(MPI_COMM_WORLD);
            curr = MPI_Wtime();
            MPI_Bcast(src, N * M, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            transpose(N, src, dst, rank, nproc);
            MPI_Barrier(MPI_COMM_WORLD);
            if (!rank) {
                for (int r = 1; r < nproc; r++) {
                    int proc = 0;
                    while (r + proc < N) {
                        MPI_Irecv(dst[r + proc], N, MPI_DOUBLE, r, proc, MPI_COMM_WORLD, &request);
                        proc += nproc;
                    }
                }
            }
            MPI_Barrier(MPI_COMM_WORLD);
            avg += MPI_Wtime() - curr;
        }
        avg /= 10;

        if (!rank) {
            printf("%d %d %f\n", N, nproc, avg);
        }
    }

    MPI_Finalize();
    return 0;
}

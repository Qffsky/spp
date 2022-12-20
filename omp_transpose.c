#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>

double omp_get_wtime(void);


double **transpose(double **src, double **dst, int N, int M, int nproc)
{
    int i, j;
    #pragma omp parallel for private(i, j) shared(src, dst) num_threads(nproc)
    for (i = 0; i < M; ++i) {
	    for (j = 0; j < N; ++j) {
            dst[j][i] = src[i][j];
	    }
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

void Matrix_fill(double **A, int N, int M) 
{
    for (int i = 0; i < M; ++i) {
        for (int j = 0; j < N; ++j) {
	    A[i][j] = (i) * N + j + 1;
        }
    }
}

int main(void)
{
    printf("Hello world!\n");
    int N = 3, M = 2;
    double **a, **b, avg, curr;
    int arrN[4] = {1, 2, 4, 8};
    int V[11] = {100, 500, 1000, 2000, 4000, 5000, 6000, 8000, 10000, 12000, 15000};

    for (int volume = 0; volume < 11; ++volume) {
        N = V[volume];
        M = V[volume];
     
        a = Matrix_alloc(N, M);
        b = Matrix_alloc(M, N);

        Matrix_fill(a, N, M);

        //Matrix_print(a, N, M);
        for (int proc = 0; proc < 4; ++proc) {
            avg = 0;
            for (int attempt = 0; attempt < 10; ++attempt) {
                curr = omp_get_wtime();
                transpose(a, b, N, M, arrN[proc]);
                avg = omp_get_wtime() - curr;
            }
            avg /= 4;
            printf("%d %d %lf\n", V[volume], arrN[proc], avg);
        }
        //Matrix_print(b, M, N);
        Matrix_free(a, N);
        Matrix_free(b, M);
    }
    return 0;
}

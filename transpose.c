#include <stdio.h>
#include <stdlib.h>
#include <time.h>

double omp_get_wtime(void);


double **transpose(double **src, double **dst, int N, int M)
{
    for (int i = 0; i < M; ++i) {
	    for (int j = 0; j < N; ++j) {
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
    int arrN[1, 2, 4, 8]; 
    double **a, **b;

    a = Matrix_alloc(N, M);
    b = Matrix_alloc(M, N);

    Matrix_fill(a, N, M);

    Matrix_print(a, N, M);

    transpose(a, b, N, M);
    Matrix_print(b, M, N);

    Matrix_free(a, N);
    Matrix_free(b, M);
    
    return 0;
}

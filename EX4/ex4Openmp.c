#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

#define MAX 10000

int main()
{
    double *x, *y;
    double dot = 0.0;
    int n = MAX;

    x = (double *)malloc(n * sizeof(double));
    y = (double *)malloc(n * sizeof(double));

    for (int i = 0; i < n; i++)
    {
        x[i] = 1.0;
        y[i] = 1.0;
    }

    double start_time = omp_get_wtime();

#pragma omp parallel for reduction(+ : dot)
    for (int i = 0; i < n; i++)
    {
        dot += x[i] * y[i];
    }

    double end_time = omp_get_wtime();
    double elapsed_time = end_time - start_time;

    printf("Produto escalar paralelo: %f\n", dot);
    printf("Tempo de execucao paralelo: %f segundos\n", elapsed_time);

    free(x);
    free(y);

    return 0;
}

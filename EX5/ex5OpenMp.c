#include <stdio.h>
#include <stdlib.h>
#include <omp.h> // Include OpenMP header
#define MAX 10000

int main(int argc, char *argv[])
{
    double *x, *y, *local_x, *local_y;
    int n_bar; /* = n/p */
    double dot, local_dot;
    int i, n;
    n = MAX;
    double start_time, end_time;

    x = (double *)calloc(n, sizeof(double));
    y = (double *)calloc(n, sizeof(double));

    for (i = 0; i < n; i++)
    {
        x[i] = 1.0;
        y[i] = 1.0;
    }

    // Criação dos vetores locais
    n_bar = n;
    local_x = (double *)calloc(n_bar, sizeof(double));
    local_y = (double *)calloc(n_bar, sizeof(double));

// Envio dos vetores para todos os processos
#pragma omp parallel for
    for (i = 0; i < n_bar; i++)
    {
        local_x[i] = x[i];
        local_y[i] = y[i];
    }

    // Record start time
    start_time = omp_get_wtime();

// Calcula o produto escalar local in parallel
#pragma omp parallel for reduction(+ : local_dot)
    for (i = 0; i < n_bar; i++)
    {
        local_dot += local_x[i] * local_y[i];
    }

    // grava o tempo
    end_time = omp_get_wtime();

    free(x);
    free(y);
    free(local_x);
    free(local_y);

    // Coleta resultados
    dot = local_dot;

    printf("valor final = %f \n", dot);
    printf("Tempo de execução: %f segundos\n", end_time - start_time);

    return 0;
} /* main */

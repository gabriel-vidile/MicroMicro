#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <mpi.h>
#define MAX 10000

int main(int argc, char *argv[])
{

    double *x, *y, *local_x, *local_y;
    int n_bar; /* = n/p */
    double dot, local_dot;
    int p, my_rank, i, n;
    n = MAX;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    double start_time, end_time;

    if (my_rank == 0)
    {
        x = (double *)calloc(n, sizeof(double));
        y = (double *)calloc(n, sizeof(double));

        for (i = 0; i < n; i++)
        {
            x[i] = 1.0;
            y[i] = 1.0;
        }
    }
    // Envio do valor de n para todos os processos
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // Criação dos vetores locais
    n_bar = n / p;
    local_x = (double *)calloc(n_bar, sizeof(double));
    local_y = (double *)calloc(n_bar, sizeof(double));

    // Envio dos vetores para todos os processos
    MPI_Scatter(x, n_bar, MPI_DOUBLE, local_x, n_bar, MPI_DOUBLE, 0,
                MPI_COMM_WORLD);
    MPI_Scatter(y, n_bar, MPI_DOUBLE, local_y, n_bar, MPI_DOUBLE, 0,
                MPI_COMM_WORLD);

    // Record start time
    start_time = MPI_Wtime();

    // Calcula o produto escalar local
    local_dot = 0.0;
    for (i = 0; i < n_bar; i++)
    {
        local_dot += local_x[i] * local_y[i];
    }

    // grava o tempo
    end_time = MPI_Wtime();

    free(local_x);
    free(local_y);

    // Coleta resultados
    MPI_Reduce(&local_dot, &dot, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    if (my_rank == 0)
    {
        printf("valor final = %f \n", dot);
        printf("Tempo de execução: %f segundos\n", end_time - start_time);
    }

    MPI_Finalize();

    return 0;
} /* main */

#include <stdio.h>
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

    int mpi_error;
    mpi_error = MPI_Init(&argc, &argv);
    if (mpi_error != MPI_SUCCESS)
    {
        fprintf(stderr, "Problema em inicializar o MPi\n");
        MPI_Abort(MPI_COMM_WORLD, mpi_error);
    }

    MPI_Comm_size(MPI_COMM_WORLD, &p);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    double start_time, end_time;

    if (my_rank == 0)
    {
        x = (double *)malloc(n * sizeof(double));
        y = (double *)malloc(n * sizeof(double));

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
    local_x = (double *)malloc(n_bar * sizeof(double));
    local_y = (double *)malloc(n_bar * sizeof(double));

    // Envio dos vetores para todos os processos
    MPI_Scatter(x, n_bar, MPI_DOUBLE, local_x, n_bar, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Scatter(y, n_bar, MPI_DOUBLE, local_y, n_bar, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // Record start time
    start_time = MPI_Wtime();

    // Calcula o produto escalar local
    local_dot = 0.0;
    for (i = 0; i < n_bar; i += 4)
    {
        local_dot += local_x[i] * local_y[i] +
                     local_x[i + 1] * local_y[i + 1] +
                     local_x[i + 2] * local_y[i + 2] +
                     local_x[i + 3] * local_y[i + 3];
    }

    // Record end time
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

    if (my_rank == 0)
    {
        free(x);
        free(y);
    }

    MPI_Finalize();

    return 0;
}

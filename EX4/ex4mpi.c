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
    int tag = 0;
    MPI_Status status;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    double start_time, end_time;

    if (my_rank == 0)
    {
        n = MAX;
        x = (double *)calloc(n, sizeof(double));
        y = (double *)calloc(n, sizeof(double));

        for (i = 0; i < n; i++)
        {
            x[i] = 1.0;
            y[i] = 1.0;
        }
    }

    start_time = MPI_Wtime(); // Início da medição de tempo

    // Envio do valor de n para todos os processos
    if (my_rank == 0)
    {
        for (i = 1; i < p; i++)
        {
            MPI_Send(&n, 1, MPI_INT, i, tag, MPI_COMM_WORLD);
        }
    }
    else
    {
        MPI_Recv(&n, 1, MPI_INT, 0, tag, MPI_COMM_WORLD, &status);
    }

    // printf("Meu id: %d Received: %d \n", my_rank, n);

    // Criação dos vetores locais
    n_bar = n / p;
    local_x = (double *)calloc(n_bar, sizeof(double));
    local_y = (double *)calloc(n_bar, sizeof(double));

    // Envio dos vetores para todos os processos
    if (my_rank == 0)
    {
        for (i = 1; i < p; i++)
        {
            MPI_Send(&(x[i * n_bar]), n_bar, MPI_DOUBLE, i, tag, MPI_COMM_WORLD);
            MPI_Send(&(y[i * n_bar]), n_bar, MPI_DOUBLE, i, tag, MPI_COMM_WORLD);
        }
        for (i = 0; i < n_bar; i++)
        {
            local_x[i] = x[i];
            local_y[i] = y[i];
        }
    }
    else
    {
        MPI_Recv(local_x, n_bar, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD, &status);
        MPI_Recv(local_y, n_bar, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD, &status);
    }

    // Calcula o produto escalar local
    local_dot = 0.0;
    for (i = 0; i < n_bar; i++)
    {
        local_dot += local_x[i] * local_y[i];
    }

    free(local_x);
    free(local_y);

    // Coleta resultados

    if (my_rank != 0)
    {
        MPI_Send(&local_dot, 1, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD);
    }
    else
    {
        for (i = 1; i < p; i++)
        {
            MPI_Recv(&dot, 1, MPI_DOUBLE, i, tag, MPI_COMM_WORLD, &status);
        }

        local_dot += dot;
    }

    end_time = MPI_Wtime(); // Fim da medição de tempo

    if (my_rank == 0)
    {
        printf("valor final = %f \n", local_dot);
        printf("Tempo de execução: %f segundos\n", end_time - start_time);
    }

    MPI_Finalize();
    return 0;
}

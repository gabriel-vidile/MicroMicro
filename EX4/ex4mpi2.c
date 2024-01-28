#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#define MAX 10000

void initialize_vectors(double *x, double *y, int n)
{
    for (int i = 0; i < n; i++)
    {
        x[i] = 1.0;
        y[i] = 1.0;
    }
}

void scatter_vectors(double *x, double *y, double *local_x, double *local_y, int n, int n_bar, int my_rank, int p)
{
    MPI_Scatter(x, n_bar, MPI_DOUBLE, local_x, n_bar, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Scatter(y, n_bar, MPI_DOUBLE, local_y, n_bar, MPI_DOUBLE, 0, MPI_COMM_WORLD);
}

double calculate_local_dot(double *local_x, double *local_y, int n_bar)
{
    double local_dot = 0.0;
    for (int i = 0; i < n_bar; i++)
    {
        local_dot += local_x[i] * local_y[i];
    }
    return local_dot;
}

int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);

    int p, my_rank, n, n_bar;
    double *x, *y, *local_x, *local_y, dot, local_dot;
    int tag = 0;
    MPI_Status status;

    MPI_Comm_size(MPI_COMM_WORLD, &p);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    double start_time, end_time;

    if (my_rank == 0)
    {
        n = MAX;
        x = (double *)calloc(n, sizeof(double));
        y = (double *)calloc(n, sizeof(double));
        initialize_vectors(x, y, n);
    }

    start_time = MPI_Wtime(); // Início da medição de tempo

    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

    n_bar = n / p;
    local_x = (double *)malloc(n_bar * sizeof(double));
    local_y = (double *)malloc(n_bar * sizeof(double));

    scatter_vectors(x, y, local_x, local_y, n, n_bar, my_rank, p);

    local_dot = calculate_local_dot(local_x, local_y, n_bar);

    free(local_x);
    free(local_y);

    if (my_rank != 0)
    {
        MPI_Send(&local_dot, 1, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD);
    }
    else
    {
        for (int i = 1; i < p; i++)
        {
            MPI_Recv(&dot, 1, MPI_DOUBLE, i, tag, MPI_COMM_WORLD, &status);
            local_dot += dot;
        }
    }

    end_time = MPI_Wtime(); // Fim da medição de tempo

    if (my_rank == 0)
    {
        printf("valor final = %f \n", local_dot);
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

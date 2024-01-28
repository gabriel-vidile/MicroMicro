#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define NI 200
#define NJ 200
#define NSTEPS 500

int main(int argc, char *argv[])
{
    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int i, j, n, im, ip, jm, jp, ni, nj, nsum, isum;
    int **old, **new;
    float x;

    // Dividir o tabuleiro igualmente entre os processos
    int rows_per_process = NI / size;
    ni = rows_per_process + 2; // Adicionar células fantasmas
    nj = NJ + 2;

    // Alocação dos arrays
    old = malloc(ni * sizeof(int *));
    new = malloc(ni * sizeof(int *));
    for (i = 0; i < ni; i++)
    {
        old[i] = malloc(nj * sizeof(int));
        new[i] = malloc(nj * sizeof(int));
    }

    // Inicialização dos elementos
    srand(time(NULL) + rank);
    for (i = 1; i <= rows_per_process; i++)
    {
        for (j = 1; j <= NJ; j++)
        {
            x = rand() / ((float)RAND_MAX + 1);
            old[i][j] = x < 0.5 ? 0 : 1;
        }
    }

    MPI_Barrier(MPI_COMM_WORLD); // Sincronizar todos os processos
    double start_time = MPI_Wtime();

    // Passos de tempo
    for (n = 0; n < NSTEPS; n++)
    {
        // Atualizar células fantasmas aqui (usar MPI_Sendrecv)

        // Atualizar o estado do tabuleiro
        for (i = 1; i <= rows_per_process; i++)
        {
            for (j = 1; j <= NJ; j++)
            {
                im = i - 1;
                ip = i + 1;
                jm = j - 1;
                jp = j + 1;

                nsum = old[im][jp] + old[i][jp] + old[ip][jp] + old[im][j] + old[ip][j] + old[im][jm] + old[i][jm] + old[ip][jm];

                switch (nsum)
                {
                case 3:
                    new[i][j] = 1;
                    break;
                case 2:
                    new[i][j] = old[i][j];
                    break;
                default:
                    new[i][j] = 0;
                }
            }
        }

        // Trocar os ponteiros old e new
        int **temp = old;
        old = new;
        new = temp;

        MPI_Barrier(MPI_COMM_WORLD); // Sincronizar todos os processos após cada iteração
    }

    double end_time = MPI_Wtime();

    // Calcular o número de células vivas no segmento de cada processo
    isum = 0;
    for (i = 1; i <= rows_per_process; i++)
    {
        for (j = 1; j <= NJ; j++)
        {
            isum += new[i][j];
        }
    }

    int total_sum;
    MPI_Reduce(&isum, &total_sum, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

    if (rank == 0)
    {
        printf("\nNúmero total de células vivas = %d\n", total_sum);
        printf("Tempo de execução: %f segundos\n", end_time - start_time);
    }

    // Limpeza
    for (i = 0; i < ni; i++)
    {
        free(old[i]);
        free(new[i]);
    }
    free(old);
    free(new);

    MPI_Finalize();
    return 0;
}
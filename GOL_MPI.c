#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <mpi.h>
#define NI 200        
#define NJ 200
#define NSTEPS 500 

// mpicc -o GOL_MPI GOL_MPI.c 
void celulas(int **grid, int ni, int nj) {
    for (int i = 0; i < ni; i++) {
        for (int j = 0; j < nj; j++) {
            float x = rand() / ((float)RAND_MAX + 1);
            grid[i][j] = x < 0.5 ? 0 : 1;
        }
    }
}

void atualiza_celulas(int **old, int **new, int ni, int nj) {
    for (int i = 1; i <= ni; i++) {
        for (int j = 1; j <= nj; j++) {
            int im = i - 1;
            int ip = i + 1;
            int jm = j - 1;
            int jp = j + 1;

            int nsum = old[im][jp] + old[i][jp] + old[ip][jp] +
                       old[im][j] + old[ip][j] +
                       old[im][jm] + old[i][jm] + old[ip][jm];

            switch (nsum) {
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
}

int main(int argc, char *argv[]) {
    int i, j, n, im, ip, jm, jp, ni, nj, nsum, isum;
    int **old, **new;

    int rank, num_procs;

    struct timeval start, end;
    gettimeofday(&start, NULL);

    // Inicializar MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

    ni = NI / num_procs;
    nj = NJ;

    old = malloc((ni + 2) * sizeof(int*));
    new = malloc((ni + 2) * sizeof(int*));

    for (i = 0; i < ni + 2; i++) {
        old[i] = malloc(nj * sizeof(int));
        new[i] = malloc(nj * sizeof(int));
    }

    celulas(old, ni + 2, nj);

    for (n = 0; n < NSTEPS; n++) {
        if (rank > 0) {
            MPI_Send(old[1], nj, MPI_INT, rank - 1, 0, MPI_COMM_WORLD);
            MPI_Recv(old[0], nj, MPI_INT, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
        if (rank < num_procs - 1) {
            MPI_Send(old[ni], nj, MPI_INT, rank + 1, 0, MPI_COMM_WORLD);
            MPI_Recv(old[ni + 1], nj, MPI_INT, rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }

        atualiza_celulas(old, new, ni, nj);

       //copia o novo estado para o antigo
        for (i = 1; i <= ni; i++) {
            for (j = 0; j < nj; j++) {
                old[i][j] = new[i][j];
            }
        }
    }

    // soma as celulas
    isum = 0;
    for (i = 1; i <= ni; i++) {
        for (j = 0; j < nj; j++) {
            isum = isum + new[i][j];
        }
    }
    int total_isum = 0;
    MPI_Reduce(&isum, &total_isum, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

    // Finalizar MPI
    MPI_Finalize();

    return 0;
}
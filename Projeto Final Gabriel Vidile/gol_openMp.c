#include <stdio.h>
#include <stdlib.h>
#include <omp.h> // Incluir a biblioteca OpenMP

#define NI 200 /* Tamanhos dos arrays */
#define NJ 200
#define NSTEPS 500 /* Número de passos de tempo */

int main(int argc, char *argv[])
{
    int i, j, n, im, ip, jm, jp, ni, nj, nsum, isum;
    int **old, **new;
    float x;
    double start_time, end_time; // Variáveis para medir o tempo

    /* Alocar arrays */
    ni = NI + 2; /* Adicionar 2 para células fantasmas à esquerda e à direita */
    nj = NJ + 2;
    old = malloc(ni * sizeof(int *));
    new = malloc(ni * sizeof(int *));

    for (i = 0; i < ni; i++)
    {
        old[i] = malloc(nj * sizeof(int));
        new[i] = malloc(nj * sizeof(int));
    }

/* Inicializar elementos de old com 0 ou 1 */
#pragma omp parallel for private(j, x) // Paralelizar este loop
    for (i = 1; i <= NI; i++)
    {
        for (j = 1; j <= NJ; j++)
        {
            x = rand() / ((float)RAND_MAX + 1);
            if (x < 0.5)
            {
                old[i][j] = 0;
            }
            else
            {
                old[i][j] = 1;
            }
        }
    }

    start_time = omp_get_wtime(); // Iniciar medição de tempo

    /* Passos de tempo */
    for (n = 0; n < NSTEPS; n++)
    {

        /* Condições de contorno nos cantos */
        old[0][0] = old[NI][NJ];
        old[0][NJ + 1] = old[NI][1];
        old[NI + 1][NJ + 1] = old[1][1];
        old[NI + 1][0] = old[1][NJ];

        /* Condições de contorno esquerda-direita */
        for (i = 1; i <= NI; i++)
        {
            old[i][0] = old[i][NJ];
            old[i][NJ + 1] = old[i][1];
        }

        /* Condições de contorno topo-fundo */
        for (j = 1; j <= NJ; j++)
        {
            old[0][j] = old[NI][j];
            old[NI + 1][j] = old[1][j];
        }

/* Atualizar estado com base em vizinhos */
#pragma omp parallel for private(j, im, ip, jm, jp, nsum) // Paralelizar este loop
        for (i = 1; i <= NI; i++)
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

/* Copiar novo estado para o estado antigo */
#pragma omp parallel for private(j) // Paralelizar este loop
        for (i = 1; i <= NI; i++)
        {
            for (j = 1; j <= NJ; j++)
            {
                old[i][j] = new[i][j];
            }
        }
    }

    end_time = omp_get_wtime(); // Fim da medição de tempo

    /* Iterações concluídas; somar o número de células vivas */
    isum = 0;
#pragma omp parallel for reduction(+ : isum) private(j) // Paralelizar este loop
    for (i = 1; i <= NI; i++)
    {
        for (j = 1; j <= NJ; j++)
        {
            isum = isum + new[i][j];
        }
    }

    printf("\nNúmero de células vivas = %d\n", isum);
    printf("Tempo de execução: %f segundos\n", end_time - start_time);

    /* Liberar memória */
    for (i = 0; i < ni; i++)
    {
        free(old[i]);
        free(new[i]);
    }
    free(old);
    free(new);

    return 0;
}

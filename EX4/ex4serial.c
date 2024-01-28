#include <stdio.h>
#include <stdlib.h>
#include <time.h>

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

    clock_t start_time, end_time;

    start_time = clock(); // Início da medição de tempo

    // Calcula o produto escalar serial
    for (int i = 0; i < n; i++)
    {
        dot += x[i] * y[i];
    }

    end_time = clock(); // Fim da medição de tempo

    printf("Produto escalar serial: %f\n", dot);
    printf("Tempo de execucao serial: %f segundos\n", (double)(end_time - start_time) / CLOCKS_PER_SEC);

    free(x);
    free(y);

    return 0;
}

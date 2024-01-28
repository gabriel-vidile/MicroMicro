#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define MAX 10000

int main()
{
    double *x, *y;
    int n, i;
    double dot = 0.0;

    n = MAX;

    x = (double *)calloc(n, sizeof(double));
    y = (double *)calloc(n, sizeof(double));

    for (i = 0; i < n; i++)
    {
        x[i] = 1.0;
        y[i] = 1.0;
    }

    // Record start time
    clock_t start_time = clock();

    // Calcula o produto escalar local
    for (i = 0; i < n; i++)
    {
        dot += x[i] * y[i];
    }

    // Record end time
    clock_t end_time = clock();

    free(x);
    free(y);

    // Calculate elapsed time
    double elapsed_time = (double)(end_time - start_time) / CLOCKS_PER_SEC;

    printf("valor final = %f \n", dot);
    printf("Tempo de execução: %f segundos\n", elapsed_time);

    return 0;
}

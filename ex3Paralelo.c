#include <stdlib.h>
#include <stdio.h>
#include <omp.h>

int main(int argc, char *argv[])
{
    double sum = 0.0;
    double a[256], b[256];
    int n = 256;
    int i;

// Inicialize os arrays a e b paralelamente
#pragma omp parallel for
    for (i = 0; i < n; i++)
    {
        a[i] = i * 0.5;
        b[i] = i * 2.0;
    }

// Realize a soma paralelamente
#pragma omp parallel for reduction(+ : sum)
    for (i = 1; i < n; i++)
    {
        sum = sum + a[i] * b[i];
    }

    printf("sum = %f\n", sum);

    return 0;
}

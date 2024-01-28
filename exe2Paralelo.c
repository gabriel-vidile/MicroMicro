#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

int num_segments = 1000; // Variável global para o número de segmentos

int main(int argc, char *argv[])
{
    int i;
    double mypi, h, sum, x;

    h = 1.0 / (double)num_segments;
    sum = 0.0;

#pragma omp parallel for private(x) reduction(+ : sum)
    for (i = 1; i <= num_segments; i++)
    {
        x = h * ((double)i - 0.5);
        sum += 4.0 / (1.0 + x * x);
    }

    mypi = h * sum;
    printf("pi is approximately %.16f\n", mypi);

    return 0;
}

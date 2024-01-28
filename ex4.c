#include <stdlib.h>
#include <stdio.h>
#include <omp.h>
#define MAX 10000

// gcc paralelo_d2_pe.c -o paralelo_d2_pe -fopenmp

int main(int argc, char *argv[])
{

    double prod_esc;
    double a[MAX], b[MAX];

    // preenche vetores
    for (int i = 0; i < MAX; i++)
    {
        a[i] = i * 0.5;
        b[i] = i * 2.0;
    }
    prod_esc = 0;

#pragma omp parallel sections reduction(+ : prod_esc)
    {
#pragma omp section
        {
            double parcial_prod_esc = 0;

            for (int i = 0; i < (MAX / 2); i++)
            {
                parcial_prod_esc = parcial_prod_esc + a[i] * b[i];
            }

#pragma omp critical
            {
                prod_esc += parcial_prod_esc;
            }
        }

#pragma omp section
        {
            int parcial_prod_esc = 0;

            for (int i = (MAX / 2); i < MAX; i++)
            {
                parcial_prod_esc = parcial_prod_esc + a[i] * b[i];
            }

#pragma omp critical
            {
                prod_esc += parcial_prod_esc;
            }
        }
    }
    printf("produto escalar = %f\n", prod_esc);
}
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <omp.h>

#define PI2 6.28318530718
#define R_ERROR 0.01

//gcc -fopenmp FFTW_OMP.c -o FFTW_OMP -lm

int DFT(int idft, double *xr, double *xi, double *Xr_o, double *Xi_o, int N);
int preenche_entrada(double *xr, double *xi, int N);
int entrada_zero(double *Xr_o, double *Xi_o, int N);
int verifica_resultados(double *xr, double *xi, double *xr_check, double *xi_check, double *Xr_o, double *Xi_r, int N);
int imprime_resultados(double *xr, double *xi, int N);

int main(int argc, char *argv[]){

     // Tamanho do vetor de entrada
    int N = 8000;
    printf("Cálculo da DFTW com N = %d \n", N);

    // Alocar espaço para o vetor de entrada
    double *xr = (double *)malloc(N * sizeof(double));
    double *xi = (double *)malloc(N * sizeof(double));
    preenche_entrada(xr, xi, N);

    double *xr_check = (double *)malloc(N * sizeof(double));
    double *xi_check = (double *)malloc(N * sizeof(double));
    entrada_zero(xr_check, xi_check, N);

    // Alocar espaço para o vetor de saída
    double *Xr_o = (double *)malloc(N * sizeof(double));
    double *Xi_o = (double *)malloc(N * sizeof(double));
    entrada_zero(Xr_o, Xi_o, N);

    //inicia o contador
    double start_time = omp_get_wtime();

    int idft = 1;
    DFT(idft, xr, xi, Xr_o, Xi_o, N);
    idft = -1;
    DFT(idft, Xr_o, Xi_o, xr_check, xi_check, N);

    //para o contador
    double run_time = omp_get_wtime() - start_time;
    printf("Cálculo DFTW em %f segundos\n", run_time);

    verifica_resultados(xr, xi, xr_check, xi_check, Xr_o, Xi_o, N);

    #ifdef DEBUG
        imprime_resultados(Xr_o, Xi_o, N);
    #endif

    // libera memória
    free(xr);
    free(xi);
    free(Xi_o);
    free(Xr_o);
    free(xr_check);
    free(xi_check);

    return 1;
}

int DFT(int idft, double *xr, double *xi, double *Xr_o, double *Xi_o, int N){

    #pragma omp parallel for
    for (int k = 0; k < N; k++){
        for (int n = 0; n < N; n++){
            Xr_o[k] += xr[n] * cos(n * k * PI2 / N) + idft * xi[n] * sin(n * k * PI2 / N);
            Xi_o[k] += -idft * xr[n] * sin(n * k * PI2 / N) + xi[n] * cos(n * k * PI2 / N);
        }
    }

    if (idft == -1){
        #pragma omp parallel for
        for (int n = 0; n < N; n++){
            Xr_o[n] /= N;
            Xi_o[n] /= N;
        }
    }
    return 1;
}

int preenche_entrada(double *xr, double *xi, int N){
    srand(time(0));
    for (int n = 0; n < 100000; n++){
        rand();
    }
    for (int n = 0; n < N; n++){
        xr[n] = 1.0;
        xi[n] = 0.0;
    }
    return 1;
}

int entrada_zero(double *Xr_o, double *Xi_o, int N){
    #pragma omp parallel for
    for (int n = 0; n < N; n++){
        Xr_o[n] = 0.0;
        Xi_o[n] = 0.0;
    }
    return 1;
}

int verifica_resultados(double *xr, double *xi, double *xr_check, double *xi_check, double *Xr_o, double *Xi_r, int N){
    for (int n = 0; n < N; n++){
        if (fabs(xr[n] - xr_check[n]) > R_ERROR){
            printf("ERROR - x[%d] = %f, inv(X)[%d]=%f \n", n, xr[n], n, xr_check[n]);
        }
        if (fabs(xi[n] - xi_check[n]) > R_ERROR){
            printf("ERROR - x[%d] = %f, inv(X)[%d]=%f \n", n, xi[n], n, xi_check[n]);
        }
    }
    printf("Xre[0] = %f \n", Xr_o[0]);
    return 1;
}

int imprime_resultados(double *xr, double *xi, int N){
    for (int n = 0; n < N; n++){
        printf("Xre[%d] = %f, Xim[%d] = %f \n", n, xr[n], n, xi[n]);
    }
    return 1;
}

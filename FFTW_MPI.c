//bibliotecas
#include "stdio.h"    
#include "stdlib.h"   
#include "time.h"     
#include "math.h"     
#include <mpi.h>
#define PI2 6.28318530718
#define R_ERROR 0.01

//mpicc -o FFTW_MPI FFTW_MPI.c -lm

int DFT(int idft, double *xr, double *xi, double *Xr_o, double *Xi_o, int N, int myrank, int num_procs);
int preenche_entrada(double *xr, double *xi, int N);
int saida_zero(double *Xr_o, double *Xi_o, int N);
int verifica_resultados(double *xr, double *xi, double *xr_check, double *xi_check, double *Xr_o, double *Xi_r, int N);
int imprime_resultados(double *xr, double *xi, int N);

int main(int argc, char *argv[]){
    // Inicializar o MPI
    MPI_Init(&argc, &argv);

    // Obter o rank e o número de processos
    int myrank, num_procs;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

    // Tamanho do vetor de entrada
    int N = 8000; 
    printf("Cálculo da DFTW com N = %d \n", N);

    // Alocar espaço para o vetor de entrada
    double *xr = (double *)malloc(N * sizeof(double));
    double *xi = (double *)malloc(N * sizeof(double));
    preenche_entrada(xr, xi, N);

    double *xr_check = (double *)malloc(N * sizeof(double));
    double *xi_check = (double *)malloc(N * sizeof(double));
    saida_zero(xr_check, xi_check, N);

    // Alocar espaço para o vetor de saída
    double *Xr_o = (double *)malloc(N * sizeof(double));
    double *Xi_o = (double *)malloc(N * sizeof(double));
    saida_zero(Xr_o, Xi_o, N);

    // Inicia contador
    double start_time = MPI_Wtime();

    // DFT
    int idft = 1;
    DFT(idft, xr, xi, Xr_o, Xi_o, N, myrank, num_procs);
    // IDFT
    idft = -1;
    DFT(idft, Xr_o, Xi_o, xr_check, xi_check, N, myrank, num_procs);

    // Pausa o contador
    double run_time = MPI_Wtime() - start_time;
    printf("Cálculo DFTW em %f segundos\n", run_time);

    verifica_resultados(xr, xi, xr_check, xi_check, Xr_o, Xi_o, N);

#ifdef DEBUG
    imprime_resultados(Xr_o, Xi_o, N);
#endif

    // libera memoria
    free(xr);
    free(xi);
    free(Xi_o);
    free(Xr_o);
    free(xr_check);
    free(xi_check);

    // Finaliza o MPI
    MPI_Finalize();

    return 1;
}

int DFT(int idft, double *xr, double *xi, double *Xr_o, double *Xi_o, int N, int myrank, int num_procs)
{
    int start = (N / num_procs) * myrank;
    int end = (N / num_procs) * (myrank + 1);

    for (int k = start; k < end; k++){
        for (int n = 0; n < N; n++){
            // real
            Xr_o[k] += xr[n] * cos(n * k * PI2 / N) + idft * xi[n] * sin(n * k * PI2 / N);
            // imaginario
            Xi_o[k] += -idft * xr[n] * sin(n * k * PI2 / N) + xi[n] * cos(n * k * PI2 / N);
        }
    }

    if (idft == -1){
        for (int n = start; n < end; n++){
            Xr_o[n] /= N;
            Xi_o[n] /= N;
        }
    }

    // Gather results from all processes
    MPI_Allgather(MPI_IN_PLACE, N / num_procs, MPI_DOUBLE, Xr_o, N / num_procs, MPI_DOUBLE, MPI_COMM_WORLD);
    MPI_Allgather(MPI_IN_PLACE, N / num_procs, MPI_DOUBLE, Xi_o, N / num_procs, MPI_DOUBLE, MPI_COMM_WORLD);

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

int saida_zero(double *Xr_o, double *Xi_o, int N){
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
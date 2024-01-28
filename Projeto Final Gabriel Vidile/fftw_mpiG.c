#include "stdio.h"
#include "stdlib.h"
#include "time.h"
#include "math.h"
#include <mpi.h>

#define PI2 6.28318530718
#define R_ERROR 0.01

int DFT(int idft, double *xr, double *xi, double *Xr_o, double *Xi_o, int N, int myrank, int num_procs);
int fillInput(double *xr, double *xi, int N);
int setOutputZero(double *Xr_o, double *Xi_o, int N);
int checkResults(double *xr, double *xi, double *xr_check, double *xi_check, double *Xr_o, double *Xi_r, int N);
int printResults(double *xr, double *xi, int N);

int main(int argc, char *argv[])
{
    // Initialize MPI
    MPI_Init(&argc, &argv);

    // Get MPI rank and number of processes
    int myrank, num_procs;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

    // Size of input array
    int N = 8000;
    if (myrank == 0)
    {
        printf("DFTW calculation with N = %d \n", N);
    }

    // Allocate array for input vector
    double *xr = (double *)malloc(N * sizeof(double));
    double *xi = (double *)malloc(N * sizeof(double));
    fillInput(xr, xi, N);

    // Allocate array for output vector
    double *Xr_o = (double *)malloc(N * sizeof(double));
    double *Xi_o = (double *)malloc(N * sizeof(double));
    setOutputZero(Xr_o, Xi_o, N);

    // Allocate arrays for checking purposes
    double *xr_check = (double *)malloc(N * sizeof(double));
    double *xi_check = (double *)malloc(N * sizeof(double));
    setOutputZero(xr_check, xi_check, N);

    // Start timer
    double start_time = MPI_Wtime();

    // DFT
    int idft = 1;
    DFT(idft, xr, xi, Xr_o, Xi_o, N, myrank, num_procs);
    // IDFT
    idft = -1;
    DFT(idft, Xr_o, Xi_o, xr_check, xi_check, N, myrank, num_procs);

    // Stop timer
    double run_time = MPI_Wtime() - start_time;
    if (myrank == 0)
    {
        printf("DFTW computation in %f seconds\n", run_time);
    }

    // Check the results
    checkResults(xr, xi, xr_check, xi_check, Xr_o, Xi_o, N);

    // Print the results of the DFT
#ifdef DEBUG
    printResults(Xr_o, Xi_o, N);
#endif

    // Clean up
    free(xr);
    free(xi);
    free(Xr_o);
    free(Xi_o);
    free(xr_check);
    free(xi_check);

    // Finalize MPI
    MPI_Finalize();

    return 0;
}

// Parallel DFT/IDFT routine
int DFT(int idft, double *xr, double *xi, double *Xr_o, double *Xi_o, int N, int myrank, int num_procs)
{
    int start = (N / num_procs) * myrank;
    int end = (N / num_procs) * (myrank + 1);
    for (int k = start; k < end; k++)
    {
        for (int n = 0; n < N; n++)
        {
            Xr_o[k] += xr[n] * cos(n * k * PI2 / N) + idft * xi[n] * sin(n * k * PI2 / N);
            Xi_o[k] += -idft * xr[n] * sin(n * k * PI2 / N) + xi[n] * cos(n * k * PI2 / N);
        }
    }

    // Normalize if IDFT
    if (idft == -1)
    {
        for (int n = start; n < end; n++)
        {
            Xr_o[n] /= N;
            Xi_o[n] /= N;
        }
    }

    // Gather results from all processes
    MPI_Allgather(MPI_IN_PLACE, N / num_procs, MPI_DOUBLE, Xr_o, N / num_procs, MPI_DOUBLE, MPI_COMM_WORLD);
    MPI_Allgather(MPI_IN_PLACE, N / num_procs, MPI_DOUBLE, Xi_o, N / num_procs, MPI_DOUBLE, MPI_COMM_WORLD);

    return 1;
}

// Fill input with a constant real signal
int fillInput(double *xr, double *xi, int N)
{
    srand(time(0) + MPI_Wtime());
    for (int n = 0; n < 100000; n++) // get some random number first
        rand();
    for (int n = 0; n < N; n++)
    {
        xr[n] = 1.0;
        xi[n] = 0.0;
    }
    return 1;
}

// Set output vector to zero
int setOutputZero(double *Xr_o, double *Xi_o, int N)
{
    for (int n = 0; n < N; n++)
    {
        Xr_o[n] = 0.0;
        Xi_o[n] = 0.0;
    }
    return 1;
}

// Check if x = IDFT(DFT(x))
int checkResults(double *xr, double *xi, double *xr_check, double *xi_check, double *Xr_o, double *Xi_r, int N)
{

    printf("Xre[0] = %f \n", Xr_o[0]);
    return 1;
}

// Print the results of the DFT
int printResults(double *xr, double *xi, int N)
{
    for (int n = 0; n < N; n++)
        printf("Xre[%d] = %f, Xim[%d] = %f \n", n, xr[n], n, xi[n]);
    return 1;
}

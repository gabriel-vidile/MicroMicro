#include <stdio.h>
#include <stdlib.h>
#include <omp.h>      
#define NI 200        
#define NJ 200
#define NSTEPS 500    

//gcc -fopenmp GOL_OMP.c -o GOL_OMP

int main(int argc, char *argv[]){
  int i, j, n, im, ip, jm, jp, ni, nj, nsum, isum;
  int **old, **new;  
  float x;

  // alocação de vetores

  ni = NI + 2; 
  nj = NJ + 2; 
  old = malloc(ni*sizeof(int*));
  new = malloc(ni*sizeof(int*));

  for(i=0; i<ni; i++){
    old[i] = malloc(nj*sizeof(int));
    new[i] = malloc(nj*sizeof(int));
  }


  for(i=1; i<=NI; i++){
    for(j=1; j<=NJ; j++){
      x = rand()/((float)RAND_MAX + 1);
      if(x<0.5){
        old[i][j] = 0;
      } else {
        old[i][j] = 1;
      }
    }
  }

  for(n=0; n<NSTEPS; n++){

    // condições de borda
    old[0][0] = old[NI][NJ];
    old[0][NJ+1] = old[NI][1];
    old[NI+1][NJ+1] = old[1][1];
    old[NI+1][0] = old[1][NJ];

    // bordas direita/esquerda
#pragma omp parallel for private(i)
    for(i=1; i<=NI; i++){
      old[i][0] = old[i][NJ];
      old[i][NJ+1] = old[i][1];
    }

    // bordas cima/baixo
#pragma omp parallel for private(j)
    for(j=1; j<=NJ; j++){
      old[0][j] = old[NI][j];
      old[NI+1][j] = old[1][j];
    }

#pragma omp parallel for private(i, j, im, ip, jm, jp, nsum)
    for(i=1; i<=NI; i++){
      for(j=1; j<=NJ; j++){
        im = i-1;
        ip = i+1;
        jm = j-1;
        jp = j+1;

        nsum =  old[im][jp] + old[i][jp] + old[ip][jp]
              + old[im][j ]              + old[ip][j ] 
              + old[im][jm] + old[i][jm] + old[ip][jm];

        switch(nsum){

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

//copia o novo estado para o antigo
#pragma omp parallel for private(i, j)
    for(i=1; i<=NI; i++){
      for(j=1; j<=NJ; j++){
        old[i][j] = new[i][j];
      }
    }
  }

  //soma as celulas
  isum = 0;
  for(i=1; i<=NI; i++){
    for(j=1; j<=NJ; j++){
      isum = isum + new[i][j];
    }
  }
  printf("\nnúmero de células = %d\n", isum);

  return 0;
}

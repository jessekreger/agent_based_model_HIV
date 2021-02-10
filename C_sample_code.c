//
// main.c
//
// Agent-based model for HIV infection, including cell-to-cell transmission and recombination. Fitness landscapes can be varied by changing probtoinfect, wtsuccess, imsuccess, and dmsuccess. IMPORTANT NOTE: In order to work on a wide variety of systems, this program uses rand() and drand48() as placeholders for the rngs (without seeding). The user should update this with the best scientific-computing appropriate rng available on their system. The Mersenne Twister was used for our simulations.
//

// Effect of synaptic cell-to-cell transmission and recombination on the evolution of double mutants in HIV
// Jesse Kreger, Natalia L. Komarova, and Dominik Wodarz
// Journal of the Royal Society Interface

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

    double lambda = 1;
    double d = 0.01;
    double beta = 0;
    double a = 0.02;
    double S = 3;
    double mutationrate = 0.00003;
    double recombinationrate = 0.2;
    double c = 0;
    double alpha = 0;
    double probtoinfect = 0.1;
    double wtsuccess = 1;
    double imsuccess = 1;
    double dmsuccess = 1;
    const int N = 100;
    int T = 10000;
    int YW = 7824;
    int YU = 2000;
    int matcher;
    double sa = 1;
    double sb = 1;
    double sw = 1;
    double sr = 1;
    int numberofloops = 10;
    int radius = 1;
    int loop = 0;
    int t;
    int R;
    double gridW[N][N], gridA[N][N], gridB[N][N], gridR[N][N];
    int i, j;
    int n, w;
    double totinf = 0;
    double ninfR = 0;
    double randomprob;
    int randomnumber1, randomnumber2;
    int G1 = 0, G2 = 0;
    double fractionA, fractionB, fractionW, fractionR;
    int randomlocation1, randomlocation2, locationi = 0, locationj = 0;
    double ratio = 0;

    double cellbirth() {
        gridW[i][j] = 1;
        gridA[i][j] = 1;
        gridB[i][j] = 1;
        gridR[i][j] = 1;
        return (gridW[i][j],gridA[i][j],gridB[i][j],gridR[i][j]);
    }

    double celldeath() {
        gridW[i][j] = 0;
        gridA[i][j] = 0;
        gridB[i][j] = 0;
        gridR[i][j] = 0;
        return (gridW[i][j],gridA[i][j],gridB[i][j],gridR[i][j]);
    }

    double pickgenomes() {
      randomprob = drand48();
      fractionA = ((gridA[i][j]-1)*sa)/((gridA[i][j]-1)*sa+(gridB[i][j]-1)*sb+(gridW[i][j]-1)*sw+(gridR[i][j]-1)*sr);
      fractionB = ((gridB[i][j]-1)*sb)/((gridA[i][j]-1)*sa+(gridB[i][j]-1)*sb+(gridW[i][j]-1)*sw+(gridR[i][j]-1)*sr);
      fractionW = ((gridW[i][j]-1)*sw)/((gridA[i][j]-1)*sa+(gridB[i][j]-1)*sb+(gridW[i][j]-1)*sw+(gridR[i][j]-1)*sr);
      fractionR = ((gridR[i][j]-1)*sr)/((gridA[i][j]-1)*sa+(gridB[i][j]-1)*sb+(gridW[i][j]-1)*sw+(gridR[i][j]-1)*sr);
      if (randomprob < fractionW) {
          G1 = 0;
      }
      else if (randomprob < (fractionW + fractionA)) {
          G1 = 1;
      }
      else if (randomprob < (fractionW + fractionA + fractionB)) {
          G1 = 2;
      }
      else if (randomprob < (fractionW + fractionA + fractionB + fractionR)) {
          G1 = 3;
      }
      randomprob = drand48();
      if (randomprob < recombinationrate) {
          randomprob = drand48();
          if (randomprob < fractionW) {
              G2 = 0;
          }
          else if (randomprob < (fractionW + fractionA)) {
              G2 = 1;
          }
          else if (randomprob < (fractionW + fractionA + fractionB)) {
              G2 = 2;
          }
          else if (randomprob < (fractionW + fractionA + fractionB + fractionR)) {
              G2 = 3;
          }
          if ((G1 == 1 && G2 == 2) || (G1 == 2 && G2 == 1)) {
              randomprob = drand48();
              if (randomprob < 0.5) {
                  G1 = 0;
              }
              else {
                  G1 = 3;
              }
          }
          else if ((G1 == 0 && G2 == 3) || (G1 == 3 && G2 == 0)) {
              randomprob = drand48();
              if (randomprob < 0.5) {
                  G1 = 1;
              }
              else {
                  G1 = 2;
              }
          }
      }
      return (G1,G2);
    }

    double mutationAandB() {
        randomprob = drand48();
        if (randomprob < mutationrate) {
            if (G1 == 0) {
                G1 = 1;
            }
            else if (G1 == 1) {
                G1 = 0;
            }
            else if (G1 == 2) {
                G1 = 3;
            }
            else if (G1 == 3) {
                G1 = 2;
            }
        }
        randomprob = drand48();
        if (randomprob < mutationrate) {
            if (G1 == 0) {
                G1 = 2;
            }
            else if (G1 == 1) {
                G1 = 3;
            }
            else if (G1 == 2) {
                G1 = 0;
            }
            else if (G1 == 3) {
                G1 = 1;
            }
        }
        return(G1);
    }

    double infectcelltocell() {
        if (G1 == 0) {
            gridW[i+locationi][j+locationj] = gridW[i+locationi][j+locationj] + 1;
        }
        else if (G1 == 1) {
            gridA[i+locationi][j+locationj] = gridA[i+locationi][j+locationj] + 1;
        }
        else if (G1 == 2) {
            gridB[i+locationi][j+locationj] = gridB[i+locationi][j+locationj] + 1;
        }
        else if (G1 == 3) {
            gridR[i+locationi][j+locationj] = gridR[i+locationi][j+locationj] + 1;
        }
        return(gridW[i+locationi][j+locationj],gridA[i+locationi][j+locationj],gridB[i+locationi][j+locationj],gridR[i+locationi][j+locationj]);
    }

    double infectfreevirus() {
        if (G1 == 0) {
            gridW[randomnumber1][randomnumber2] = gridW[randomnumber1][randomnumber2] + 1;
        }
        else if (G1 == 1) {
            gridA[randomnumber1][randomnumber2] = gridA[randomnumber1][randomnumber2] + 1;
        }
        else if (G1 == 2) {
            gridB[randomnumber1][randomnumber2] = gridB[randomnumber1][randomnumber2] + 1;
        }
        else if (G1 == 3) {
            gridR[randomnumber1][randomnumber2] = gridR[randomnumber1][randomnumber2] + 1;
        }
        return(gridW[randomnumber1][randomnumber2],gridA[randomnumber1][randomnumber2],gridB[randomnumber1][randomnumber2],gridR[randomnumber1][randomnumber2]);
    }

    double celltocelltransmission() {
        R = 1;
        for (w = 0; w < R; w++) {
            randomlocation1 = rand() % (2 * radius + 1);
            randomlocation2 = rand() % (2 * radius + 1);
            locationi = randomlocation1 - radius;
            locationj = randomlocation2 - radius;
            if (locationi == 0 && locationj == 0) {
                R = R + 1;
            }
            if ((0 > (i + locationi)) || ((i + locationi) > N-1) || (0 > (j + locationj)) || ((j + locationj) > N-1)) {
                R = R + 1;
            }
        }
        if (gridW[i+locationi][j+locationj] > 0 && gridA[i+locationi][j+locationj] > 0 && gridB[i+locationi][j+locationj] > 0 && gridR[i+locationi][j+locationj] > 0) {
            for (w = 0; w < S; w++) {
                pickgenomes();
                if (G1 == 0) {
                    randomprob = drand48();
                    if (randomprob < wtsuccess) {
                        mutationAandB();
                        infectcelltocell();
                    }
                }
                else if (G1 == 1 || G1 == 2) {
                    randomprob = drand48();
                    if (randomprob < imsuccess) {
                        mutationAandB();
                        infectcelltocell();
                    }
                }
                else if (G1 == 3) {
                    randomprob = drand48();
                    if (randomprob < dmsuccess) {
                        mutationAandB();
                        infectcelltocell();
                    }
                }
            }
        }
        return (gridW[i+locationi][j+locationj],gridA[i+locationi][j+locationj],gridB[i+locationi][j+locationj],gridR[i+locationi][j+locationj]);
    }

    double freevirustransmission() {
        randomnumber1 = rand() % N;
        randomnumber2 = rand() % N;
        if (gridW[randomnumber1][randomnumber2] > 0 && gridA[randomnumber1][randomnumber2] > 0 && gridB[randomnumber1][randomnumber2] > 0 && gridR[randomnumber1][randomnumber2] > 0) {
            pickgenomes();
            if (G1 == 0) {
                randomprob = drand48();
                if (randomprob < wtsuccess) {
                    mutationAandB();
                    infectfreevirus();
                }
            }
            else if (G1 == 1 || G1 == 2) {
                randomprob = drand48();
                if (randomprob < imsuccess) {
                    mutationAandB();
                    infectfreevirus();
                }
            }
            else if (G1 == 3) {
                randomprob = drand48();
                if (randomprob < dmsuccess) {
                    mutationAandB();
                    infectfreevirus();
                }
            }
        }
        return (gridW[randomnumber1][randomnumber2],gridA[randomnumber1][randomnumber2],gridB[randomnumber1][randomnumber2],gridR[randomnumber1][randomnumber2]);
    }

// main program
int main(int argc, const char * argv[]) {
    
    FILE *file_handle;

    for (loop = 0; loop < numberofloops; loop++) {
        
        // initial condition
        for (i = 0; i < N; i++) {
            for (j = 0; j < N; j++) {
                celldeath();
            }
        }
        matcher = 0;
        for (i = 0; i < (YW + matcher); i++) {
            randomnumber1 = rand() % N;
            randomnumber2 = rand() % N;
            if (gridW[randomnumber1][randomnumber2] == 0) {
                gridW[randomnumber1][randomnumber2] = 2;
                gridA[randomnumber1][randomnumber2] = 1;
                gridB[randomnumber1][randomnumber2] = 1;
                gridR[randomnumber1][randomnumber2] = 1;
            }
            else {
                matcher = matcher + 1;
            }
        }
        matcher = 0;
        for (i = 0; i < (YU + matcher); i++) {
            randomnumber1 = rand() % N;
            randomnumber2 = rand() % N;
            if (gridW[randomnumber1][randomnumber2] == 0) {
                gridW[randomnumber1][randomnumber2] = 1;
                gridA[randomnumber1][randomnumber2] = 1;
                gridB[randomnumber1][randomnumber2] = 1;
                gridR[randomnumber1][randomnumber2] = 1;
            }
            else {
                matcher = matcher + 1;
            }
        }
        
        // time loop
        t = 0;
        while (t < T) {
            t = t + 1;
            for (n = 0; n < (N * N); n++) {
                i = rand() % N;
                j = rand() % N;
                randomprob = drand48();
                if (gridW[i][j] == 0 && gridA[i][j] == 0 && gridB[i][j] == 0 && gridR[i][j] == 0) {
                    if (randomprob < lambda) {
                        cellbirth();
                    }
                }
                else if (gridW[i][j] == 1 && gridA[i][j] == 1 && gridB[i][j] == 1 && gridR[i][j] == 1) {
                    if (randomprob < d) {
                        celldeath();
                    }
                }
                else {
                    if (randomprob < a) {
                        celldeath();
                    }
                    else if (randomprob < a + probtoinfect) {
                        randomprob = drand48();
                        if (randomprob < (1 - beta)) {
                            celltocelltransmission();
                        }
                        else {
                            freevirustransmission();
                        }
                    }
                }
            }
        }
        
        // calculate populations and print to file
        ninfR = 0; totinf = 0; ratio = 0;
        for (i = 0; i < N; i++) {
            for (j = 0; j < N; j++) {
                if (gridR[i][j] > 1) {
                    ninfR = ninfR + 1;
                }
                if (gridW[i][j] > 1 || gridA[i][j] > 1 || gridB[i][j] > 1 || gridR[i][j] > 1) {
                    totinf = totinf + 1;
                }
            }
        }
        ratio = ninfR / totinf;
        printf("Loop number %d \n", loop);
        file_handle = fopen("filename.csv", "a");
        fprintf(file_handle, "%d, %lf, %lf, %lf \n", t, ninfR, totinf, ratio);
        fclose(file_handle);
    }
    return 0;
}




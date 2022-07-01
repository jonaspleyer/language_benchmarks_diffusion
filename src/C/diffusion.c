#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))


void reaction_terms(
    int rows,
    int cols,
    double** dX,
    double** dY,
    double** X,
    double** Y,
    double t
){
    double k1 = 10.0;
    double k2 = 0.1;
    double k3 = 4.938271604938272e-07;
    double k4 = 80.0;

    int i,j;
    for (i=0; i<rows; i++) {
        for (j=0; j<cols; j++) {
            dX[i][j] += k1 - k2*X[i][j] + k3*X[i][j]*X[i][j]*Y[i][j];
            dY[i][j] += k4              - k3*X[i][j]*X[i][j]*Y[i][j];
        }
    }
}


void diffusion_terms(
    int rows,
    int cols,
    double** dX,
    double** dY,
    double** X,
    double** Y,
    double t,
    double diffusion_constant1,
    double diffusion_constant2,
    double dx,
    double dy
){
    int i, j;
    double r1x = diffusion_constant1/dx/dx;
    double r1y = diffusion_constant1/dy/dy;
    double r2x = diffusion_constant2/dx/dx;
    double r2y = diffusion_constant2/dy/dy;
    // Solve for the rest of the diffusion process
    for (i=0; i<rows; i++) {
        for (j=0; j<cols; j++) {
            dX[i][j] += 
                 r1x * (X[MIN(i+1, rows-1)][j] + X[MAX(i-1, 0)][j] - 2.0*X[i][j])
                +r1y * (X[i][MIN(j+1, cols-1)] + X[i][MAX(0, j-1)] - 2.0*X[i][j]);
            dY[i][j] += 
                 r2x * (Y[MIN(i+1, rows-1)][j] + Y[MAX(i-1, 0)][j] - 2.0*Y[i][j])
                +r2y * (Y[i][MIN(j+1, cols-1)] + Y[i][MAX(0, j-1)] - 2.0*Y[i][j]);
        }
    }
}


void output_to_file(int rows, int cols, double** X, double** Y, char* filename1, char* filename2) {
    FILE *fp1 = fopen(filename1, "w+");
    FILE *fp2 = fopen(filename2, "w+");
    for (int i=0; i<rows; i++) {
        for (int j=0; j<cols;j++) {
            if (j==cols-1) {
                fprintf(fp1, "%f\n",X[i][j]);
                fprintf(fp2, "%f\n",Y[i][j]);
            } else {
                fprintf(fp1, "%f,", X[i][j]);
                fprintf(fp2, "%f,", Y[i][j]);
            }
        }
    }
    fclose(fp1);
    fclose(fp2);
}


int main() {
    double t = 0.0;
    double t_max = 1000.0;
    double dt = 0.01;

    double x_min = -1000;
    double x_max = 1000;
    int n_x = 100;
    double dx = (x_max-x_min)/n_x;

    double y_min = -1000;
    double y_max = 1000;
    int n_y = 100;
    double dy = (y_max-y_min)/n_y;

    double diffusion_constant1 = 100.0;
    double diffusion_constant2 = 5000.0;

    int i;

    double **X;
    double **dX;
    double **Y;
    double **dY;

    /* allocate the array */
    X = malloc(n_x * sizeof *X);
    dX = malloc(n_x * sizeof *dX);
    Y = malloc(n_y * sizeof *Y);
    dY = malloc(n_y * sizeof *dY);

    for (i=0; i<n_x; i++) {
        X[i] = malloc(n_y * sizeof *X[i]);
        dX[i] = malloc(n_y * sizeof *dX[i]);
        Y[i] = malloc(n_y * sizeof *Y[i]);
        dY[i] = malloc(n_y * sizeof *dY[i]);
    }
    
    /* Fill with initial values */
    int seed = 20220701;
    srand(seed);
    double deltaX = 10.0;
    double deltaY = 5.0;
    for (int i=0; i<n_x; i++) {
        for (int j=0; j<n_y; j++) {
            X[i][j] = 900 + deltaX * (double)rand() / (double)RAND_MAX;
            dX[i][j] = 0.0;
            Y[i][j] = 200 + deltaY * (double)rand() / (double)RAND_MAX;
            dY[i][j] = 0.0;
        }
    }
    
    /* Start main solving loop */
    while (t < t_max) {
        diffusion_terms(n_x, n_y, dX, dY, X, Y, t, diffusion_constant1, diffusion_constant2, dx, dy);
        reaction_terms(n_x, n_y, dX, dY, X, Y, t);
        for (int i=0; i<n_x; i++) {
            for (int j=0; j<n_y; j++) {
                X[i][j] += dt*dX[i][j];
                Y[i][j] += dt*dY[i][j];
                dX[i][j] = 0.0;
                dY[i][j] = 0.0;
            }
        }
        t += dt;
    }
    
    /* Save the results to a file */
    output_to_file(n_x, n_y, X, Y, "results_C_1.csv", "results_C_2.csv");

    /* deallocate the array */
    for (i=0; i<n_x; i++) {
        free(X[i]);
        free(Y[i]);
    }

    free(X);
    free(Y);
    
    return 0;
}

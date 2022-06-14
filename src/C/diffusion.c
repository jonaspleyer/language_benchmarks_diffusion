#include <stdio.h>
#include <stdlib.h>
#include <math.h>


double k1 = 0.01;
double k2 = 0.1;
double k3 = 0.3;
double k4 = 0.2;

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))


void update_diffusion(int rows, int cols, double** X, double** Y, double t, double dt, double diffusion_constant, double dx, double dy)
{
    int i, j;
    double rx = diffusion_constant/dx/dx;
    double ry = diffusion_constant/dy/dy;

    for (i=0; i<rows; i++) {
        for (j=0; j<cols; j++) {
            X[i][j] += dt*(k1 - k2*X[i][j]  + k3*X[i][j]*X[i][j]*Y[i][j])
                + dt*rx*(X[MIN(i+1, rows-1)][j] - 2*X[i][j] + X[MAX(i-1, 0)][j])
                + dt*rx*(X[i][MIN(j+1, cols-1)] - 2*X[i][j] + X[i][MAX(0, j-1)]);
            Y[i][j] += dt*(k4               - k3*X[i][j]*X[i][j]*Y[i][j])
                + dt*ry*(Y[MIN(i+1, rows-1)][j] - 2*Y[i][j] + Y[MAX(i-1, 0)][j])
                + dt*ry*(Y[i][MIN(j+1, cols-1)] - 2*Y[i][j] + Y[i][MAX(0, j-1)]);
        }
    }
}


int main() {
    double t = 0.0;
    double t_max = 100.0;
    double dt = 0.01;

    double x_min = -200;
    double x_max = 200;
    int n_x = 40;
    double dx = (x_max-x_min)/dx;

    double y_min = -200;
    double y_max = 200;
    int n_y = 40;
    double dy = (y_max-y_min)/dy;

    double diffusion_constant = 300.0;

    int i;
    int rows = 1000;
    int cols = 1000;

    double **X;
    double **Y;

    /* obtain values for rows & cols */

    /* allocate the array */
    X = malloc(rows * sizeof *X);
    Y = malloc(rows * sizeof *Y);

    for (i=0; i<rows; i++) {
        X[i] = malloc(cols * sizeof *X[i]);
        Y[i] = malloc(cols * sizeof *Y[i]);
    }


    while (t < t_max) {
        update_diffusion(rows, cols, X, Y, t, dt, diffusion_constant, dx, dy);
        t += dt;
    }

    printf("%f %f\n",
        X[0][0],
        Y[0][0]
    );

    /* deallocate the array */
    for (i=0; i<rows; i++) {
        free(X[i]);
        free(Y[i]);
    }

    free(X);
    free(Y);
    
    return 0;
}
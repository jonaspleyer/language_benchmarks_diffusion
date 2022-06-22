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
    double t_max = 500.0;
    double dt = 0.1;

    double x_min = -600;
    double x_max = 600;
    int n_x = 30;
    double dx = (x_max-x_min)/dx;

    double y_min = -600;
    double y_max = 600;
    int n_y = 30;
    double dy = (y_max-y_min)/dy;

    double diffusion_constant = 100.0;

    int i;

    double **X;
    double **Y;

    /* obtain values for rows & cols */

    /* allocate the array */
    X = malloc(n_x * sizeof *X);
    Y = malloc(n_y * sizeof *Y);

    for (i=0; i<n_x; i++) {
        X[i] = malloc(n_y * sizeof *X[i]);
        Y[i] = malloc(n_y * sizeof *Y[i]);
    }

    for (int i=0; i<n_x; i++) {
        for (int j=0; j<n_y; j++) {
            X[i][j] = 0.0001;
            Y[i][j] = 0.0001;
        }
    }


    while (t < t_max) {
        update_diffusion(n_x, n_y, X, Y, t, dt, diffusion_constant, dx, dy);
        t += dt;
    }

    printf("%f %f\n",
        X[0][0],
        Y[0][0]
    );
    printf("%f %f\n",
        X[5][5],
        Y[5][5]
    );
    printf("%f %f\n",
        X[15][15],
        Y[15][15]
    );

    /* deallocate the array */
    for (i=0; i<n_x; i++) {
        free(X[i]);
        free(Y[i]);
    }

    free(X);
    free(Y);
    
    return 0;
}

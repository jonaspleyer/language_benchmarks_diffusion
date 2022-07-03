#include <stdio.h>
#include <stdlib.h>
#include <math.h>


// Define common types for matrices
// and boundary values

typedef struct {
    int size_x;
    int size_y;
    double** values;
} Matrix2D;


typedef struct {
    int size_x;
    int size_y;
    int size_z;
    double*** values;
} Matrix3D;


typedef struct {
    char type[256];
    Matrix2D values;
} BoundaryValues;


typedef struct {
    Matrix3D a;
    Matrix3D b;
    Matrix3D c;
    Matrix3D cp;
} ThomasConstants;


// Define function prototypes. Specify them down below
void print_Matrix2D(Matrix2D M);
void print_Matrix3D(Matrix3D N);
Matrix2D get_2D_random_initial_values(double min, double max, int size_x, int size_y);
Matrix3D get_initial_values(int n_components, double min[], double max[], int size_x, int size_y);
void add_reaction_terms(Matrix3D &dN, const Matrix3D &N);
ThomasConstants define_thomas_constants(const int N, const double* D, const double &dx, const double &dy, const double &dt);


int main(int argc, char* argv[]) {
    // Initialize the random seed
    srand(2022);
    
    // Get sizes from arguments of function
    // int size_x = atoi(argv[1]);
    // int size_y = atoi(argv[2]);
    int n_components = 2;
    int size_x = 10;
    int size_y = 10;

    // Define values for diffusion constants
    double D[n_components] = {100.0, 5000.0};

    // Define values for spatial separation
    double dx = 20.0;
    double dy = 20.0;

    // Define value for time interval to solve over
    double t0 = 0.0;
    double tmax = 100.0;
    double dt = 0.1;
    
    // Define initial values
    double mins[n_components] = {200, 900};
    double maxs[n_components] = {250, 950};
    double nots[n_components] = {0.0, 0.0};
    Matrix3D N = get_initial_values(2, mins, maxs, size_x, size_y);
    Matrix3D dN = get_initial_values(2, nots, nots, size_x, size_y);
    
    // For testing of reaction terms.
    add_reaction_terms(dN, N);

    ThomasConstants TC = define_thomas_constants(n_components, D, dx, dy, dt);


    // add_diffusion_terms(dN, N, TC);

    print_Matrix3D(dN);

    return 0;
}


// Implements the tridiagonal thomas algorithm with constant 
// coefficients a, b, c and c_prime
void thomas_algorithm(const int N, const double* a, const double* b, const double* cp, const double* d, double* x) {
    double dp[N];

    for (int i=0; i<N; i++) {
        if (i==0) {
            dp[i] = d[i]/b[i];
        } else {
            dp[i] = (d[i] - a[i]*dp[i])/(b[i]-a[i]*cp[i-1]);
        }
    }
    for (int i=N-1; i>=0; i--) {
        if (i==N-1) {
            x[i] = dp[i];
        } else {
            x[i] = dp[i] - cp[i]*x[i+1];
        }
    }
}


void add_diffusion_terms(Matrix3D &dN, const Matrix3D &N, const ThomasConstants &ct) {
    for (int i=0; i<N.size_x; i++) {
        for (int j=0; j<N.size_y; j++) {
            for (int k=0; k<N.size_z; k++) {
                dN.values[i][j][k] += 0.0;
            }
        }
    }
}


double* calculate_thomas_cp(const int N, const double* a, const double* b, const double* c) {
    double* cp = (double*)malloc(N * sizeof(double));
    for (int i=0; i<N; i++) {
        if (i==0) {
            cp[i] = c[i]/b[i];
        } else {
            cp[i] = c[i]/(b[i]-a[i]*cp[i]);
        }
    }
    return cp;
}


void add_reaction_terms(Matrix3D &dN, const Matrix3D &N) {
    double k1 = 10.0;
    double k2 = 0.1;
    double k3 = 4.938271604938272e-07;
    double k4 = 80.0;

    for (int i=0; i<N.size_y; i++) {
        for (int j=0; j<N.size_z; j++) {
            dN.values[0][i][j] += k1 - k2 * N.values[0][i][j] + k3 * pow(N.values[0][i][j], 2.0) * N.values[1][i][j];
            dN.values[1][i][j] += k4                          - k3 * pow(N.values[0][i][j], 2.0) * N.values[1][i][j];
        }
    }
}


// TODO this is not working yet!
ThomasConstants define_thomas_constants(const int N, const double* D, const double &dx, const double &dy, const double &dt) {

}


// Simple helper function to print a 2D matrix
void print_Matrix2D(Matrix2D M, int dim=1) {
    if (dim) {
        printf("2D Matrix:\n");
        printf("  Size x: %i\n", M.size_x);
        printf("  Size y: %i", M.size_y);
    }
    for (int i=0; i<M.size_x; i++) {
        if (i==0) {
            printf("\n  [[");
        } else {
            printf("\n   [");
        }
        for (int j=0; j<M.size_y; j++) {
            if (j==M.size_y-1) {
                printf("%4.3e]", M.values[i][j]);
            } else {
                printf("%4.3e ", M.values[i][j]);
            }
        }
    }
    printf("]\n");
}


// Simple helper function to print a 3d matrix
void print_Matrix3D(Matrix3D N) {
    printf("[");
    for (int i=0; i<N.size_x; i++) {
        Matrix2D M = {
            N.size_y,
            N.size_z,
            N.values[i]
        };
        print_Matrix2D(M, 0);
        if (i!=N.size_x-1) {
            printf(",");
        } else {
            printf("]\n");
        }
    }
}


Matrix2D get_2D_random_initial_values(double min, double max, int size_x, int size_y) {
    // Check that min and max are well defined
    if (min>max) {
        printf("Minimum has to be smaller than maximum value! Changing min and max!\n");
        int intermediate = min;
        min = max;
        max = min;
    }

    // initialize the values for Matrix
    double** vals;
    vals = (double**) malloc(size_x * sizeof *vals);

    // Allocate and define initial values
    for (int i=0; i<size_x; i++) {
        // Allocate 
        vals[i] = (double*) malloc(size_y * sizeof *vals[i]);
        // define initial values vie random sampling
        for (int j=0; j<size_y; j++) {
            vals[i][j] = (max-min) * (double)rand()/(double)RAND_MAX + min;
        }
    }
    // Store everything in a matrix and return it
    Matrix2D initial_vals = {
        size_x,
        size_y,
        vals
    };
    return initial_vals;
}


Matrix3D get_initial_values(int n_components, double min[], double max[], int size_x, int size_y) {
    // Allocate values and write them
    double*** values;
    values = (double*** ) malloc(n_components * sizeof *values);
    
    // Allocate more and define initial values 
    for (int i=0; i<n_components; i++) {
        Matrix2D D = get_2D_random_initial_values(min[i], max[i], size_x, size_y);
        values[i] = (double**) malloc( sizeof(D.values));
        values[i] = D.values;
    }
    
    // Store everything in a matrix of 3 dimensions and return it
    Matrix3D initials = {
        n_components,
        size_x,
        size_y,
        values
    };
    return initials;
}

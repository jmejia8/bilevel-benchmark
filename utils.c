#include <stdlib.h>

void randm(double a, double b, double *x, int D){
    int i;
    for (i = 0; i < D; ++i) {
        x[i] = a + (double) (rand() / (double) (RAND_MAX) * (b - a));
    }
}


double* array(int D){
    return (double *) malloc(D * sizeof( double ));
}
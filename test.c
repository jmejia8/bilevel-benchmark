#include <stdio.h>
#include <time.h>

#include "blb18_op.c"

void randm(double a, double b, double *x, int D){
    int i;
    for (i = 0; i < D; ++i) {
        x[i] = a + (double) (rand() / (double) (RAND_MAX) * (b - a));
    }
}


double* array(int D){
    return (double *) malloc(D * sizeof( double ));
}

int main(int argc, char const *argv[])
{
    srand(time(NULL));
    int i, id;
 
    // population size
    int N = 5;

    // upper and lower level dimension
    int D = 2;

    // allocate vectors
    double *x = array(N*D);
    double *y = array(N*D);

    // upper level
    double *F = array(N);
    // lower lower
    double *f = array(N);

    // random initialization
    randm(0, 1, x, N*D);
    randm(0, 1, y, N*D);


    printf("run \t i \t F \t \t f \n");
    for (id = 1; id <= 8; ++id) {
        // evaluate
        blb18_leader_cop(N, D, x, y, F, id);
        blb18_follower_cop(N, D, x, y, f, id);

        for (i = 0; i < N; ++i) {
            printf("%i \t %i \t %e \t %e\n", id, i+1, F[i], f[i] );
        }
        
    }

    return 0;
}
#include <stdio.h>
#include <time.h>

#include "blb18_cop.c"

void randm(double a, double b, double *x, int D){
    int i;
    for (i = 0; i < D; ++i) {
        x[i] = a + (double) (rand() / (double) (RAND_MAX) * (b - a));
    }
}

int main(int argc, char const *argv[])
{
    // population size
    int N = 17;

    // upper and lower lever dimension
    int D = 2;

    // number of test function
    int id = 5;

    double *x = array(N*D);
    double *y = array(N*D);

    // upper level
    double *F = array(N);
    // upper lower
    double *f = array(N);

    srand(time(NULL));

    randm(-1, 1, x, N*D);
    randm(-1, 1, y, N*D);

    // evaluate
    blb18_leader_cop(N, D, x, y, F, id);
    blb18_follower_cop(N, D, x, y, f, id);

    int i;
    for (i = 0; i < N; ++i) {
        printf("%i \t %e \t %e\n",i+1, F[i], f[i] );
    }

    return 0;
}
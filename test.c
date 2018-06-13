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

int test(){
    int p=3, q=3, r=2, s=0;
    int D = p + r;

    double *x = array(D);
    double *y = array(D);
    double F[1], Fv = 0;

    int i = 0;

    for (i = 0; i < D; ++i){
            x[i] = 0; y[i] = 0;
    }

    SDM1_follower(p, q, r, x, y, F);
    if (F[0] != 0) {
        printf("1 --> %e\n", F[0]);
        Fv += F[0];
    }
    
    SDM3_follower(p, q, r, x, y, F);
    if (F[0] != 0) {
        printf("3 --> %e\n", F[0]);
        Fv += F[0];
    }

    SDM4_follower(p, q, r, x, y, F);
    if (F[0] != 0) {
        printf("4 --> %e\n", F[0]);
        Fv += F[0];
    }

    SDM6_follower(p, 1, r, 2, x, y, F);
    if (F[0] != 0) {
        printf("6 --> %e\n", F[0]);
        Fv += F[0];
    }

    SDM8_follower(p, q, r, x, y, F);
    if (F[0] != 0) {
        printf("8 --> %e\n", F[0]);
        Fv += F[0];
    }

    for (i = q; i < D; ++i) y[i] = 1;

    SDM2_follower(p, q, r, x, y, F);
    if (F[0] != 0) {
        printf("2 --> %e\n", F[0]);
        Fv += F[0];
    }

    SDM7_follower(p, q, r, x, y, F);
    if (F[0] != 0) {
        printf("7 --> %e\n", F[0]);
        Fv += F[0];
    }


    for (i = 0; i < D; ++i) y[i] = 0;
    for (i = 0; i < q; ++i) y[i] = 1;    
    
    SDM5_follower(p, q, r, x, y, F);
    if (F[0] != 0) {
        printf("5 --> %e\n", F[0]);
        Fv += F[0];
    }

    if (Fv != 0){
        printf("Error\n");
        exit(0);
    }

    return 1;
}

int main(int argc, char const *argv[])
{
    test();
    srand(time(NULL));
    int i, id;
 
    // population size
    int N = 5;

    // upper and lower level dimension
    int D_upper = 10;
    int D_lower = 10;

    // allocate vectors
    double *x = array(N*D_upper);
    double *y = array(N*D_lower);

    // upper level
    double *F = array(N);
    // lower lower
    double *f = array(N);

    // random initialization
    randm(0, 1, x, N*D_upper);
    randm(0, 1, y, N*D_lower);


    printf("Problem  i \t F \t \t f \n");
    for (id = 1; id <= 8; ++id) {
        // evaluate
        blb18_leader_cop(N, D_upper, D_lower, x, y, F, id);
        blb18_follower_cop(N, D_upper, D_lower, x, y, f, id);

        for (i = 0; i < N; ++i) {
            printf("%i \t %i \t %e \t %e\n", id, i+1, F[i], f[i] );
        }
        printf("---------------------------------------------\n");
        
    }

    return 0;
}
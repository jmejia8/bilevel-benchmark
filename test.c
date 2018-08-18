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
    int D_upper = p + r;
    int D_lower = q + r;

    double *x = array(D_upper);
    double *y = array(D_lower);
    double F[1], f[1], Fv = 0;

    F[0] = 1000; f[0] = 1000;

    int i = 0;

    ///////////////////////////////////////////////// 
    ///////////////////////////////////////////////// 
    blb18_cop_solutions(D_upper, D_lower, x, y, 1);

    SMD1_leader(p, q, r, x, y, F);
    SMD1_follower(p, q, r, x, y, f);
    if (F[0] + f[0] != 0) {
        printf("1 --> %e\n", F[0]);
        Fv += F[0];
    }
    
    SMD3_leader(p, q, r, x, y, F);
    SMD3_follower(p, q, r, x, y, f);
    if (F[0] + f[0] != 0) {
        printf("3 --> %e\n", F[0]);
        Fv += F[0];
    }

    SMD4_leader(p, q, r, x, y, F);
    SMD4_follower(p, q, r, x, y, f);
    if (F[0] + f[0] != 0) {
        printf("4 --> %e\n", F[0]);
        Fv += F[0];
    }

    SMD6_leader(p, 1, r, 2, x, y, F);
    SMD6_follower(p, 1, r, 2, x, y, f);
    if (F[0] + f[0] != 0) {
        printf("6 --> %e\n", F[0]);
        Fv += F[0];
    }

    ///////////////////////////////////////////////// 
    ///////////////////////////////////////////////// 
    blb18_cop_solutions(D_upper, D_lower, x, y, 2);

    SMD2_leader(p, q, r, x, y, F);
    SMD2_follower(p, q, r, x, y, f);
    if (F[0] + f[0] != 0) {
        printf("2 --> %e\n", F[0]);
        Fv += F[0];
    }

    SMD7_leader(p, q, r, x, y, F);
    SMD7_follower(p, q, r, x, y, f);
    if (F[0] + f[0] != 0) {
        printf("7 --> %e\n", F[0]);
        Fv += F[0];
    }


    ///////////////////////////////////////////////// 
    ///////////////////////////////////////////////// 
    blb18_cop_solutions(D_upper, D_lower, x, y, 5);
    
    SMD5_leader(p, q, r, x, y, F);
    SMD5_follower(p, q, r, x, y, f);
    if (F[0] + f[0] != 0) {
        printf("5 --> %e\n", F[0]);
        Fv += F[0];
    }


    SMD8_leader(p, q, r, x, y, F);
    SMD8_follower(p, q, r, x, y, f);
    if (F[0] + f[0] != 0) {
        printf("8 --> %e\n", F[0]);
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
        blb18_leader_cop(N, D_upper, D_lower, x, y, F, NULL, id);
        blb18_follower_cop(N, D_upper, D_lower, x, y, f, NULL, id);

        for (i = 0; i < N; ++i) {
            printf("%i \t %i \t %e \t %e\n", id, i+1, F[i], f[i] );
        }
        printf("---------------------------------------------\n");
        
    }

    return 0;
}
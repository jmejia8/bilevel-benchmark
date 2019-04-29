#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "blb18_op.c"
#include "utils.h"

#define DEBUG 1

double PMM_test1(int m, int n){
    int FNUN = 10, fnum;

    double *x = array(n);
    double *y = array(n);
    double F[1], f[1], f_sum = 0;
    double G[2], g[2];

    randm(-1.0, 1, x, 1*n);
    randm(-1, 1.0, y, 1*n);

    for (fnum = 1; fnum <= FNUN; ++fnum){
        PMM_leader(m, n, x, y, F,G, fnum);
        PMM_follower(m, n, x, y, f, g, fnum);

        f_sum += fabs( F[0] + f[0] );
        // printf("PMM%d \t F = %.4e \t f = %.4e\n", fnum, F[0], f[0]);
    } 
    return f_sum;
}

double PMM_test2(int m, int n){
    int FNUN = 10, fnum, i;

    double *x = array(n);
    double *y = array(n);
    double F[1], f[1], f_sum = 0;
    double G[2], g[2];


    for (fnum = 1; fnum <= FNUN; ++fnum){
        randm(-10.0, 10.0, x, 1*n);
        PMM_Psi(n,n,m, x, y, fnum);

        PMM_leader(m, n, x, y, F,G, fnum);
        PMM_follower(m, n, x, y, f, g, fnum);

        f_sum += fabs( f[0]);
        
        if (fnum > 5 && DEBUG) {
            printf("PMM%d \t F = %.4e \t f = %.4e \t G = %.4e \t g = %.4e\n", fnum, F[0], f[0], G[0], g[0]);
        }else if (DEBUG){
            printf("PMM%d \t F = %.4e \t f = %.4e\n", fnum, F[0], f[0]);
        }

    }

    if (DEBUG) printf("--------------------------\n");

    for (i = 0; i < n; ++i) {
        x[i] = 0.0; y[i] = 0.0;
    }

    for (fnum = 1; fnum <= FNUN; ++fnum){

        PMM_leader(m, n, x, y, F,G, fnum);
        PMM_follower(m, n, x, y, f, g, fnum);

        if (fnum <= 5) {
            f_sum += fabs( F[0] + f[0] );
            if (DEBUG) printf("PMM%d \t F = %.4e \t f = %.4e\n", fnum, F[0], f[0]);
        }else{
            double ss = fabs( F[0] + f[0]);
            ss += G[0] > 0 ? G[0]: 0.0;
            ss += g[0] > 0 ? g[0]: 0.0;
            f_sum += ss;
            if (DEBUG) printf("PMM%d \t F = %.4e \t f = %.4e \t G = %.4e \t g = %.4e\n", fnum, F[0], f[0], G[0], g[0]);
        }

    }


    return f_sum;
}

double PMM_test(int m, int n){
    double a = PMM_test2(m, n);
    return a;
}


int test(){
    int settings[6];
    int D_ul = 5,  D_ll = 5;
    int id;

    blb18_cop_settings(D_ul, D_ll, settings, id);

    int p = settings[0], q = settings[1], r = settings[2], s = settings[3];
    int lenG, leng;

    double *x = array(D_ul);
    double *y = array(D_ll);
    double F[1], f[1], f_sum = 0;
    double *G, *g;



    ///////////////////////////////////////////////// 
    /////////////////////////////////////////////////
    for (id = 1; id <= 12; ++id){

        if (id > 8) {
            blb18_cop_settings(D_ul, D_ll, settings, id);
            lenG = settings[4];
            leng = settings[5];
            G = array(lenG); g = array(leng);

        }


        blb18_cop_solutions(D_ul, D_ll, x, y, id);
        blb18_leader_cop(1, D_ul, D_ll, x, y, F, G, id);
        blb18_follower_cop(1, D_ul, D_ll, x, y, f, g, id);

        if (id < 10 && abs(F[0]) + abs(f[0]) > 1e-8) {
            printf("%d --> F = %.4e \t f = %.4e\n", id, F[0], f[0]);
            f_sum += F[0] + f[0];
        }


         
    }

    // SMD1_leader(p, q, r, x, y, F);
    // SMD1_follower(p, q, r, x, y, f);


    if (abs(f_sum) > 1e-8){
        printf("Error\n");
        exit(0);
    }

    return 1;
}

int main(int argc, char const *argv[])
{  
    srand(time(NULL));
    double r = PMM_test(5, 10);
    if (r > 1e-16) {
        printf("Errors in PMM: %g\n", r);
        exit(0);
    }else{
        printf("PMM: OK\n");
    }
    test();
    int i, j, id, settings[6];

 
    // population size
    int N = 2;

    // upper and lower level dimension
    int D_ul = 5;
    int D_ll = 5;
    
    // blb18_cop_settings(D_ul, D_ll, settings, fnum);

    // allocate vectors
    double *x = array(N*D_ul);
    double *y = array(N*D_ll);

    // upper level
    double *F = array(N);
    // lower lower
    double *f = array(N);

    // random initialization
    randm(0, 1, x, N*D_ul);
    randm(0, 1, y, N*D_ll);

    int lenG = 0, leng = 0;
    double *G, *g;

    printf("Problem  i \t F \t \t f  \n");
    for (id = 1; id <= 12; ++id) {
        
        if (id > 8) {
            // get settings
            blb18_cop_settings(D_ul, D_ll, settings, id);

            // free(G);
            // free(g);

            // numbre or restrictions
            lenG = settings[4]; leng = settings[5];
            G = array(lenG*N); g = array(leng*N);
        }

        // evaluate
        blb18_leader_cop(N, D_ul, D_ll, x, y, F, G, id);
        blb18_follower_cop(N, D_ul, D_ll, x, y, f, g, id);

        for (i = 0; i < N; ++i) {
            printf("%i \t %i \t %.2e \t %.2e \n", id, i+1, F[i], f[i] );

            printf("  G = [");
            for (j = 0; j < lenG; ++j) printf("%.2e, ", G[i*lenG + j]);
            printf("]\n  g = [");
            for (j = 0; j < leng; ++j) printf("%.2e, ", g[i*leng + j]);
            printf("]\n");

        }
        printf("---------------------------------------------\n");

        
    }

    return 0;
}
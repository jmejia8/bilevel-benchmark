#include <stdio.h>
#include <time.h>

#include "blb18_op.c"

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
            printf("%d --> F = %e \t f = %e\n", id, F[0], f[0]);
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
    PMM_test(5, 10);
    return 0;
    test();
    srand(time(NULL));
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
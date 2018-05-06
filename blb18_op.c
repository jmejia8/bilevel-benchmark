// #include <WINDOWS.H>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
// #include <iostream>
// using namespace std;


#define INF 1.0e99
#define EPS 1.0e-14
#define E  2.7182818284590452353602874713526625
#define PI 3.1415926535897932384626433832795029

// F (x_u , x_l ) = F_1 (x_u1 ) + F_2 (x_l1 ) + F_3 (x_u2 , x_l2 )
// f (x_u , x_l ) = f_1 (x_u1 , x_u2 ) + f_2 (x_l1 ) + f_3 (x_u2 , x_l2 )
// where
// x_u = x = (x_u1 , x_u2 )
// x_l = y = (x_l1 , x_l2 )

double sphere(double *x, int D){
    int i; double f = 0.0;

    for (i = 0; i < D; ++i) f += x[i] * x[i];

    return f;
}

double* array(int D){
    return (double *) malloc(D * sizeof( double ));
}

void SDM1_leader(int p, int q, int r, double *x, double *y, double *F){
    int i;

    double *x_u1 = x, *x_u2 = &x[p];
    double *x_l1 = y, *x_l2 = &y[q];

    double F1 = sphere(x_u1, p);
    double F2 = sphere(x_l1, q);
    double F3 = 0.0;

    for (i = 0; i < r; ++i) {
        F3 += pow(x_u2[i], 2) + pow( x_u2[i] - tan( x_l2[i] ), 2);
    }

    F[0] = F1 + F2 + F3;
}

void SDM1_follower(int p, int q, int r, double *x, double *y, double *f){
    int i;

    double* x_u1 = x, *x_u2 = &x[p];
    double* x_l1 = y, *x_l2 = &y[q];

    double f1 = sphere(x_u1, p);
    double f2 = sphere(x_l1, q);
    double f3 = 0.0;

    for (i = 0; i < r; ++i)
        f3 += pow( x_u2[i] - tan( x_l2[i] ), 2);

    f[0] = f1 + f2 + f3;
}

void SDM2_leader(int p, int q, int r, double *x, double *y, double *F){
    int i;

    double *x_u1 = x, *x_u2 = &x[p];
    double *x_l1 = y, *x_l2 = &y[q];

    double F1 =  sphere(x_u1, p);
    double F2 = -sphere(x_l1, q);
    double F3 = 0.0;

    for (i = 0; i < r; ++i) {
        F3 += pow(x_u2[i], 2) - pow( x_u2[i] - log( x_l2[i] ), 2);
    }

    F[0] = F1 + F2 + F3;
}

void SDM2_follower(int p, int q, int r, double *x, double *y, double *f){
    int i;

    double* x_u1 = x, *x_u2 = &x[p];
    double* x_l1 = y, *x_l2 = &y[q];

    double f1 = sphere(x_u1, p);
    double f2 = sphere(x_l1, q);
    double f3 = 0.0;

    for (i = 0; i < r; ++i)
        f3 += pow( x_u2[i] - log( x_l2[i] ), 2);

    f[0] = f1 + f2 + f3;
}

void SDM3_leader(int p, int q, int r, double *x, double *y, double *F){
    int i;

    double *x_u1 = x, *x_u2 = &x[p];
    double *x_l1 = y, *x_l2 = &y[q];

    double F1 = sphere(x_u1, p);
    double F2 = sphere(x_l1, q);
    double F3 = 0.0;

    double aux;
    for (i = 0; i < r; ++i) {
        aux = x_u2[i] * x_u2[i];
        F3 += aux + pow( aux - tan( x_l2[i] ), 2);
    }

    F[0] = F1 + F2 + F3;
}

void SDM3_follower(int p, int q, int r, double *x, double *y, double *f){
    int i;

    double* x_u1 = x, *x_u2 = &x[p];
    double* x_l1 = y, *x_l2 = &y[q];

    double f1 = sphere(x_u1, p);
    double f2 = (double) q;
    double f3 = 0.0;

    double pi2 = 2*PI;
    for (i = 0; i < q; ++i)
        f2 += pow( x_l1[i]*x_l1[i] - cos( pi2* x_l1[i] ), 2);

    for (i = 0; i < r; ++i)
        f3 += pow( x_u2[i]*x_u2[i] - tan( x_l2[i] ), 2);

    f[0] = f1 + f2 + f3;
}

void SDM4_leader(int p, int q, int r, double *x, double *y, double *F){
    int i;

    double *x_u1 = x, *x_u2 = &x[p];
    double *x_l1 = y, *x_l2 = &y[q];

    double F1 =  sphere(x_u1, p);
    double F2 = -sphere(x_l1, q);
    double F3 = 0.0;

    for (i = 0; i < r; ++i) {
        F3 += x_u2[i] * x_u2[i] + pow( abs(x_u2[i]) - log( 1.0 + x_l2[i] ), 2);
    }

    F[0] = F1 + F2 + F3;
}

void SDM4_follower(int p, int q, int r, double *x, double *y, double *f){
    int i;

    double* x_u1 = x, *x_u2 = &x[p];
    double* x_l1 = y, *x_l2 = &y[q];

    double f1 = sphere(x_u1, p);
    double f2 = (double) q;
    double f3 = 0.0;

    double pi2 = 2*PI;
    for (i = 0; i < q; ++i)
        f2 += pow( x_l1[i]*x_l1[i] - cos( pi2* x_l1[i] ), 2);

    for (i = 0; i < r; ++i)
        f3 += pow( abs(x_u2[i]) - log( 1.0 + x_l2[i] ), 2);

    f[0] = f1 + f2 + f3;
}

void SDM5_leader(int p, int q, int r, double *x, double *y, double *F){
    int i;

    double *x_u1 = x, *x_u2 = &x[p];
    double *x_l1 = y, *x_l2 = &y[q];

    double F1 =  sphere(x_u1, p);
    double F2 = 0.0;
    double F3 = 0.0;

    for (i = 0; i < q-1; ++i) {
        F2 -= ( x_l1[i+1] - x_l1[i]*x_l1[i] ) + pow(x_l1[i] - 1.0, 2);
    }

    for (i = 0; i < r; ++i) {
        F3 += x_u2[i] * x_u2[i] - pow( abs(x_u2[i]) - pow( x_l2[i], 2 ), 2);
    }

    F[0] = F1 + F2 + F3;
}

void SDM5_follower(int p, int q, int r, double *x, double *y, double *f){
    int i;

    double* x_u1 = x, *x_u2 = &x[p];
    double* x_l1 = y, *x_l2 = &y[q];

    double f1 = sphere(x_u1, p);
    double f2 = 0.0;
    double f3 = 0.0;

    for (i = 0; i < q-1; ++i) {
        f2 += ( x_l1[i+1] - x_l1[i]*x_l1[i] ) + pow(x_l1[i] - 1.0, 2);
    }

    for (i = 0; i < r; ++i)
        f3 += pow( abs(x_u2[i]) - pow( x_l2[i], 2 ), 2);

    f[0] = f1 + f2 + f3;
}

void SDM7_leader(int p, int q, int r, double *x, double *y, double *F){
    int i;

    double *x_u1 = x, *x_u2 = &x[p];
    double *x_l1 = y, *x_l2 = &y[q];

    double F1 = 1.0 + 0.0025*sphere(x_u1, p);
    double F2 = -sphere(x_l1, q);
    double F3 = 0.0;

    double aux = 1.0;
    for (i = 0; i < p; ++i) {
        aux *= cos(x_u1[i] / sqrt(i + 1.0)) ;
    }

    F1 -= aux;

    for (i = 0; i < r; ++i) {
        F3 += x_u2[i] * x_u2[i] - pow( x_u2[i] - log(x_l2[i]), 2);
    }

    F[0] = F1 + F2 + F3;
}

void SDM7_follower(int p, int q, int r, double *x, double *y, double *f){
    int i;

    double* x_u1 = x, *x_u2 = &x[p];
    double* x_l1 = y, *x_l2 = &y[q];

    double f1 = 0.0;
    double f2 = sphere(x_l1, q);
    double f3 = 0.0;

    for (i = 0; i < p; ++i) f1 += pow(x_u1[i], 3);

    for (i = 0; i < r; ++i)
        f3 += pow( x_u2[i] - log( x_l2[i]), 2 );

    f[0] = f1 + f2 + f3;
}

void SDM8_leader(int p, int q, int r, double *x, double *y, double *F){
    int i;

    double *x_u1 = x, *x_u2 = &x[p];
    double *x_l1 = y, *x_l2 = &y[q];

    double F1 = 20.0 + E - 20.0*exp(-0.2 * sqrt((1.0 / (double) p) * sphere(x_u1, p)));
    double F2 = 0.0;
    double F3 = 0.0;

    double aux = 0.0;
    for (i = 0; i < p; ++i) {
        aux += cos(2*PI*x_u1[i]);
    }

    F1 -= exp( (1.0 / (double) p) * aux);

    for (i = 0; i < q-1; ++i) {
        F2 -= x_l1[i+1] - x_l1[i]*x_l1[i] + pow( x_l1[i] - 1.0, 2);
    }

    for (i = 0; i < r; ++i) {
        F3 += x_u2[i] * x_u2[i] - pow( x_u2[i] - pow( x_l2[i], 3), 2);
    }

    F[0] = F1 + F2 + F3;
}

void SDM8_follower(int p, int q, int r, double *x, double *y, double *f){
    int i;

    double* x_u1 = x, *x_u2 = &x[p];
    double* x_l1 = y, *x_l2 = &y[q];

    double f1 = 0.0;
    double f2 = sphere(x_l1, q);
    double f3 = 0.0;

    for (i = 0; i < p; ++i)
        f1 += abs(x_u1[i]);
    
    for (i = 0; i < q-1; ++i) 
        f2 += x_l1[i+1] - x_l1[i]*x_l1[i] + pow( x_l1[i] - 1.0, 2);

    for (i = 0; i < r; ++i)
        f3 += pow( x_u2[i] - pow( x_l2[i], 3), 2);

    f[0] = f1 + f2 + f3;
}

void blb18_leader_cop(int N, int D, double *x, double *y, double *F, int id){
    int i, j, p, q, r;

    if (id >= 1 && id <= 10 ) {
        p = 1; q = 1; r = 1;
        D = 2;
    }

    for (i = 0; i < N; ++i) {
        j = i*D;

        switch(id) {
            case 1:
                SDM1_leader(p, q, r, &x[j], &y[j], &F[i]);
                break;
            case 2:
                SDM2_leader(p, q, r, &x[j], &y[j], &F[i]);
                break;
            case 3:
                SDM3_leader(p, q, r, &x[j], &y[j], &F[i]);
                break;
            case 4:
                SDM4_leader(p, q, r, &x[j], &y[j], &F[i]);
                break;
            case 5:
                SDM5_leader(p, q, r, &x[j], &y[j], &F[i]);
                break;
            case 6:
                SDM5_leader(p, q, r, &x[j], &y[j], &F[i]);
                break;
            case 7:
                SDM7_leader(p, q, r, &x[j], &y[j], &F[i]);
                break;
            case 8:
                SDM8_leader(p, q, r, &x[j], &y[j], &F[i]);
                break;
            // case 9:
            //     SDM9_leader(p, q, r, &x[j], &y[j], &F[i]);
            //     break;
            // case 10:
            //     SDM10_leader(p, q, r, &x[j], &y[j], &F[i]);
            //     break;
            default:
                printf("Error\n");
                break;
        }
    }
}

void blb18_follower_cop(int N, int D, double *x, double *y, double *f, int id){
    int i, j, p, q, r;

    if (id >= 1 && id <= 10 ) {
        p = 1; q = 1; r = 1;
    }

    for (i = 0; i < N; ++i) {
        j = i*D;

        switch(id) {
            case 1:
                SDM1_follower(p, q, r, &x[j], &y[j], &f[i]);
                break;
            case 2:
                SDM2_follower(p, q, r, &x[j], &y[j], &f[i]);
                break;
            case 3:
                SDM3_follower(p, q, r, &x[j], &y[j], &f[i]);
                break;
            case 4:
                SDM4_follower(p, q, r, &x[j], &y[j], &f[i]);
                break;
            case 5:
                SDM5_follower(p, q, r, &x[j], &y[j], &f[i]);
                break;
            case 6:
                SDM5_follower(p, q, r, &x[j], &y[j], &f[i]);
                break;
            case 7:
                SDM7_follower(p, q, r, &x[j], &y[j], &f[i]);
                break;
            case 8:
                SDM8_follower(p, q, r, &x[j], &y[j], &f[i]);
                break;
            // case 9:
            //     SDM9_follower(p, q, r, &x[j], &y[j], &f[i]);
            //     break;
            // case 10:
            //     SDM10_follower(p, q, r, &x[j], &y[j], &f[i]);
            //     break;
            default:
                printf("Error\n");
                break;
        }
    }
}
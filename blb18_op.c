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

void SMD1_leader(int p, int q, int r, double *x, double *y, double *F){
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

void SMD1_follower(int p, int q, int r, double *x, double *y, double *f){
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

void SMD2_leader(int p, int q, int r, double *x, double *y, double *F){
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

void SMD2_follower(int p, int q, int r, double *x, double *y, double *f){
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

void SMD3_leader(int p, int q, int r, double *x, double *y, double *F){
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

void SMD3_follower(int p, int q, int r, double *x, double *y, double *f){
    int i;

    double* x_u1 = x, *x_u2 = &x[p];
    double* x_l1 = y, *x_l2 = &y[q];

    double f1 = sphere(x_u1, p);
    double f2 = (double) q;
    double f3 = 0.0;

    double pi2 = 2*PI;
    for (i = 0; i < q; ++i)
        f2 += x_l1[i]*x_l1[i] - cos( pi2* x_l1[i] );

    for (i = 0; i < r; ++i)
        f3 += pow( x_u2[i]*x_u2[i] - tan( x_l2[i] ), 2);

    f[0] = f1 + f2 + f3;
}

void SMD4_leader(int p, int q, int r, double *x, double *y, double *F){
    int i;

    double *x_u1 = x, *x_u2 = &x[p];
    double *x_l1 = y, *x_l2 = &y[q];

    double F1 =  sphere(x_u1, p);
    double F2 = -sphere(x_l1, q);
    double F3 =  sphere(x_u2, r);

    for (i = 0; i < r; ++i) {
        F3 -= pow( fabs(x_u2[i]) - log( 1.0 + x_l2[i] ), 2);
    }

    F[0] = F1 + F2 + F3;
}

void SMD4_follower(int p, int q, int r, double *x, double *y, double *f){
    int i;

    double* x_u1 = x, *x_u2 = &x[p];
    double* x_l1 = y, *x_l2 = &y[q];

    double f1 = sphere(x_u1, p);
    double f2 = (double) q;
    double f3 = 0.0;

    double pi2 = 2*PI;
    for (i = 0; i < q; ++i)
        f2 += x_l1[i]*x_l1[i] - cos( pi2* x_l1[i] );

    for (i = 0; i < r; ++i)
        f3 += pow( fabs(x_u2[i]) - log( 1.0 + x_l2[i] ), 2);

    f[0] = f1 + f2 + f3;
}

void SMD5_leader(int p, int q, int r, double *x, double *y, double *F){
    int i;

    double *x_u1 = x, *x_u2 = &x[p];
    double *x_l1 = y, *x_l2 = &y[q];

    double F1 =  sphere(x_u1, p);
    double F2 = 0.0;
    double F3 = 0.0;

    for (i = 0; i < q-1; ++i) {
        F2 -= pow( x_l1[i+1] - x_l1[i]*x_l1[i], 2) + pow(x_l1[i] - 1.0, 2);
    }

    for (i = 0; i < r; ++i) {
        F3 += x_u2[i] * x_u2[i] - pow( fabs(x_u2[i]) - pow( x_l2[i], 2 ), 2);
    }

    F[0] = F1 + F2 + F3;
}

void SMD5_follower(int p, int q, int r, double *x, double *y, double *f){
    int i;

    double* x_u1 = x, *x_u2 = &x[p];
    double* x_l1 = y, *x_l2 = &y[q];

    double f1 = sphere(x_u1, p);
    double f2 = 0.0;
    double f3 = 0.0;

    for (i = 0; i < q-1; ++i) {
        f2 += pow( x_l1[i+1] - x_l1[i]*x_l1[i], 2) + pow(x_l1[i] - 1.0, 2);
    }

    for (i = 0; i < r; ++i)
        f3 += pow( fabs(x_u2[i]) - pow( x_l2[i], 2 ), 2);

    f[0] = f1 + f2 + f3;
}

void SMD6_leader(int p, int q, int r, int s, double *x, double *y, double *F){
    int i;

    double *x_u1 = x, *x_u2 = &x[p];
    double *x_l1 = y, *x_l2 = &y[q+s];

    double F1 = sphere(x_u1, p);
    double F2 = -sphere(x_l1, q) + sphere(&x_l1[q], s) ;
    double F3 = sphere(x_u2, r);

    for (i = 0; i < r; ++i) {
        F3 -= pow( x_u2[i] - x_l2[i], 2);
    }

    F[0] = F1 + F2 + F3;
}

void SMD6_follower(int p, int q, int r, int s, double *x, double *y, double *f){
    int i;

    double* x_u1 = x, *x_u2 = &x[p];
    double* x_l1 = y, *x_l2 = &y[q+s];

    double f1 = sphere(x_u1, p);
    double f2 = sphere(x_l1, q);
    double f3 = 0.0;

    for (i = q; i < q+s-1; i += 2) {
        f2 += pow(x_l1[i+1] - x_l1[i], 2);
    }

    for (i = 0; i < r; ++i)
        f3 += pow( x_u2[i] - x_l2[i], 2);

    f[0] = f1 + f2 + f3;
}

void SMD7_leader(int p, int q, int r, double *x, double *y, double *F){
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

void SMD7_follower(int p, int q, int r, double *x, double *y, double *f){
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

void SMD8_leader(int p, int q, int r, double *x, double *y, double *F){
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
        F2 -= pow(x_l1[i+1] - x_l1[i]*x_l1[i], 2) + pow( x_l1[i] - 1.0, 2);
    }

    for (i = 0; i < r; ++i) {
        F3 += x_u2[i] * x_u2[i] - pow( x_u2[i] - pow( x_l2[i], 3), 2);
    }

    F[0] = F1 + F2 + F3;
}

void SMD8_follower(int p, int q, int r, double *x, double *y, double *f){
    int i;

    double* x_u1 = x, *x_u2 = &x[p];
    double* x_l1 = y, *x_l2 = &y[q];

    double f1 = 0.0;
    double f2 = 0.0;
    double f3 = 0.0;

    for (i = 0; i < p; ++i)
        f1 += fabs(x_u1[i]);
    
    for (i = 0; i < q-1; ++i) 
        f2 += pow(x_l1[i+1] - x_l1[i]*x_l1[i], 2) + pow( x_l1[i] - 1.0, 2);

    for (i = 0; i < r; ++i)
        f3 += pow( x_u2[i] - pow( x_l2[i], 3), 2);

    f[0] = f1 + f2 + f3;
}

void blb18_leader_cop(int N, int D_upper, int D_lower, double *x, double *y, double *F, int id){
    int i, u, l, p, q, r, s;

    r = D_upper / 2;
    p = D_upper - r;
    
    if (id == 6) {
        q = (int) floor( (D_lower - r) / 2);
        s = (int)  ceil( EPS + (D_lower - r) / 2);
    }else{
        r = D_upper / 2;
        q = D_lower - r;
    }

    printf("%d %d %d %d\n",p,q,r,s );

    for (i = 0; i < N; ++i) {
        u = i*D_upper;
        l = i*D_lower;

        switch(id) {
            case 1:
                SMD1_leader(p, q, r, &x[u], &y[l], &F[i]);
                break;
            case 2:
                SMD2_leader(p, q, r, &x[u], &y[l], &F[i]);
                break;
            case 3:
                SMD3_leader(p, q, r, &x[u], &y[l], &F[i]);
                break;
            case 4:
                SMD4_leader(p, q, r, &x[u], &y[l], &F[i]);
                break;
            case 5:
                SMD5_leader(p, q, r, &x[u], &y[l], &F[i]);
                break;
            case 6:
                SMD6_leader(p, q, r, s, &x[u], &y[l], &F[i]);
                break;
            case 7:
                SMD7_leader(p, q, r, &x[u], &y[l], &F[i]);
                break;
            case 8:
                SMD8_leader(p, q, r, &x[u], &y[l], &F[i]);
                break;
            default:
                printf("Error\n");
                break;
        }
    }
}

void blb18_follower_cop(int N, int D_upper, int D_lower, double *x, double *y, double *f, int id){
    int i, u, l, p, q, r, s;

    r = D_upper / 2;
    p = D_upper - r;
    
    if (id == 6) {
        q = (int) floor( (D_lower - r) / 2);
        s = (int)  ceil( EPS + (D_lower - r) / 2);
    }else{
        r = D_upper / 2;
        q = D_lower - r;
    }

    for (i = 0; i < N; ++i) {
        u = i*D_upper;
        l = i*D_lower;

        switch(id) {
            case 1:
                SMD1_follower(p, q, r, &x[u], &y[l], &f[i]);
                break;
            case 2:
                SMD2_follower(p, q, r, &x[u], &y[l], &f[i]);
                break;
            case 3:
                SMD3_follower(p, q, r, &x[u], &y[l], &f[i]);
                break;
            case 4:
                SMD4_follower(p, q, r, &x[u], &y[l], &f[i]);
                break;
            case 5:
                SMD5_follower(p, q, r, &x[u], &y[l], &f[i]);
                break;
            case 6:
                SMD6_follower(p, q/2, r, q/2, &x[u], &y[l], &f[i]);
                break;
            case 7:
                SMD7_follower(p, q, r, &x[u], &y[l], &f[i]);
                break;
            case 8:
                SMD8_follower(p, q, r, &x[u], &y[l], &f[i]);
                break;
            default:
                printf("Error\n");
                break;
        }
    }
}
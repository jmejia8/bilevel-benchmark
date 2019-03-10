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

#include "smd.c"

void blb18_cop_settings(int D_ul, int D_ll, int *settings, int fnum){
    int p, q, r, s;

    r = D_ul / 2;
    p = D_ul - r;
    
    if (fnum == 6) {
        q = (int) floor( (D_ll - r) / 2);
        s = (int)  ceil( EPS + (double) (D_ll - r) / (double) 2.0);
    }else{
        r = D_ul / 2;
        q = D_ll - r;
    }

    settings[0] = p; settings[1] = q; settings[2] = r; settings[3] = s;

    if (fnum == 9){   
        settings[4] = 1;
        settings[5] = 1;
    } else if (fnum == 10){
        settings[4] = p+r;
        settings[5] = q;
    } else if (fnum == 11){
        settings[4] = r;
        settings[5] = 1;
    } else if (fnum == 12){
        settings[4] = 2*r + p;
        settings[5] = q+1;
    }else{
        settings[4] = 0;
        settings[5] = 0;
    }
}

void blb18_leader_cop(int N, int D_ul, int D_ll, double *x, double *y, double *F, double *G, int id){
    int i, u, l, settings[6];

    blb18_cop_settings(D_ul, D_ll, settings, id);

    int p = settings[0], q = settings[1], r = settings[2], s = settings[3];
    int lenG = settings[4]; 

    for (i = 0; i < N; ++i) {
        u = i*D_ul;
        l = i*D_ll;

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
            case 9:
                SMD9_leader(p, q, r, &x[u], &y[l], &F[i], &G[i*lenG]);
                break;
            case 10:
                SMD10_leader(p, q, r, &x[u], &y[l], &F[i], &G[i*lenG]);
                break;
            case 11:
                SMD11_leader(p, q, r, &x[u], &y[l], &F[i], &G[i*lenG]);
                break;
            case 12:
                SMD12_leader(p, q, r, &x[u], &y[l], &F[i], &G[i*lenG]);
                break;
            default:
                printf("Error: Invalid function id (1,2,...,12).\n");
                break;
        }
    }
}

void blb18_follower_cop(int N, int D_ul, int D_ll, double *x, double *y, double *f, double *g, int id){
    int i, u, l, settings[6];

    blb18_cop_settings(D_ul, D_ll, settings, id);

    int p = settings[0], q = settings[1], r = settings[2], s = settings[3];
    int leng = settings[5];

    for (i = 0; i < N; ++i) {
        u = i*D_ul;
        l = i*D_ll;

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
                SMD6_follower(p, q, r, s, &x[u], &y[l], &f[i]);
                break;
            case 7:
                SMD7_follower(p, q, r, &x[u], &y[l], &f[i]);
                break;
            case 8:
                SMD8_follower(p, q, r, &x[u], &y[l], &f[i]);
                break;
            case 9:
                SMD9_follower(p, q, r, &x[u], &y[l], &f[i], &g[i*leng]);
                break;
            case 10:
                SMD10_follower(p, q, r, &x[u], &y[l], &f[i], &g[i*leng]);
                break;
            case 11:
                SMD11_follower(p, q, r, &x[u], &y[l], &f[i], &g[i*leng]);
                break;
            case 12:
                SMD12_follower(p, q, r, &x[u], &y[l], &f[i], &g[i*leng]);
                break;
            default:
                printf("Error: Invalid function id (1,2,...,12).\n");
                break;
        }
    }
}

void blb18_cop_ranges(int D_ul, int D_ll, double *bounds_ul, double *bounds_ll, int fnum){
    int settings[6];

    blb18_cop_settings(D_ul, D_ll, settings, fnum);
    
    int p = settings[0], q = settings[1], r = settings[2], s = settings[3];


    int i;
    double ul1_a, ul1_b, ul2_a, ul2_b, ll1_a, ll1_b, ll2_a, ll2_b;
    if (fnum == 1 || fnum == 3){
        ul1_a = -5.0; ul1_b = 10.0;
        ul2_a = -5.0; ul2_b = 10.0;

        ll1_a = -5.0;  ll1_b = 10.0;
        ll2_a = -PI/2 + EPS; ll2_b = PI/2 - EPS;
    } else if (fnum == 2 || fnum == 7){
        ul1_a = -5.0; ul1_b = 10.0;
        ul2_a = -5.0; ul2_b = 1.0;

        ll1_a = -5.0;  ll1_b = 10.0;
        ll2_a =  EPS; ll2_b = E ;   
    } else if (fnum == 4){   
        ul1_a = -5.0; ul1_b = 10.0;
        ul2_a = -1.0; ul2_b = 1.0;

        ll1_a = -5.0;  ll1_b = 10.0;
        ll2_a =  0; ll2_b = E ;
    } else if (fnum == 5  || fnum == 6 || fnum == 8){   
        ul1_a = -5.0; ul1_b = 10.0;
        ul2_a = -5.0; ul2_b = 10.0;

        ll1_a = -5.0; ll1_b = 10.0;
        ll2_a = -5.0; ll2_b = 10.0;
    } else if (fnum == 9) {
        ul1_a = -5.0; ul1_b = 10.0;
        ul2_a = -5.0; ul2_b = 1.0;

        ll1_a = -5.0; ll1_b = 10.0;
        ll2_a = -1.0 + EPS; ll2_b = -1.0 + E;
    } else if (fnum == 10) {
        ul1_a = -5.0; ul1_b = 10.0;
        ul2_a = -5.0; ul2_b = 10.0;

        ll1_a = -5.0; ll1_b = 10.0;
        ll2_a = -PI / 2.0 + EPS; ll2_b = PI / 2.0 - EPS;
    } else if (fnum == 11) {
        ul1_a = -5.0; ul1_b = 10.0;
        ul2_a = -1.0; ul2_b = 1.0;

        ll1_a = -5.0; ll1_b = 10.0;
        ll2_a = 1.0/E; ll2_b = E;
    } else if (fnum == 12) {
        ul1_a = -5.0; ul1_b = 10.0;
        ul2_a = -14.10; ul2_b = 14.10;

        ll1_a = -5.0; ll1_b = 10.0;
        ll2_a = -1.5 + EPS; ll2_b = 1.5 - EPS;
    }

    for (i = 0; i < p; ++i) {
        bounds_ul[i] = ul1_a;
        bounds_ul[D_ul + i] = ul1_b;
    }
    
    for (i = p; i < p + r; ++i){
        bounds_ul[i] = ul2_a;
        bounds_ul[D_ul + i] = ul2_b;
    }

    for (i = 0; i < q; ++i){
        bounds_ll[i] = ll1_a;
        bounds_ll[D_ul + i] = ll1_b;
    }

    for (i = q; i < q+r+s; ++i){
        bounds_ll[i] = ll2_a;
        bounds_ll[D_ul + i] = ll2_b;
    }


}

void blb18_cop_solutions(int D_ul, int D_ll, double *x, double *y, int fnum){

    int settings[6], i;
    blb18_cop_settings(D_ul, D_ll, settings, fnum);
    int p = settings[0], q = settings[1], r = settings[2], s = settings[3];

    for (i = 0; i < D_ul; ++i)x[i] = 0;
    for (i = 0; i < D_ll; ++i) y[i] = 0;

    if (fnum == 2 || fnum == 7){
        for (i = q; i < D_ll; ++i) y[i] = 1;
    
    } else if (fnum == 5 || fnum == 8) {
        for (i = 0; i < D_ll; ++i) y[i] = 0;
        for (i = 0; i < q; ++i) y[i] = 1; 
    
    } else if (fnum == 10) {
        for (i = 0; i < D_ul; ++i) x[i] = 1.0 / sqrt(p + r - 1);
        for (i = 0; i < q; ++i) y[i] = 1.0 / sqrt(q-1);
        for (i = 0; i < r; ++i) y[q+i] = atan(x[p + i]);
    
    } else if (fnum == 11) {
        for (i = 0; i < r; ++i) y[q + i] = exp(-1.0 / sqrt(r));
    
    } else if (fnum == 12) {
        for (i = 0; i < D_ul; ++i) x[i] = 1.0 / sqrt(p+r-1);
        for (i = 0; i < q; ++i) y[i] = 1.0 / sqrt(q-1);
        for (i = 0; i < r; ++i) y[q+i] = atan( (1.0 / sqrt(p+r-1)) - (1.0 / sqrt(r)));
        
    }


}

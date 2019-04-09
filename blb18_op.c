#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#include "constants.h"
#include "smd.c"
#include "pmm.c"
#include "tp.c"

void blb18_cop_settings(int D_ul, int D_ll, int *settings, int fnum){
    int p, q, r, s;

    r = D_ul / 2;
    p = D_ul - r;
    
    if (fnum == 6) {
        q = (int) floor( - EPS + (D_ll - r) /  (double) 2.0 );
        s = (int)  ceil( EPS + (double) (D_ll - r) / (double) 2.0);
    }else{
        r = D_ul / 2;
        q = D_ll - r;
        s = 0;
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
    SMD_ranges(D_ul, D_ll, bounds_ul, bounds_ll, fnum);
}

void blb18_cop_solutions(int D_ul, int D_ll, double *x, double *y, int fnum){

    int settings[6], i;
    blb18_cop_settings(D_ul, D_ll, settings, fnum);
    int p = settings[0], q = settings[1], r = settings[2];

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

void TP_leader(int D_ul, int D_ll, double *x, double *y, double *F, double *G, int fnum){
    if ( fnum == 1 ) {
        TP1_leader(D_ul, D_ll, x, y, F, G);
    }else if (fnum == 2) {
        TP2_leader(D_ul, D_ll, x, y, F, G);
    }else if (fnum == 3) {
        TP3_leader(D_ul, D_ll, x, y, F, G);
    }else if (fnum == 4) {
        TP4_leader(D_ul, D_ll, x, y, F, G);
    }else if (fnum == 5) {
        TP5_leader(D_ul, D_ll, x, y, F, G);
    }else if (fnum == 6) {
        TP6_leader(D_ul, D_ll, x, y, F, G);
    }else if (fnum == 7) {
        TP7_leader(D_ul, D_ll, x, y, F, G);
    }else if (fnum == 8) {
        TP8_leader(D_ul, D_ll, x, y, F, G);
    }else if (fnum == 9) {
        TP9_leader(D_ul, D_ll, x, y, F);
    }else if (fnum == 10) {
        TP10_leader(D_ul, D_ll, x, y, F);
    }
}

void TP_follower(int D_ul, int D_ll, double *x, double *y, double *f, double *g, int fnum){
    if ( fnum == 1 ) {
        TP1_follower(D_ul, D_ll, x, y, f);
    }else if (fnum == 2) {
        TP2_follower(D_ul, D_ll, x, y, f, g);
    }else if (fnum == 3) {
        TP3_follower(D_ul, D_ll, x, y, f, g);
    }else if (fnum == 4) {
        TP4_follower(D_ul, D_ll, x, y, f, g);
    }else if (fnum == 5) {
        TP5_follower(D_ul, D_ll, x, y, f, g);
    }else if (fnum == 6) {
        TP6_follower(D_ul, D_ll, x, y, f, g);
    }else if (fnum == 7) {
        TP7_follower(D_ul, D_ll, x, y, f, g);
    }else if (fnum == 8) {
        TP8_follower(D_ul, D_ll, x, y, f, g);
    }else if (fnum == 9) {
        TP9_follower(D_ul, D_ll, x, y, f);
    }else if (fnum == 10) {
        TP10_follower(D_ul, D_ll, x, y, f);
    }
}

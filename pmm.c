#include "constants.h"

void PMM_config(int D_ul, int D_ll, int *settings, int fnum){
    int m = D_ul < D_ll? D_ul : D_ll;
    m /= 2;

    settings[0] = D_ul;
    settings[1] = D_ll;
    settings[2] = m;

    if ( 6 <= fnum && fnum <= 8) {
        settings[3] = 1; settings[4] = 1;
    }else if (fnum == 9 || fnum == 10) {
        settings[3] = 2; settings[4] = 2;
    }else{
        settings[3] = 0; settings[4] = 0;
    }

}

void PMM_Psi(int D_ul, int D_ll, int m, double *x, double *y, int fnum){
    int i = 1;

    if (fnum == 1 || fnum == 6){
       for (i = 0; i < D_ll; ++i) { y[i] = i < m ? 0.01*pow(x[i], 3) : 0.0; }
    }else if (fnum == 2) {
       for (i = 0; i < D_ll; ++i) { y[i] = i < m ? x[i]*sin(x[i]) : 0.0; }
    }else if (fnum == 3) {
       for (i = 0; i < D_ll; ++i) { y[i] = i < m ? (10.0/(1.0 + 2.5*x[i]*x[i])) - 10.0 : 0.0; }
    }else if (fnum == 4) {
       for (i = 0; i < D_ll; ++i) { y[i] = i < m ? 0.01*pow(x[i], 3) + sin(2.0*PI*x[i]) : 0.0; }
    }else if (fnum == 5 || fnum == 8 || fnum == 9) {
       for (i = 0; i < D_ll; ++i) { y[i] = i < m ? x[i] : 0.0; }
    }else if (fnum == 7) {
       for (i = 0; i < D_ll; ++i) { y[i] = i < m ? 10.0*exp(-0.01*x[i]*x[i])*sin(x[i]) : 0.0; }
    }else if (fnum == 10) {
       for (i = 0; i < D_ll; ++i) { y[i] = i < m ? -x[i]*fabs(sin(PI*x[i])) : 0.0; }
    }
}

void PMM1_leader(int D_ul, int D_ll, int k, double *x, double *y, double *F){
    int i;

    double p1  = 0.0, p2  = 0.0, q = 0.0, r = 0.0;

    r += x[0]*x[0];
    for (i = 1; i < k; ++i){
        q += y[i]*y[i];
        r += x[i]*x[i];
    }


    q = 10.0 + y[0] + 1e6*q;

    p1 = q  - r;

    q  = r = 0.0;
    for (i = k; i < D_ul; ++i){
        q += pow(y[i] - (pow(x[i], 3) / 100.0), 2);
        r += x[i]*x[i];
    }

    p2 = q + r;

    double Q = fabs(p1);
    F[0] = Q + p2;
}

void PMM1_follower(int D_ul, int D_ll, int k, double *x, double *y, double *f){
    int i;

    double p1  = 0.0, p2  = 0.0, q = 0.0, r = 0.0;

    r += x[0]*x[0];
    for (i = 1; i < k; ++i){
        q += y[i]*y[i];
        r += x[i]*x[i];
    }


    q = 10.0 + y[0] + 1e6*q;

    p1 = q  - r;

    q  = r = 0.0;
    for (i = k; i < D_ul; ++i){
        q += pow(y[i] - (pow(x[i], 3) / 100.0), 2);
        r += x[i]*x[i];
    }

    p2 = q + r;

    f[0] = p1 + p2;
}


void PMM2_leader(int D_ul, int D_ll, int k, double *x, double *y, double *F){
    int i;

    double p1  = 0.0, p2  = 0.0, q = 0.0, r = 0.0;

    r += x[0]*x[0];
    for (i = 1; i < k; ++i){
        q += y[i]*y[i];
        r += x[i]*x[i];
    }


    q = 10.0 + y[0] + 1e6*q;

    p1 = q  - r;

    q  = r = 0.0;
    for (i = k; i < D_ul; ++i){
        q += pow(y[i] - 0.01*(pow(x[i], 3)), 2);
        r += x[i]*x[i];
    }

    p2 = q - r;

    F[0] = abs(p1) + abs(p2);
}

void PMM2_follower(int D_ul, int D_ll, int k, double *x, double *y, double *f){
    int i;

    double p1  = 0.0, p2  = 0.0, q = 0.0, r = 0.0;

    r += x[0]*x[0];
    for (i = 1; i < k; ++i){
        q += y[i]*y[i];
        r += x[i]*x[i];
    }


    q = 10.0 + y[0] + 1e6*q;

    p1 = q  - r;

    q  = r = 0.0;
    for (i = k; i < D_ul; ++i){
        q += pow(y[i] - (pow(x[i], 3) / 100.0), 2);
        r += x[i]*x[i];
    }

    p2 = q + r;

    f[0] = p1 + p2;
}


void PMM_leader(int D_ul, int D_ll, double *x, double *y, double *F, double *G, int fnum){
    int m = D_ul < D_ll? D_ul : D_ll;
    m /= 2;
    if ( fnum == 1 ) {
        PMM1_leader(D_ul, D_ll, m, x, y, F);
    }else if (fnum == 2) {
        PMM2_leader(D_ul, D_ll, m, x, y, F);
    }else if (fnum == 3) {
        PMM3_leader(D_ul, D_ll, m, x, y, F);
    }else if (fnum == 4) {
        PMM4_leader(D_ul, D_ll, m, x, y, F);
    }else if (fnum == 5) {
        PMM5_leader(D_ul, D_ll, m, x, y, F);
    }else if (fnum == 6) {
        PMM6_leader(D_ul, D_ll, m, x, y, F, G);
    }else if (fnum == 7) {
        PMM7_leader(D_ul, D_ll, m, x, y, F, G);
    }else if (fnum == 8) {
        PMM8_leader(D_ul, D_ll, m, x, y, F, G);
    }else if (fnum == 9) {
        PMM9_leader(D_ul, D_ll, m, x, y, F, G);
    }else if (fnum == 10) {
        PMM10_leader(D_ul, D_ll, m, x, y, F, G);
    }
}

void PMM_follower(int D_ul, int D_ll, double *x, double *y, double *f, double *g, int fnum){
    int m = D_ul < D_ll? D_ul : D_ll;
    m /= 2;
    if ( fnum == 1 ) {
        PMM1_follower(D_ul, D_ll, m, x, y, f);
    }else if (fnum == 2) {
        PMM2_follower(D_ul, D_ll, m, x, y, f);
    }else if (fnum == 3) {
        PMM3_follower(D_ul, D_ll, m, x, y, f);
    }else if (fnum == 4) {
        PMM4_follower(D_ul, D_ll, m, x, y, f);
    }else if (fnum == 5) {
        PMM5_follower(D_ul, D_ll, m, x, y, f);
    }else if (fnum == 6) {
        PMM6_follower(D_ul, D_ll, m, x, y, f, g);
    }else if (fnum == 7) {
        PMM7_follower(D_ul, D_ll, m, x, y, f, g);
    }else if (fnum == 8) {
        PMM8_follower(D_ul, D_ll, m, x, y, f, g);
    }else if (fnum == 9) {
        PMM9_follower(D_ul, D_ll, m, x, y, f, g);
    }else if (fnum == 10) {
        PMM10_follower(D_ul, D_ll, m, x, y, f, g);
    }
}


#include "constants.h"


/*
 * save in `settings` the following values: upper level dimension, lower level dimension and `k`.
*/
void PMM_config(int D_ul, int D_ll, int *settings, int fnum){
    int k = D_ul < D_ll? D_ul : D_ll;
    k /= 2;

    settings[0] = D_ul;
    settings[1] = D_ll;
    settings[2] = k;

    if ( 7 <= fnum && fnum <= 10) {
        settings[3] = 1; settings[4] = 1;
    }else if (fnum == 9 || fnum == 10) {
        settings[3] = 2; settings[4] = 2;
    }else{
        settings[3] = 0; settings[4] = 0;
    }

}

/*
 * D_ul: upper level dimension
 * D_ll: lower level dimension
 * x: upper level decision vector
 * fnum: function number, i.e., PMM<fnum>
 * Save in y the lower level optimal solution, i.e., y in argmin f(x, y)
*/
void PMM_Psi(int D_ul, int D_ll, int k, double *x, double *y, int fnum){
    int i = 1;

    if (fnum == 1 ){
        for (i = 0; i < D_ll; ++i) { y[i] = i < k ? 0.01*pow(x[i], 3) : 0.0; }
        y[k] = -10.0; 
    }else if (fnum == 2) {
       for (i = 0; i < D_ll; ++i) { y[i] = i < k ? x[i]*sin(x[i]) : sqrt((double) i+1 ); }
    }else if (fnum == 3) {
        for (i = 0; i < D_ll; ++i) { y[i] = i < k ? 0.01*pow(x[i], 3) : 1.0; }
    }else if (fnum == 4) {
       for (i = 0; i < D_ll; ++i) { y[i] = i < k ? (10.0/(1.0 + 2.5*x[i]*x[i])) : 0.0; }
    }else if (fnum == 5 || fnum == 6) {
       for (i = 0; i < D_ll; ++i) { y[i] = i < k ? x[i] : 0.0; }
    }else if (fnum == 7) {
       for (i = 0; i < D_ll; ++i) { y[i] = i < k ? 10.0*exp(-0.01*x[i]*x[i])*sin(x[i]) : 0.0; }
    }else if (fnum == 10) {
       for (i = 0; i < D_ll; ++i) { y[i] = i < k ? -x[i]*fabs(sin(PI*x[i])) : 0.0; }
    }
}

void PMM1_leader(int D_ul, int D_ll, int k, double *x, double *y, double *F){
    int i;

    double p1  = 0.0, p2  = 0.0, q = 0.0, r = 0.0;

    for (i = 0; i < k; ++i){
        q += pow(y[i] - (pow(x[i], 3) / 100.0), 2);
        r += x[i]*x[i];
    }


    p1 = q  - r;

    q  = r = 0.0;
    for (i = k+1; i < D_ll; ++i)
        q += y[i]*y[i];
    
    for (i = k; i < D_ul; ++i)
        r += x[i]*x[i];
 
    q = 10.0 + y[k] + 1e6*q;

    p2 = q - r;

    F[0] = fabs(p1) + fabs(p2);
}

void PMM1_follower(int D_ul, int D_ll, int k, double *x, double *y, double *f){
    int i;

    double p1  = 0.0, p2  = 0.0, q = 0.0, r = 0.0;

    for (i = 0; i < k; ++i){
        q += pow(y[i] - (pow(x[i], 3) / 100.0), 2);
        r += x[i]*x[i];
    }


    p1 = q  - r;

    q  = r = 0.0;
    for (i = k+1; i < D_ll; ++i)
        q += y[i]*y[i];
    
    for (i = k; i < D_ul; ++i)
        r += x[i]*x[i];
 
    q = 10.0 + y[k] + 1e6*q;

    p2 = q - r;

    f[0] = p1 + p2;
}


void PMM2_leader(int D_ul, int D_ll, int k, double *x, double *y, double *F){
    int i;

    double p1  = 0.0, p2  = 0.0, q = 0.0, r = 0.0;

    double a = 0.0, s = 0.0;

    q = pow(y[0] - x[0]*sin(x[0]), 2);
    r += x[0]*x[0];
    for (i = 1; i < k; ++i){
        q += 1e6 * pow(y[i] - x[i]*sin(x[i]), 2);
        r += x[i]*x[i];
    }


    p1 = q  - r;

    q = r = 0.0;
    for (i = k; i < D_ll; ++i){
        a = y[i] - sqrt( (double) i + 1.0 );
        q += a*a;
        s += a;
    }

    q += pow(0.5*s, 2) + pow(0.5*s, 4);

    r = x[k]*x[k];
    for (i = k+1; i < D_ul; ++i){
        r += 1e6*x[i]*x[i];
    }

    p2 = q - r;

    F[0] = fabs(p1) + fabs(p2);
}

void PMM2_follower(int D_ul, int D_ll, int k, double *x, double *y, double *f){
    int i;

    double p1  = 0.0, p2  = 0.0, q = 0.0, r = 0.0;

    double a = 0.0, s = 0.0;

    q = pow(y[0] - x[0]*sin(x[0]), 2);
    r += x[0]*x[0];
    for (i = 1; i < k; ++i){
        q += 1e6 * pow(y[i] - x[i]*sin(x[i]), 2);
        r += x[i]*x[i];
    }


    p1 = q  - r;

    q = r = 0.0;
    for (i = k; i < D_ll; ++i){
        a = y[i] - sqrt( (double) i + 1.0 );
        q += a*a;
        s += a;
    }

    q += pow(0.5*s, 2) + pow(0.5*s, 4);

    r = x[k]*x[k];
    for (i = k+1; i < D_ul; ++i){
        r += 1e6*x[i]*x[i];
    }

    p2 = q - r;

    f[0] = p1 + p2;
}

void PMM3_leader(int D_ul, int D_ll, int k, double *x, double *y, double *F){
    int i;

    double p1  = 0.0, p2  = 0.0, q = 0.0, r = 0.0;

    for (i = 0; i < k; ++i){
        q += pow(y[i] - (pow(x[i], 3) / 100.0), 2);
        r += x[i]*x[i];
    }


    p1 = q  - r;

    q  = r = 0.0;
    for (i = k; i < D_ll-1; ++i){
        q += 100*pow(y[i]*y[i] - y[i+1], 2 ) + pow(y[i] - 1, 2);
    }
    
    for (i = k; i < D_ul; ++i){
        r += x[i]*x[i];
    }

    p2 = q - r;

    F[0] = fabs(p1) + fabs(p2);
}

void PMM3_follower(int D_ul, int D_ll, int k, double *x, double *y, double *f){
    int i;

    double p1  = 0.0, p2  = 0.0, q = 0.0, r = 0.0;

    for (i = 0; i < k; ++i){
        q += pow(y[i] - (pow(x[i], 3) / 100.0), 2);
        r += x[i]*x[i];
    }


    p1 = q  - r;

    q  = r = 0.0;
    for (i = k; i < D_ll-1; ++i){
        q += 100*pow(y[i]*y[i] - y[i+1], 2 ) + pow(y[i] - 1, 2);
    }
    
    for (i = k; i < D_ul; ++i){
        r += x[i]*x[i];
    }

    p2 = q - r;

    f[0] = p1 + p2;
}

void PMM4_leader(int D_ul, int D_ll, int k, double *x, double *y, double *F){
    int i;

    double p1  = 0.0, p2  = 0.0, q = 0.0, r = 0.0;

    for (i = 0; i < k; ++i){
        q += pow( y[i] - 10.0  / (1.0 + 2.5*x[i]*x[i]),  2);
        r += x[i]*x[i];
    }

    p1 = q  - r;

    q = r = 0.0;
    q = 10.0*(D_ll - k);
    for (i = k; i < D_ll; ++i){
        q += y[i]*y[i]- 10.0*cos( 2*PI *y[i]);
    }

    for (i = k; i < D_ul; ++i){
        r += x[i]*x[i];
    }

    p2 = q - r;

    F[0] = fabs(p1) + fabs(p2);
}

void PMM4_follower(int D_ul, int D_ll, int k, double *x, double *y, double *f){
    int i;

    double p1  = 0.0, p2  = 0.0, q = 0.0, r = 0.0;

    for (i = 0; i < k; ++i){
        q += pow( y[i] - 10.0  / (1.0 + 2.5*x[i]*x[i]),  2);
        r += x[i]*x[i];
    }


    p1 = q  - r;

    q = r = 0.0;
    q = 10.0*(D_ll - k);
    for (i = k; i < D_ll; ++i){
        q += y[i]*y[i]- 10.0*cos( 2*PI *y[i]);
    }

    for (i = k; i < D_ul; ++i){
        r += x[i]*x[i];
    }

    p2 = q - r;

    f[0] = p1 + p2;
}

void PMM5_leader(int D_ul, int D_ll, int k, double *x, double *y, double *F){
    int i;

    double p1  = 0.0, p2  = 0.0, q = 0.0, r = 0.0;

    q = 10.0*k;
    for (i = 0; i < k; ++i){
        q += pow( y[i] - x[i], 2 ) - 10.0*cos(2.0*PI*(y[i] - x[i]));
        r += x[i]*x[i];
    }


    p1 = q  - r;

    q = r = 0.0;
    for (i = k; i < D_ll; ++i){
        q += pow(fabs(y[i]), i - k + 2);
    }
    for (i = k; i < D_ul; ++i){
        r += x[i]*x[i];
    }

    p2 = q - r;

    F[0] = fabs(p1) + fabs(p2);
}

void PMM5_follower(int D_ul, int D_ll, int k, double *x, double *y, double *f){
    int i;

    double p1  = 0.0, p2  = 0.0, q = 0.0, r = 0.0;

    q = 10.0*k;
    for (i = 0; i < k; ++i){
        q += pow( y[i] - x[i], 2 ) - 10.0*cos(2.0*PI*(y[i] - x[i]));
        r += x[i]*x[i];
    }


    p1 = q  - r;

    q = r = 0.0;
    for (i = k; i < D_ll; ++i){
        q += pow(fabs(y[i]), i - k + 2);
    }
    for (i = k; i < D_ul; ++i){
        r += x[i]*x[i];
    }

    p2 = q - r;

    f[0] = p1 + p2;
}

void PMM6_leader(int D_ul, int D_ll, int k, double *x, double *y, double *F){
    int i;

    double p1  = 0.0, p2  = 0.0, q = 0.0, r = 0.0;

    double prod = 1.0;
    for (i = 0; i < k; ++i){
        q += pow( y[i] - x[i], 2 );
        prod *= cos( 10.0*(y[i] - x[i]) / sqrt( i + 1 ) );
        r += x[i]*x[i];
    }

    q = 1.0 + ( 1.0 / 4.0 )*q - prod;

    p1 = q  - r;

    q = 10.0*(D_ll - k);
    r = 10.0*(D_ul - k);
    for (i = k; i < D_ll; ++i){
        q += y[i]*y[i]- 10.0*cos( 2*PI *y[i]);
    }
    for (i = k; i < D_ul; ++i){
        r += x[i]*x[i]- 10.0*cos( 2*PI *x[i]);
    }

    p2 = q - r;

    F[0] = fabs(p1) + fabs(p2);
}

void PMM6_follower(int D_ul, int D_ll, int k, double *x, double *y, double *f){
    int i;

    double p1  = 0.0, p2  = 0.0, q = 0.0, r = 0.0;

    double prod = 1.0;
    for (i = 0; i < k; ++i){
        q += pow( y[i] - x[i], 2 );
        prod *= cos( 10.0*(y[i] - x[i]) / sqrt( i + 1 ) );
        r += x[i]*x[i];
    }

    q = 1.0 + ( 1.0 / 4.0 )*q - prod;

    p1 = q  - r;

    q = 10.0*(D_ll - k);
    r = 10.0*(D_ul - k);
    for (i = k; i < D_ll; ++i){
        q += y[i]*y[i]- 10.0*cos( 2*PI *y[i]);
    }
    for (i = k; i < D_ul; ++i){
        r += x[i]*x[i]- 10.0*cos( 2*PI *x[i]);
    }

    p2 = q - r;

    f[0] = p1 + p2;
}



/*
 * D_ul: upper level dimension
 * D_ll: lower level dimension
 * x: upper level decision vector
 * y: lower level decision vector
 * fnum: function number, i.e., PMM<fnum>
 * Save in F, G the function and constraints values.
*/
void PMM_leader(int D_ul, int D_ll, double *x, double *y, double *F, double *G, int fnum){
    int k = D_ul < D_ll? D_ul : D_ll;
    k /= 2;
    if ( fnum == 1 ) {
        PMM1_leader(D_ul, D_ll, k, x, y, F);
    }else if (fnum == 2) {
        PMM2_leader(D_ul, D_ll, k, x, y, F);
    }else if (fnum == 3) {
        PMM3_leader(D_ul, D_ll, k, x, y, F);
    }else if (fnum == 4) {
        PMM4_leader(D_ul, D_ll, k, x, y, F);
    }else if (fnum == 5) {
        PMM5_leader(D_ul, D_ll, k, x, y, F);
    }else if (fnum == 6) {
        PMM6_leader(D_ul, D_ll, k, x, y, F);
    }/*else if (fnum == 7) {
        PMM7_leader(D_ul, D_ll, k, x, y, F, G);
    }else if (fnum == 8) {
        PMM8_leader(D_ul, D_ll, k, x, y, F, G);
    }else if (fnum == 9) {
        PMM9_leader(D_ul, D_ll, k, x, y, F, G);
    }else if (fnum == 10) {
        PMM10_leader(D_ul, D_ll, k, x, y, F, G);
    }*/
}


/*
 * D_ul: upper level dimension
 * D_ll: lower level dimension
 * x: upper level decision vector
 * y: lower level decision vector
 * fnum: function number, i.e., PMM<fnum>
 * Save in f, g the function and constraints values.
*/
void PMM_follower(int D_ul, int D_ll, double *x, double *y, double *f, double *g, int fnum){
    int k = D_ul < D_ll? D_ul : D_ll;
    k /= 2;
    if ( fnum == 1 ) {
        PMM1_follower(D_ul, D_ll, k, x, y, f);
    }else if (fnum == 2) {
        PMM2_follower(D_ul, D_ll, k, x, y, f);
    }else if (fnum == 3) {
        PMM3_follower(D_ul, D_ll, k, x, y, f);
    }else if (fnum == 4) {
        PMM4_follower(D_ul, D_ll, k, x, y, f);
    }else if (fnum == 5) {
        PMM5_follower(D_ul, D_ll, k, x, y, f);
    }else if (fnum == 6) {
        PMM6_follower(D_ul, D_ll, k, x, y, f);
    }/*else if (fnum == 7) {
        PMM7_follower(D_ul, D_ll, k, x, y, f, g);
    }else if (fnum == 8) {
        PMM8_follower(D_ul, D_ll, k, x, y, f, g);
    }else if (fnum == 9) {
        PMM9_follower(D_ul, D_ll, k, x, y, f, g);
    }else if (fnum == 10) {
        PMM10_follower(D_ul, D_ll, k, x, y, f, g);
    }*/
}


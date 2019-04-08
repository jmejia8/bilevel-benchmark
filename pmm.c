#include "constants.h"

void PMM_config(int *settings, int fnum){
    if ( 6 <= fnum && fnum <= 8) {
        settings[0] = 1; settings[1] = 1;
    }else if (fnum == 9 || fnum == 10) {
        settings[0] = 2; settings[1] = 2;
    }else{
        settings[0] = 0; settings[1] = 0;
    }
}

void PMM_Psi(int m, int n, double *x, double *y, int fnum){
    int i = 1;

    if (fnum == 2 || fnum == 7){
       for (i = 0; i < n; ++i) { y[i] = i < m ? pow(x[i], 3) : 0.0; }
    }else if (1 <= fnum && fnum <= 10) {
       for (i = 0; i < n; ++i) { y[i] = i < m ? x[i] : 0.0; }
    }
}

void PMM1_leader(int m, int n, double *x, double *y, double *F){
    int i;

    double P = 0.0, q  = 0.0, Q  = 0.0;

    double Q2 = 0.0;
    for (i = 0; i < m; ++i) {
        q += fabs(x[i]);
        Q += pow(y[i]-x[i],2);
        Q2 += y[i]-x[i];
    }

    P += q;
    for (i = m; i < n; ++i){P += fabs(x[i]);}

    Q += pow(0.5*Q2,2) + pow(0.5*Q2,4);

    F[0] = P - q*Q;
}

void PMM1_follower(int m, int n, double *x, double *y, double *f){
    int i;

    double p = 0.0, q  = 0.0, Q  = 0.0;

    double Q2 = 0.0;
    for (i = 0; i < m; ++i) {
        q += fabs(x[i]);
        Q += pow(y[i]-x[i],2);
        Q2 += y[i]-x[i];
    }

    for (i = m; i < n; ++i){p += pow(y[i], 2);}

    Q += pow(0.5*Q2,2) + pow(0.5*Q2,4);

    f[0] = p + q*Q;
}

void PMM2_leader(int m, int n, double *x, double *y, double *F){
    int i;

    double P = 0.0, q  = 0.0, Q  = 0.0;

    for (i = 0; i < m; ++i) {
        q += pow(x[i], 2);
        Q += pow(pow(x[i], 3) - y[i],2);
    }

    for (i = m; i < n; ++i){P += pow(x[i], 2);}

    P = q + 1e6*P;

    F[0] = P - q*Q;
}

void PMM2_follower(int m, int n, double *x, double *y, double *f){
    int i;

    double p = 0.0, q  = 0.0, Q  = 0.0;

    for (i = 0; i < m; ++i) {
        q += pow(x[i], 2);
        Q += pow(pow(x[i], 3) - y[i],2);
    }

    for (i = m; i < n; ++i){p += pow(y[i], 2);}

    f[0] = p + q*Q;
}

void PMM3_leader(int m, int n, double *x, double *y, double *F){
    int i;

    double P = 0.0, q  = 0.0, Q  = 0.0;

    Q = pow(x[0]-y[0], 2);
    q = pow(x[0], 2);
    for (i = 1; i < m; ++i) {
        q += pow(x[i], 2);
        Q += 1e6*pow(x[i] - y[i], 2);
    }

    P = q;
    for (i = m; i < n; ++i){
        P += fabs(x[i]);
    }


    F[0] = P - q*Q;
}

void PMM3_follower(int m, int n, double *x, double *y, double *f){
    int i;

    double p = 0.0, q  = 0.0, Q  = 0.0;

    Q = pow(x[0]-y[0], 2);
    q = pow(x[0], 2);
    for (i = 1; i < m; ++i) {
        q += pow(x[i], 2);
        Q += 1e6*pow(x[i] - y[i], 2);
    }

    for (i = m; i < n; ++i){p += pow(y[i], 2);}

    f[0] = p + q*Q;
}

void PMM4_leader(int m, int n, double *x, double *y, double *F){
    int i;

    double P = 0.0, q  = 0.0, Q  = 0.0;

    Q = (double) 10*m;
    for (i = 0; i < m; ++i) {
        q += pow(x[i], 2);
        Q += pow(x[i] - y[i], 2) + 10*cos(PI*fabs(x[i] - y[i]) / ( 0.001 + pow(x[m], 2) ));
    }

    for (i = m; i < n; ++i){
        P += pow(x[i], 2);
    }

    P += 1e6*q ;

    F[0] = P - q*Q;
}

void PMM4_follower(int m, int n, double *x, double *y, double *f){
    int i;

    double p = 0.0, q  = 0.0, Q  = 0.0;

    Q = (double) 10*m;
    for (i = 0; i < m; ++i) {
        q += pow(x[i], 2);
        Q += pow(x[i] - y[i], 2) + 10*cos(PI*fabs(x[i] - y[i]) / ( 0.001 + pow(x[m], 2) ));
    }

    for (i = m; i < n; ++i){
        p += pow(y[i], 2);
    }

    f[0] = p + q*Q;
}

void PMM5_leader(int m, int n, double *x, double *y, double *F){
    int i;

    double P = 0.0, q  = 0.0, Q  = 0.0;

    q += pow(x[m], 2);
    for (i = 0; i < m-1; ++i) {
        q += pow(x[i], 2);
        double o = pow(x[i]-y[i], 2) + pow(x[i+1]-y[i+1], 2);
        Q += (pow(sin( sqrt(o) ), 2) - 0.5) / pow(1.0 + 0.001*(o) , 2);
    }

    P = 10*q + (double) 10.0*(n-m);
    for (i = m; i < n; ++i){
        P += pow(x[i], 2) + 10*cos(2*PI*x[i]);
    }

    F[0] = P - q*Q;
}

void PMM5_follower(int m, int n, double *x, double *y, double *f){
    int i;

    double p = 0.0, q  = 0.0, Q  = 0.0;

    q += pow(x[m], 2);
    for (i = 0; i < m-1; ++i) {
        q += pow(x[i], 2);
        double o = pow(x[i]-y[i], 2) + pow(x[i+1]-y[i+1], 2);
        Q += (pow(sin( sqrt(o) ), 2) - 0.5) / pow(1.0 + 0.001*(o) , 2);
    }

    for (i = m; i < n; ++i){
        p += pow(y[i], 2);
    }

    f[0] = p + q*Q;
}

void PMM6_leader(int m, int n, double *x, double *y, double *F, double *G){
    int i;

    double P = 0.0, q  = 0.0, Q  = 0.0;

    double Q2 = 0.0;
    for (i = 0; i < m; ++i) {
        q += fabs(x[i]);
        Q += pow(y[i]-x[i],2);
        Q2 += y[i]-x[i];
    }

    for (i = m; i < n; ++i) {
        P =(P < fabs(x[i])) ? fabs(x[i]) : P ;
    }

    Q += pow(0.5*Q2,2) + pow(0.5*Q2,4);

    F[0] = P - q*Q;
    ////////////////////////////////////////////////////////////////////////////
    ///// Constraints
    ////////////////////////////////////////////////////////////////////////////
    G[0] = 0.0;
    double min = fabs(y[0]);
    for (int i = 0; i < n; ++i) {
        G[0] += pow(x[i], 2) + 50*cos(2*PI*x[i]);
        min = (min > fabs(y[i])) ? fabs(y[i]) : min ;
    }
    G[0] -= min; 
    ////////////////////////////////////////////////////////////////////////////
}

void PMM6_follower(int m, int n, double *x, double *y, double *f, double *g){
    int i;

    double p = 0.0, q  = 0.0, Q  = 0.0;

    double Q2 = 0.0;
    for (i = 0; i < m; ++i) {
        q += fabs(x[i]);
        Q += pow(y[i]-x[i],2);
        Q2 += y[i]-x[i];
    }

    for (i = m; i < n; ++i){p += pow(y[i], 2);}

    Q += pow(0.5*Q2,2) + pow(0.5*Q2,4);

    f[0] = p + q*Q;
    ////////////////////////////////////////////////////////////////////////////
    ///// Constraints
    ////////////////////////////////////////////////////////////////////////////
    g[0] = 1.0;
    for (int i = 0; i < n; ++i) {
        g[0] *= y[i] - x[i];
    }
    ////////////////////////////////////////////////////////////////////////////
}

void PMM7_leader(int m, int n, double *x, double *y, double *F, double *G){
    int i;

    double P = 0.0, q  = 0.0, Q  = 0.0;

    for (i = 0; i < m; ++i) {
        q += pow(x[i], 2);
        Q += pow(pow(x[i], 3) - y[i],2);
        P += fabs(x[i]);
    }

    for (i = m; i < n; ++i){P += fabs(x[i]);}

    F[0] = P - q*Q;

    ////////////////////////////////////////////////////////////////////////////
    ///// Constraints
    ////////////////////////////////////////////////////////////////////////////
    G[0] = 1.0;
    q = 0.0;
    for (int i = 0; i < n; ++i) {
        G[0] *= x[i];
        q += fabs(y[i]);
    }
    G[0] -= q; 
    ////////////////////////////////////////////////////////////////////////////
}

void PMM7_follower(int m, int n, double *x, double *y, double *f, double *g){
    int i;

    double p = 0.0, q  = 0.0, Q  = 0.0;

    for (i = 0; i < m; ++i) {
        q += pow(x[i], 2);
        Q += pow(pow(x[i], 3) - y[i],2);
    }

    for (i = m; i < n; ++i){p += pow(y[i], 2);}

    f[0] = p + q*Q;

    ////////////////////////////////////////////////////////////////////////////
    ///// Constraints
    ////////////////////////////////////////////////////////////////////////////
    g[0] = 0.0;
    double min = fabs(x[0]);
    for (int i = 0; i < n; ++i) {
        g[0] += pow(y[i], 2) + 50*cos(2*PI*y[i]);
        min = min > fabs(x[i]) ? fabs(x[i]) : min;
    }
    g[0] -= min;
    ////////////////////////////////////////////////////////////////////////////
}

void PMM8_leader(int m, int n, double *x, double *y, double *F, double *G){
    int i;

    double P = 0.0, q  = 0.0, Q  = 0.0;

    Q = pow(x[0]-y[0], 2);
    P += fabs(x[0]);
    q = fabs(x[0]);
    for (i = 1; i < m; ++i) {
        P += fabs(x[i]);
        Q += 1e6*pow(x[i] - y[i], 2);
        q = (q > fabs(x[i])) ? fabs(x[i]) : q;
    }

    for (i = m; i < n; ++i){
        P += 100*pow(x[i], 2);
    }


    F[0] = P - q*Q;

    ////////////////////////////////////////////////////////////////////////////
    ///// Constraints
    ////////////////////////////////////////////////////////////////////////////
    G[0] = 0.0;

    for (int i = 0; i < n; ++i) {
        G[0] += pow(y[i], 2) + 100*cos(2*PI*x[i]);
    }
    ////////////////////////////////////////////////////////////////////////////
}

void PMM8_follower(int m, int n, double *x, double *y, double *f, double *g){
    int i;

    double p = 0.0, q  = 0.0, Q  = 0.0;

    Q = pow(x[0]-y[0], 2);
    q = fabs(x[0]);
    for (i = 1; i < m; ++i) {
        Q += 1e6*pow(x[i] - y[i], 2);
        q = (q > fabs(x[i])) ? fabs(x[i]) : q;
    }

    p = (double) 10.0*(n-m);
    for (i = m; i < n; ++i){
        p += pow(y[i], 2) - 10*cos(2*PI*y[i]);
    }

    f[0] = p + q*Q;

    ////////////////////////////////////////////////////////////////////////////
    ///// Constraints
    ////////////////////////////////////////////////////////////////////////////
    g[0] = 1.0;
    for (int i = 0; i < n; ++i) {
        g[0] *= sin(x[i] + 2*y[i]);
    }
    ////////////////////////////////////////////////////////////////////////////
}

void PMM9_leader(int m, int n, double *x, double *y, double *F, double *G){
    int i;

    double P = 0.0, q  = 0.0, Q  = 0.0;

    P = pow(x[m], 2);
    q += fabs(x[m]);
    for (i = 0; i < m-1; ++i) {
        P += pow(x[i], 2);
        q += fabs(x[m]);
        double o = pow(x[i]-y[i], 2) + pow(x[i+1]-y[i+1], 2);
        Q += (pow(sin( sqrt(o) ), 2) - 0.5) / pow(1.0 + 0.001*(o) , 2);
    }

    q *= 100;

    P *= 1e6;
    for (i = m; i < n; ++i){ P += pow(x[i], 2);}

    F[0] = P - q*Q;
    ////////////////////////////////////////////////////////////////////////////
    ///// Constraints
    ////////////////////////////////////////////////////////////////////////////
    G[0] = 0.0;
    G[1] = 0.0;
    double min = fabs(y[0]);
    for (int i = 0; i < n; ++i) {
        G[0] += pow(x[i], 2) + 100*cos(2*PI*x[i]);
        G[1] += (i>=m) ? x[i] : 0;
        min = (min > fabs(y[i])) ? fabs(y[i]) : min ;
    }

    G[0] -= min;
    ////////////////////////////////////////////////////////////////////////////
}

void PMM9_follower(int m, int n, double *x, double *y, double *f, double *g){
    int i,j;

    double p = 0.0, q  = 0.0, Q  = 0.0;

    q += fabs(x[m]);
    for (i = 0; i < m; ++i) {
        double o = pow(x[i]-y[i], 2) + pow(x[i+1]-y[i+1], 2);
        Q += (i < m-1) ? (pow(sin( sqrt(o) ), 2) - 0.5) / pow(1.0 + 0.001*(o) , 2) : 0;
        
        for (j = m; j < i; ++j) {
            p += pow(y[j], 2);        
        }
        q += fabs(x[m]);
    }

    q *= 100;


    f[0] = p + q*Q;
    ////////////////////////////////////////////////////////////////////////////
    ///// Constraints
    ////////////////////////////////////////////////////////////////////////////
    g[0] = 1.0;
    g[1] = 0.0;
    double max = pow(y[0], 2);
    for (int i = 0; i < n; ++i) {
        g[0] *= y[i] - x[i];
        g[1] += (i > m) ? pow(y[i], 2) : 0 ;
        max = (i > m && max < pow(y[i], 2)) ? pow(y[i], 2) : max;
    }
    g[1] -= (double) (n-m)  * max;
    ////////////////////////////////////////////////////////////////////////////
}

void PMM10_leader(int m, int n, double *x, double *y, double *F, double *G){
    int i;

    double P = 0.0, q  = 0.0, Q  = 0.0;

    P = (double) 10.0*n + pow(x[m], 2) + 100*cos(2*PI*x[m]);;
    q += pow(x[m], 2) + floor(fabs(x[m]));
    for (i = 0; i < m-1; ++i) {
        P += pow(x[i], 2) + 10*cos(2*PI*x[i]);;
        q += pow(x[m], 2) + floor(fabs(x[m]));
        double o = pow(x[i]-y[i], 2) + pow(x[i+1]-y[i+1], 2);
        Q += (pow(sin( sqrt(o) ), 2) - 0.5) / pow(1.0 + 0.001*(o) , 2);
    }

    P = 10*q + (double) 10.0*(n-m);
    for (i = m; i < n; ++i){
        P += pow(x[i], 2) + 10*cos(2*PI*x[i]);
    }

    F[0] = P - q*Q;

    ////////////////////////////////////////////////////////////////////////////
    ///// Constraints
    ////////////////////////////////////////////////////////////////////////////
    G[0] = 0.0;
    G[1] = 1.0;
    double min = fabs(y[0]);
    for (int i = 0; i < n; ++i) {
        G[0] += pow(x[i], 2) + 50*cos(2*PI*x[i]);
        G[1] *= (i>=m) ? x[i] : 1.0;
        min = (min > fabs(y[i])) ? fabs(y[i]) : min ;
    }

    G[0] -= min;
    ////////////////////////////////////////////////////////////////////////////
}

void PMM10_follower(int m, int n, double *x, double *y, double *f, double *g){
    int i;

    double p = 0.0, q  = 0.0, Q  = 0.0;

    q += pow(x[m], 2) + floor(fabs(x[m]));
    for (i = 0; i < m-1; ++i) {
        q += pow(x[m], 2) + floor(fabs(x[m]));
        double o = pow(x[i]-y[i], 2) + pow(x[i+1]-y[i+1], 2);
        Q += (pow(sin( sqrt(o) ), 2) - 0.5) / pow(1.0 + 0.001*(o) , 2);
    }

    for (i = m; i < n; ++i){
        p = pow(y[i], 2);
    }

    f[0] = p - q*Q;

    ////////////////////////////////////////////////////////////////////////////
    ///// Constraints
    ////////////////////////////////////////////////////////////////////////////
    g[0] = 0.0;
    g[1] = 1.0;
    double min = fabs(y[0]);
    for (int i = 0; i < n; ++i) {
        g[0] += pow(y[i], 2) + 50*cos(2*PI*y[i]);
        g[1] *= (i>=m) ? y[i] : 1.0;
        min = (min > fabs(x[i])) ? fabs(x[i]) : min ;
    }

    g[0] -= min;
    ////////////////////////////////////////////////////////////////////////////
}

void PMM_leader(int D_ul, int D_ll, double *x, double *y, double *F, double *G, int fnum){
    if ( fnum == 1 ) {
        PMM1_leader(D_ul, D_ll, x, y, F);
    }else if (fnum == 2) {
        PMM2_leader(D_ul, D_ll, x, y, F);
    }else if (fnum == 3) {
        PMM3_leader(D_ul, D_ll, x, y, F);
    }else if (fnum == 4) {
        PMM4_leader(D_ul, D_ll, x, y, F);
    }else if (fnum == 5) {
        PMM5_leader(D_ul, D_ll, x, y, F);
    }else if (fnum == 6) {
        PMM6_leader(D_ul, D_ll, x, y, F, G);
    }else if (fnum == 7) {
        PMM7_leader(D_ul, D_ll, x, y, F, G);
    }else if (fnum == 8) {
        PMM8_leader(D_ul, D_ll, x, y, F, G);
    }else if (fnum == 9) {
        PMM9_leader(D_ul, D_ll, x, y, F, G);
    }else if (fnum == 10) {
        PMM10_leader(D_ul, D_ll, x, y, F, G);
    }
}

void PMM_follower(int D_ul, int D_ll, double *x, double *y, double *f, double *g, int fnum){
    if ( fnum == 1 ) {
        PMM1_follower(D_ul, D_ll, x, y, f);
    }else if (fnum == 2) {
        PMM2_follower(D_ul, D_ll, x, y, f);
    }else if (fnum == 3) {
        PMM3_follower(D_ul, D_ll, x, y, f);
    }else if (fnum == 4) {
        PMM4_follower(D_ul, D_ll, x, y, f);
    }else if (fnum == 5) {
        PMM5_follower(D_ul, D_ll, x, y, f);
    }else if (fnum == 6) {
        PMM6_follower(D_ul, D_ll, x, y, f, g);
    }else if (fnum == 7) {
        PMM7_follower(D_ul, D_ll, x, y, f, g);
    }else if (fnum == 8) {
        PMM8_follower(D_ul, D_ll, x, y, f, g);
    }else if (fnum == 9) {
        PMM9_follower(D_ul, D_ll, x, y, f, g);
    }else if (fnum == 10) {
        PMM10_follower(D_ul, D_ll, x, y, f, g);
    }
}

double PMM_test(int m, int n){
    int FNUN = 10, fnum;

    double *x = array(n);
    double *y = array(n);
    double F[1], f[1], f_sum = 0;
    double G[2], g[2];
    int settings[2];

    PMM_config(settings, fnum);
    randm(0, 1, x, 1*n);
    randm(0, 1, y, 1*n);

    for (fnum = 0; fnum < FNUN; ++fnum){
        PMM_leader(n, n, x, y, F, G, fnum);
        PMM_follower(n, n, x, y, f, g, fnum);

        printf("PMM%d \t F = %lf \t f = %lf\n", fnum, F[0], f[0]);
    } 
    return 1;
}

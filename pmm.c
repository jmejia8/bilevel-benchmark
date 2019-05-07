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

void PMM1_leader(int D_ul, int D_ll, int m, double *x, double *y, double *F){
    int i;

    double P = 0.0, q  = 0.0, Q  = 0.0, phi=0.0;

    double Q2 = 0.0;
    for (i = 0; i < m; ++i) {
        q += fabs(x[i]);
        phi = 0.01*pow(x[i], 3);
        Q += pow(y[i]-phi,2);
        Q2 += y[i]-phi;
    }

    P += q;
    for (i = m; i < D_ul; ++i){P += fabs(x[i]);}

    Q += pow(0.5*Q2,2) + pow(0.5*Q2,4);

    F[0] = P - q*Q;
}

void PMM1_follower(int D_ul, int D_ll, int m, double *x, double *y, double *f){
    int i;

    double p = 0.0, q  = 0.0, Q  = 0.0, phi;

    double Q2 = 0.0;
    for (i = 0; i < m; ++i) {
        q += fabs(x[i]);
        phi = 0.01*pow(x[i], 3);
        Q += pow(y[i]-phi,2);
        Q2 += y[i]-phi;
    }

    for (i = m; i < D_ll; ++i){p += pow(y[i], 2);}

    Q += pow(0.5*Q2,2) + pow(0.5*Q2,4);

    f[0] = p + q*Q;
}

void PMM2_leader(int D_ul, int D_ll, int m, double *x, double *y, double *F){
    int i;

    double P = 0.0, q  = 0.0, Q  = 0.0, phi;

    for (i = 0; i < m; ++i) {
        q += pow(x[i], 2);
        phi = x[i] * sin(x[i]);
        Q += pow(y[i] - phi,2);
    }

    for (i = m; i < D_ul; ++i){P += pow(x[i], 2);}

    P = 1e6*q + P;

    F[0] = P - q*Q;
}

void PMM2_follower(int D_ul, int D_ll, int m, double *x, double *y, double *f){
    int i;

    double p = 0.0, q  = 0.0, Q  = 0.0, phi;

    for (i = 0; i < m; ++i) {
        q += pow(x[i], 2);
        phi = x[i] * sin(x[i]);
        Q += pow(y[i] - phi, 2);
    }

    for (i = m; i < D_ll; ++i){p += pow(y[i], 2);}

    f[0] = p + q*Q;
}

void PMM3_leader(int D_ul, int D_ll, int m, double *x, double *y, double *F){
    int i;

    double P = 0.0, q  = 0.0, Q  = 0.0, phi;

    phi = (10.0/ (1.0 + 2.5*pow(x[0], 2))) - 10.0;
    Q = pow(y[0]-phi, 2);
    q = pow(x[0], 2);
    for (i = 1; i < m; ++i) {
        q += pow(x[i], 2);
        phi = (10.0/ (1.0 + 2.5*pow(x[i], 2))) - 10.0;
        Q += 1e6*pow(y[i] - phi, 2);
    }

    P = q;
    for (i = m; i < D_ul; ++i){
        P += fabs(x[i]);
    }


    F[0] = P - q*Q;
}

void PMM3_follower(int D_ul, int D_ll, int m, double *x, double *y, double *f){
    int i;

    double p = 0.0, q  = 0.0, Q  = 0.0, phi;

    phi = (10.0/ (1.0 + 2.5*pow(x[0], 2))) - 10.0;
    Q = pow(y[0]-phi, 2);
    q = pow(x[0], 2);
    for (i = 1; i < m; ++i) {
        q += pow(x[i], 2);
        phi = (10.0/ (1.0 + 2.5*pow(x[i], 2))) - 10.0;
        Q += 1e6*pow(y[i] - phi, 2);
    }

    for (i = m; i < D_ll; ++i){p += pow(y[i], 2);}

    f[0] = p + q*Q;
}

void PMM4_leader(int D_ul, int D_ll, int m, double *x, double *y, double *F){
    int i;

    double P = 0.0, q  = 0.0, Q  = 0.0, phi;

    Q = (double) 10*m;
    for (i = 0; i < m; ++i) {
        q += pow(x[i], 2);
        phi = 0.01*pow(x[i], 3) + sin(2.0*PI*x[i]);
        Q += pow(y[i] - phi, 2) - 10.0*cos(PI*fabs(y[i] - phi) / ( 0.001 + pow(x[m], 2) ));
    }

    Q *= pow(x[0], 4);


    for (i = m; i < D_ul; ++i){
        P += pow(x[i], 2);
    }

    P += 1e6*q ;

    F[0] = P - q*Q;
}

void PMM4_follower(int D_ul, int D_ll, int m, double *x, double *y, double *f){
    int i;

    double p = 0.0, q  = 0.0, Q  = 0.0, phi;

    Q = (double) 10*m;
    for (i = 0; i < m; ++i) {
        q += pow(x[i], 2);
        phi = 0.01*pow(x[i], 3) + sin(2.0*PI*x[i]);
        Q += pow(y[i] - phi, 2) - 10.0*cos(PI*fabs(y[i] - phi) / ( 0.001 + pow(x[m], 2) ));
    }

    for (i = m; i < D_ll; ++i){
        p += pow(y[i], 4);
    }

    f[0] = p + q*Q;
}

void PMM5_leader(int D_ul, int D_ll, int m, double *x, double *y, double *F){
    int i;

    double P = 0.0, q  = 0.0, Q  = 0.0;
    double prod = 1.0;
    for (i = 0; i < m; ++i) {
        q += pow(fabs(x[i]), i+2);
        Q += pow(y[i] - x[i], 2);
        prod *= cos(10.0*fabs(y[i] - x[i]) / sqrt((double) i + 1));
    }

    Q = 1.0 + 0.025*Q - prod;

    P = 10.0*q + (double) 10.0*(D_ul-m);
    for (i = m; i < D_ul; ++i){
        P += pow(x[i], 2) - 10.0*cos(2.0*PI*x[i]);
    }

    F[0] = P - q*Q;
}

void PMM5_follower(int D_ul, int D_ll, int m, double *x, double *y, double *f){
    int i;

    double p = 0.0, q  = 0.0, Q  = 0.0, prod=1.0;

    q += pow(x[m-1], 2);
    for (i = 0; i < m-1; ++i) {
        q += pow(fabs(x[i]), i+2);
        Q += pow(y[i] - x[i], 2);
        prod *= cos(10.0*fabs(y[i] - x[i]) / sqrt((double) i + 1));
    }

    Q = 1.0 + 0.025*Q - prod;


    for (i = m; i < D_ll; ++i){
        p += pow(y[i], 2);
    }

    f[0] = p + q*Q;
}

void PMM6_leader(int D_ul, int D_ll, int m, double *x, double *y, double *F, double *G){
    int i;

    double P = 0.0, q  = 0.0, Q  = 0.0, phi;

    double Q2 = 0.0;
    for (i = 0; i < m; ++i) {
        q += fabs(x[i]);
        phi = 0.01*pow(x[i], 3);
        Q += pow(y[i]-phi,2);
        Q2 += y[i]-phi;
    }

    P = fabs(x[m]);
    for (i = m+1; i < D_ul; ++i) {
        P =(P < fabs(x[i])) ? fabs(x[i]) : P ;
    }

    for (i = 0; i < m; ++i) {
        P += fabs(x[i]);
    }

    Q += pow(0.5*Q2,2) + pow(0.5*Q2,4);

    F[0] = P - q*Q;
    ////////////////////////////////////////////////////////////////////////////
    ///// Constraints
    ////////////////////////////////////////////////////////////////////////////
    G[0] = 0.0;
    double min = fabs(y[0]);
    for (i = 1; i < D_ll; ++i) { min = (min > fabs(y[i])) ? fabs(y[i]) : min ; }
    for (i = 0; i < D_ul; ++i) {
        G[0] += pow(x[i], 2) - 50*sin(2*PI*x[i]);
    }
    G[0] -= min; 
    ////////////////////////////////////////////////////////////////////////////
}

void PMM6_follower(int D_ul, int D_ll, int m, double *x, double *y, double *f, double *g){
    int i;
    g[0] = 1.0;

    double p = 0.0, q  = 0.0, Q  = 0.0, phi;

    double Q2 = 0.0;
    for (i = 0; i < m; ++i) {
        q += fabs(x[i]);
        phi = 0.01*pow(x[i], 3);
        Q += pow(y[i]-phi,2);
        Q2 += y[i]-phi;
        // Constraints
        g[0] *= y[i] - phi;
    }

    for (i = m; i < D_ll; ++i){p += pow(y[i], 2);}

    Q += pow(0.5*Q2,2) + pow(0.5*Q2,4);

    f[0] = p + q*Q;
    ////////////////////////////////////////////////////////////////////////////
    ///// Constraints
    ////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////
}

void PMM7_leader(int D_ul, int D_ll, int m, double *x, double *y, double *F, double *G){
    int i;

    double P = 0.0, q  = 0.0, Q  = 0.0, phi;

    for (i = 0; i < m; ++i) {
        q += pow(x[i], 2);
        phi = 10.0*exp( -0.01*pow(x[i],2) )*sin(x[i]);
        Q += pow(y[i] - phi,2);
        P += fabs(x[i]);
    }

    for (i = m; i < D_ul; ++i){P += fabs(x[i]);}

    F[0] = P - q*Q;

    ////////////////////////////////////////////////////////////////////////////
    ///// Constraints
    ////////////////////////////////////////////////////////////////////////////
    G[0] = 1.0;
    q = 0.0;
    for (i = 0; i < D_ul; ++i) {
        G[0] *= x[i];
    }
    for (i = 0; i < D_ll; ++i) {
        q += fabs(y[i]);
    }
    G[0] -= q; 
    ////////////////////////////////////////////////////////////////////////////
}

void PMM7_follower(int D_ul, int D_ll, int m, double *x, double *y, double *f, double *g){
    int i;

    double p = 0.0, q  = 0.0, Q  = 0.0, phi;

    for (i = 0; i < m; ++i) {
        q += pow(x[i], 2);
        phi = 10.0*exp( -0.01*pow(x[i],2) )*sin(x[i]);
        Q += pow(y[i] - phi,2);
    }

    for (i = m; i < D_ll; ++i){p += pow(y[i], 2);}

    f[0] = p + q*Q;

    ////////////////////////////////////////////////////////////////////////////
    ///// Constraints
    ////////////////////////////////////////////////////////////////////////////
    g[0] = 0.0;
    double min = fabs(x[0]);
    for (i = 0; i < D_ll; ++i) {
        g[0] += pow(y[i], 2) + 50*sin(2*PI*y[i]);
        min = min > fabs(x[i]) ? fabs(x[i]) : min;
    }
    g[0] -= min;
    ////////////////////////////////////////////////////////////////////////////
}

void PMM8_leader(int D_ul, int D_ll, int m, double *x, double *y, double *F, double *G){
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

    for (i = m; i < D_ul; ++i){
        P += 100*pow(x[i], 2);
    }


    F[0] = P - q*Q;

    ////////////////////////////////////////////////////////////////////////////
    ///// Constraints
    ////////////////////////////////////////////////////////////////////////////
    G[0] = 0.0;

    for (i = 0; i < D_ll; ++i) {
        G[0] += pow(y[i], 2);
    }
    for (i = 0; i < D_ul; ++i) {
        G[0] += 100*sin(2*PI*x[i]);

    }
    ////////////////////////////////////////////////////////////////////////////
}

void PMM8_follower(int D_ul, int D_ll, int m, double *x, double *y, double *f, double *g){
    int i;

    double p = 0.0, q  = 0.0, Q  = 0.0;

    Q = pow(x[0]-y[0], 2);
    q = fabs(x[0]);
    for (i = 1; i < m; ++i) {
        Q += 1e6*pow(x[i] - y[i], 2);
        q = (q > fabs(x[i])) ? fabs(x[i]) : q;
    }

    p = (double) 10.0*(D_ll-m);
    for (i = m; i < D_ll; ++i){
        p += pow(y[i], 2) - 10.0*cos(2*PI*y[i]);
    }

    f[0] = p + q*Q;

    ////////////////////////////////////////////////////////////////////////////
    ///// Constraints
    ////////////////////////////////////////////////////////////////////////////
    g[0] = 1.0;
    for (i = 0; i < D_ll; ++i) {
        g[0] *= sin(y[i]);
    }

    q = 0.0;
    for (i = 0; i < D_ul; ++i){
        q += x[i];
    }

    g[0] -= 0.3*sin(q);
    ////////////////////////////////////////////////////////////////////////////
}

void PMM9_leader(int D_ul, int D_ll, int m, double *x, double *y, double *F, double *G){
    int i;

    double P = 0.0, q  = 0.0, Q  = 0.0, prod = 1.0;

    for (i = 0; i < m; ++i) {
        P += pow(x[i], 2);
        q += fabs(x[i]);
        Q += pow(y[i] - x[i], 2);
        prod *= cos(10.0*fabs(y[i] - x[i]) / sqrt((double) i + 1));
    }

    Q = 1.0 + 0.025*Q - prod;

    q *= 100.0;

    P *= 1e6;
    for (i = m; i < D_ul; ++i){ P += pow(x[i], 2);}

    F[0] = P - q*Q;
    ////////////////////////////////////////////////////////////////////////////
    ///// Constraints
    ////////////////////////////////////////////////////////////////////////////
    G[0] = 0.0;
    G[1] = 0.0;
    double min = fabs(y[0]);
    for (i = 0; i < D_ul; ++i) {
        G[0] += pow(x[i], 2) + 50*sin(2*PI*x[i]);
        G[1] += (i>=m) ? x[i] : 0;
    }
    for (i = i; i < D_ll; ++i) {
        min = (min > fabs(y[i])) ? fabs(y[i]) : min ;
    }

    G[0] -= min;
    ////////////////////////////////////////////////////////////////////////////
}

void PMM9_follower(int D_ul, int D_ll, int m, double *x, double *y, double *f, double *g){
    int i,j;

    double p = 0.0, q  = 0.0, Q  = 0.0, prod = 1.0;

    for (i = 0; i < m; ++i) {
        q += fabs(x[i]);
        Q += pow(y[i] - x[i], 2);
        prod *= cos(10.0*fabs(y[i] - x[i]) / sqrt((double) i + 1));
    }

    Q = 1.0 + 0.025*Q - prod;

    for (i = m+1; i < D_ll; ++i) {
        for (j = m+1; j < i; ++j) {
            p += pow(y[j], 2);        
        }
    }
    

    q *= 100;


    f[0] = p + q*Q;
    ////////////////////////////////////////////////////////////////////////////
    ///// Constraints
    ////////////////////////////////////////////////////////////////////////////
    g[0] = 1.0;
    g[1] = 0.0;
    double max = pow(y[0], 2);
    for (i = 0; i < m; ++i) {
        g[0] *= y[i] - x[i];
    }
    for (i = m; i < D_ul; ++i) {
        g[1] += pow(y[i], 2) ;
        max = (max < pow(y[i], 2)) ? pow(y[i], 2) : max;
    }
    g[1] -= (double) (D_ll -m)  * max;
    ////////////////////////////////////////////////////////////////////////////
}

void PMM10_leader(int D_ul, int D_ll, int m, double *x, double *y, double *F, double *G){
    int i;

    double P = 0.0, q  = 0.0, Q = 10.0*m, phi;


    for (i = 0; i < m; ++i) {
        P += pow(x[i], 2) - 10.0*cos(2*PI*x[i]);;
        q += 0.1*pow(x[i], 2) + floor(fabs(x[i]));
        phi = -x[i]*fabs( sin(PI*x[i]) );
        Q += pow(y[i] - phi, 2) - 10.0*cos(PI*fabs(y[i] - phi) / ( 0.001 + pow(x[m], 2) ));
    }

    for (i = m; i < D_ul ; ++i){
        P += pow(x[i], 2) - 10.0*cos(2*PI*x[i]);
    }

    P += 10.0*D_ul;
    
    F[0] = P - q*Q;

    ////////////////////////////////////////////////////////////////////////////
    ///// Constraints
    ////////////////////////////////////////////////////////////////////////////
    G[0] = 0.0;
    G[1] = 1.0;
    for (i = 0; i < D_ul; ++i) {
        G[0] += pow(x[i], 2) + 50*sin(2*PI*x[i]);
        G[1] *= (i>=m) ? x[i] : 1.0;
    }

    double min = fabs(y[0]);
    for (i = 0; i < D_ll; ++i) {
        min = (min > fabs(y[i])) ? fabs(y[i]) : min ;
    }

    G[0] -= min;
    ////////////////////////////////////////////////////////////////////////////
}

void PMM10_follower(int D_ul, int D_ll, int m, double *x, double *y, double *f, double *g){
    int i;

    double p = 0.0, q  = 0.0, Q = 10.0*m, phi;

    q += pow(x[m], 2) + floor(fabs(x[m]));
    for (i = 0; i < m; ++i) {
        q += pow(x[m], 2) + floor(fabs(x[m]));
        phi = -x[i]*fabs( sin(PI*x[i]) );
        Q += pow(y[i] - phi, 2) - 10.0*cos(PI*fabs(y[i] - phi) / ( 0.001 + pow(x[m], 2) ));
        
    }

    for (i = m; i < D_ll; ++i){
        p = pow(y[i], 2);
    }

    f[0] = p - q*Q;

    ////////////////////////////////////////////////////////////////////////////
    ///// Constraints
    ////////////////////////////////////////////////////////////////////////////
    g[0] = 0.0;
    g[1] = 1.0;
    double min = fabs(y[0]);
    for (i = 0; i < m; ++i) {
        g[0] += pow(y[i], 2) + 10*sin(2*PI*y[i]);
        min = (min > fabs(x[i])) ? fabs(x[i]) : min ;
    }

    for (i = m; i < D_ll; ++i) {
        g[1] *= y[i];
    }

    g[0] -= min;
    ////////////////////////////////////////////////////////////////////////////
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


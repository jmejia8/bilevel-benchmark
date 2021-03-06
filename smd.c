#include "constants.h"

// F (x_u , x_l ) = F_1 (x_u1 ) + F_2 (x_l1 ) + F_3 (x_u2 , x_l2 )
// f (x_u , x_l ) = f_1 (x_u1 , x_u2 ) + f_2 (x_l1 ) + f_3 (x_u2 , x_l2 )
// where
// x_u = x = (x_u1 , x_u2 )
// x_l = y = (x_l1 , x_l2 )

void SMD_settings(int D_ul, int D_ll, int *settings, int fnum){
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


void SMD_ranges(int D_ul, int D_ll, double *bounds_ul, double *bounds_ll, int fnum){
    int settings[6];

    SMD_settings(D_ul, D_ll, settings, fnum);
    
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
    } else if (fnum == 12) { // from matlab version which is inconsistent with the paper
        // but correct for experiments
        ul1_a = -5.0; ul1_b = 10.0;
        ul2_a = -1.0; ul2_b = 1.0;

        ll1_a = -5.0; ll1_b = 10.0;
        ll2_a = -PI / 4.0 + EPS; ll2_b = PI / 4.0 - EPS;
    }else{
        ul1_a = 0.0; ul1_b = 0.0; ul2_a = 0.0; ul2_b = 0.0; ll1_a = 0.0;
        ll1_b = 0.0; ll2_a = 0.0; ll2_b=0.0;
    }

    for (i = 0; i < p; ++i) {
        bounds_ul[i] = ul1_a;
        bounds_ul[D_ul + i] = ul1_b;
    }
    
    for (i = p; i < p + r; ++i){
        bounds_ul[i] = ul2_a;
        bounds_ul[D_ul + i] = ul2_b;
    }

    for (i = 0; i < q+s; ++i){
        bounds_ll[i] = ll1_a;
        bounds_ll[D_ll + i] = ll1_b;
    }

    for (i = q+s; i < q+r+s; ++i){
        bounds_ll[i] = ll2_a;
        bounds_ll[D_ll + i] = ll2_b;
    }
}


void SMD_solutions(int D_ul, int D_ll, double *x, double *y, int fnum){

    int settings[6], i;
    SMD_settings(D_ul, D_ll, settings, fnum);
    int p = settings[0], q = settings[1], r = settings[2];

    for (i = 0; i < D_ul; ++i)x[i] = 0;
    for (i = 0; i < D_ll; ++i) y[i] = 0;

    if (fnum == 2 || fnum == 7){
        for (i = q; i < D_ll; ++i) y[i] = 1;
    
    } else if (fnum == 5 || fnum == 8) {
        for (i = 0; i < D_ll; ++i) y[i] = 0;
        for (i = 0; i < q; ++i) y[i] = 1; 
    
    } else if (fnum == 10) {
        for (i = 0; i < D_ul; ++i) x[i] = 1.0 / sqrt( (double) p + r - 1);
        for (i = 0; i < q; ++i) y[i] = 1.0 / sqrt((double) q-1);
        for (i = 0; i < r; ++i) y[q+i] = atan( 1.0 / sqrt( (double) p + r - 1) );
    
    } else if (fnum == 11) {
        for (i = 0; i < r; ++i) y[q + i] = exp(-1.0 / sqrt((double) r));
    
    } else if (fnum == 12) {
        for (i = 0; i < D_ul; ++i) x[i] = 1.0 / sqrt( (double) p+r-1);
        for (i = 0; i < q; ++i) y[i] = 1.0 / sqrt( (double) q-1);
        for (i = 0; i < r; ++i) y[q+i] = atan( (1.0 / sqrt( (double) p+r-1)) - (1.0 / sqrt( (double) r)));
        
    }
}

void SMD_optimum(double* F, double* f, int fnum)
{
    if (1 <= fnum && fnum <= 9) {
        F[0] = 0.0; f[0] = 0.0;
        return;
    }else if (fnum == 10) {
        F[0] = 4.0; f[0] = 3.0;
        return;
    }else if (fnum == 11) {
        F[0] = -1.0; f[0] = 1.0;
        return;
    }else if (fnum == 12) {
        F[0] = 3.0; f[0] = 4.0;
        return;
    }

}


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

void SMD9_leader(int p, int q, int r, double *x, double *y, double *F, double *G){
    int i;
    double a = 1.0, b = 1.0;

    double *x_u1 = x, *x_u2 = &x[p];
    double *x_l1 = y, *x_l2 = &y[q];

    double F1 =  sphere(x_u1, p);
    double F2 = -sphere(x_l1, q);
    double F3 =  0.0;

    for (i = 0; i < r; ++i) {
        F3 += pow(x_u2[i], 2) - pow( x_u2[i] - log( 1.0 + x_l2[i] ), 2);
    }

    F[0] = F1 + F2 + F3;

    G[0] = (sphere(x_u1, p) + sphere(x_u2, r)) / a;
    G[0] -= floor( G[0] + 0.5 / b);
    G[0] *= -1.0;
}

void SMD9_follower(int p, int q, int r, double *x, double *y, double *f, double *g){
    int i;
    double a = 1.0, b = 1.0;

    double *x_u1 = x, *x_u2 = &x[p];
    double *x_l1 = y, *x_l2 = &y[q];

    double f1 = sphere(x_u1, p);
    double f2 = sphere(x_l1, q);
    double f3 =  0.0;

    for (i = 0; i < r; ++i) {
        f3 += pow( x_u2[i] - log( 1.0 + x_l2[i] ), 2);
    }

    f[0] = f1 + f2 + f3;

    g[0] = (sphere(x_l1, p) + sphere(x_l2, r)) / a;
    g[0] -= floor( g[0] + 0.5 / b);
    
    g[0] *= -1.0;
}

void SMD10_leader(int p, int q, int r, double *x, double *y, double *F, double *G){
    int i;

    double *x_u1 = x, *x_u2 = &x[p];
    double *x_l1 = y, *x_l2 = &y[q];

    double F1 = 0.0;
    double F2 = sphere(x_l1, q);
    double F3 = 0.0;

    for (i = 0; i < p; ++i){
        F1 += pow(x_u1[i] - 2.0, 2);
    }

    for (i = 0; i < r; ++i) {
        F3 += pow(x_u2[i] - 2.0, 2) - pow( x_u2[i] - tan( x_l2[i] ), 2);
    }

    F[0] = F1 + F2 + F3;

    double xu1_3[p];
    double xu2_3[r];

    for (i = 0; i < p; ++i) { xu1_3[i] = pow(x_u1[i], 3); }
    for (i = 0; i < r; ++i) { xu2_3[i] = pow(x_u2[i], 3); }

    double sum_xu1_3 = 0.0;
    for (i = 0; i < p; ++i) { sum_xu1_3 += xu1_3[i]; }

    double sum_xu2_3 = 0.0;
    for (i = 0; i < r; ++i) { sum_xu2_3 += xu2_3[i]; }

    for (i = 0; i < p; ++i) {
        G[i] = x_u1[i] - (sum_xu1_3 - xu1_3[i]) - sum_xu2_3;
        G[i] *= -1.0;
    }

    for (i = 0; i < r; ++i) {
        G[p + i] = x_u2[i] - (sum_xu2_3 - xu2_3[i]) - sum_xu1_3;
        G[p + i] *= -1.0;
    }

}

void SMD10_follower(int p, int q, int r, double *x, double *y, double *f, double *g){
    int i;

    double *x_u1 = x, *x_u2 = &x[p];
    double *x_l1 = y, *x_l2 = &y[q];

    double f1 = sphere(x_u1, p);
    double f2 = 0.0;
    double f3 = 0.0;

    for (i = 0; i < q; ++i){
        f2 += pow(x_l1[i] - 2.0, 2);
    }

    for (i = 0; i < r; ++i) {
        f3 += pow( x_u2[i] - tan( x_l2[i] ), 2);
    }

    f[0] = f1 + f2 + f3;

    double xl1_3[p];

    for (i = 0; i < q; ++i) { xl1_3[i] = pow(x_l1[i], 3); }

    double sum_xl1_3 = 0.0;
    for (i = 0; i < q; ++i) { sum_xl1_3 += xl1_3[i]; }


    for (i = 0; i < q; ++i) {
        g[i] = x_l1[i] - (sum_xl1_3 - xl1_3[i]);
        g[i] *= -1.0;
    }


}

void SMD11_leader(int p, int q, int r, double *x, double *y, double *F, double *G){
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

    double a = 1.0 / sqrt((double) r);
    
    for (i = 0; i < r; ++i) {
        G[i] = x_u2[i] - a - log(x_l2[i]);
        G[i] *= -1.0;
    }
}

void SMD11_follower(int p, int q, int r, double *x, double *y, double *f, double *g){
    int i;

    double *x_u1 = x, *x_u2 = &x[p];
    double *x_l1 = y, *x_l2 = &y[q];

    double f1 = sphere(x_u1, p);
    double f2 = sphere(x_l1, q);
    double f3 = 0.0;

    for (i = 0; i < r; ++i) {
        f3 += pow( x_u2[i] - log( x_l2[i] ), 2);
    }

    f[0] = f1 + f2 + f3;

    g[0] = -1.0;
    for (i = 0; i < r; ++i) {
        g[0] += pow(x_u2[i] - log(x_l2[i]) ,2);
    }
    g[0] *= -1.0;
}
void SMD12_leader(int p, int q, int r, double *x, double *y, double *F, double *G){
    int i;

    double *x_u1 = x, *x_u2 = &x[p];
    double *x_l1 = y, *x_l2 = &y[q];

    double F1 = 0.0;
    double F2 = sphere(x_l1, q);
    double F3 = 0.0;

    for (i = 0; i < p; ++i){
        F1 += pow(x_u1[i] - 2.0, 2);
    }

    for (i = 0; i < r; ++i) {
        F3 += pow(x_u2[i] - 2.0, 2) + tan(fabs(x_l2[i])) - pow( x_u2[i] - tan( x_l2[i] ), 2);
    }

    F[0] = F1 + F2 + F3;

    for (i = 0; i < r; ++i) {
        G[i] = x_u2[i] - tan(x_l2[i]);
        G[i] *= -1.0;
    }

    //////////////////////////////////////////////////////////
    double xu1_3[p];
    double xu2_3[p];

    for (i = 0; i < p; ++i) { xu1_3[i] = pow(x_u1[i], 3); }
    for (i = 0; i < r; ++i) { xu2_3[i] = pow(x_u2[i], 3); }

    double sum_xu1_3 = 0.0;
    for (i = 0; i < p; ++i) { sum_xu1_3 += xu1_3[i]; }

    double sum_xu2_3 = 0.0;
    for (i = 0; i < r; ++i) { sum_xu2_3 += xu2_3[i]; }
    //////////////////////////////////////////////////////////

    for (i = 0; i < p; ++i) {
        G[r+i] = x_u1[i] - (sum_xu1_3 - xu1_3[i]) - sum_xu2_3;
        G[r+i] *= -1.0;
    }

    p +=r ;
    for (i = 0; i < r; ++i) {
        G[p + i] = x_u2[i] - (sum_xu2_3 - xu2_3[i]) - sum_xu1_3;
        G[p + i] *= -1.0;
    }
}

void SMD12_follower(int p, int q, int r, double *x, double *y, double *f, double *g){
    int i;

    double *x_u1 = x, *x_u2 = &x[p];
    double *x_l1 = y, *x_l2 = &y[q];

    double f1 = sphere(x_u1, p);
    double f2 = 0.0;
    double f3 = 0.0;

    for (i = 0; i < q; ++i){
        f2 += pow(x_l1[i] - 2.0, 2);
    }

    for (i = 0; i < r; ++i) {
        f3 += pow( x_u2[i] - tan( x_l2[i] ), 2);
    }

    f[0] = f1 + f2 + f3;

    double xl1_3[q];

    for (i = 0; i < q; ++i) { xl1_3[i] = pow(x_l1[i], 3); }

    double sum_xl1_3 = 0.0;
    for (i = 0; i < q; ++i) { sum_xl1_3 += xl1_3[i]; }


    g[0] = f3-1.0;
    g[0] *= -1.0;
    for (i = 1; i <= q; ++i) {
        g[i] = x_l1[i] - (sum_xl1_3 - xl1_3[i]);
        g[i] *= -1.0;
    }
}

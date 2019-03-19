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
    }

    for (i = 0; i < r; ++i) {
        G[p + i] = x_u2[i] - (sum_xu2_3 - xu2_3[i]) - sum_xu1_3;
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

    g[0] = 0.0;
    for (i = 0; i < r; ++i) {
        g[0] += pow(x_u2[i] - log(x_l2[i]) ,2);
    }
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
    }

    p +=r ;
    for (i = 0; i < r; ++i) {
        G[p + i] = x_u2[i] - (sum_xu2_3 - xu2_3[i]) - sum_xu1_3;
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

    for (i = 1; i <= q; ++i) {
        g[i] = x_l1[i] - (sum_xl1_3 - xl1_3[i]);
    }
}
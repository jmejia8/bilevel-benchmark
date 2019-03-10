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
    for (i = m+1; i < n; ++i){P += fabs(x[i]);}

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

    for (i = m+1; i < n; ++i){p += pow(y[i], 2);}

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

    for (i = m+1; i < n; ++i){P += pow(x[i], 2);}

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

    for (i = m+1; i < n; ++i){p += pow(y[i], 2);}

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
    for (i = m+1; i < n; ++i){
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

    for (i = m+1; i < n; ++i){p += pow(y[i], 2);}

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

    for (i = m+1; i < n; ++i){
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

    for (i = m+1; i < n; ++i){
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
    for (i = m+1; i < n; ++i){
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

    for (i = m+1; i < n; ++i){
        p += pow(y[i], 2);
    }

    f[0] = p + q*Q;
}
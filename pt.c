#include "constants.h"

void TP1_leader(int nx, int ny, double *x, double *y, double *F, double *G){
    double x1 = x[0];
    double x2 = x[1];
    double y1 = y[0];
    double y2 = y[1];    
    //At optima x1=20,x2=5,y1=10,y2=5,fu=-225
    F[0] = pow(x[1]-30.0, 2) + pow(x[2]-20.0, 2) - 20.0*y[1] + 20.0*y[2];
        
    
    //////////////////////////////////////////////////////////////////////////
    // constraints here
    G[0] = 30.0 - x[1] - 2*x[2];
    G[1] = x[1] + x[2] - 25.0;
    //////////////////////////////////////////////////////////////////////////
}
void TP2_leader(int nx, int ny, double *x, double *y, double *F, double *G){

    //At optima x1=0,x2=30,y1=-10,y2=10,fu=0,fl=-100
    double x1 = x[0];
    double x2 = x[1];
    double y1 = y[0];
    double y2 = y[1];
    F[0] = 2*x1 + 2*x2 - 3*y1 - 3*y2 - 60.0;
        
    
    //////////////////////////////////////////////////////////////////////////
    // constraints here
    G[0] = x1 + x2 + y1 - 2*y2 - 40.0;
    //Lower level constraints included at upper level
    G[1] = 10.0 - x1 + 2*y1;
    G[2] = 10.0 - x2 + 2*y2;
    //////////////////////////////////////////////////////////////////////////
}
void TP3_leader(int nx, int ny, double *x, double *y, double *F, double *G){

    //At optima x1=0,x2=2,y1=1.8750,y2=0.9062,fu=18.6787,fl=1.0156
    double x1 = x[0];
    double x2 = x[1];
    double y1 = y[0];
    double y2 = y[1];
    F[0] = -pow(x1, 2) - 3*pow(x2, 2) - 4*y1 + pow(y2, 2);
        
    
    //////////////////////////////////////////////////////////////////////////
    // constraints here
    G[0] = pow(x1, 2) + 2*x2 - 4;
    //Lower level constraints included at upper level
    G[1] = -3-pow(x1, 2)+2*x1 - pow(x2, 2)+2*y1-y2;
    G[2] = 4-x2-3*y1+4*y2;
    //////////////////////////////////////////////////////////////////////////
}
void TP4_leader(int nx, int ny, double *x, double *y, double *F, double *G){

    //At optima x1=0.29,x2=0.70,y1=0,y2=0.27,y3=0.27,fu=29.2,fl=-3.2
    double x1 = x[0];
    double x2 = x[1];
    double y1 = y[0];
    double y2 = y[1];
    double y3 = y[2];
    F[0] = -8*x1 - 4*x2 + 4*y1 - 40*y2 - 4*y3;
        
    
    //////////////////////////////////////////////////////////////////////////
    // constraints here
    //G = [];
    //Lower level constraints included at upper level
    G[0] = y2+y3-y1-1;
    G[1] = 2*x1-y1+2*y2-0.5*y3-1;
    G[2] = 2*x2+2*y1-y2-0.5*y3-1;
    //////////////////////////////////////////////////////////////////////////
}
void TP5_leader(int nx, int ny, double *x, double *y, double *F, double *G){

    //At optima x1=2,x2=0,y1=2,y2=0,fu=3.6,fl=2
    double x1 = x[0];
    double x2 = x[1];
    double y1 = y[0];
    double y2 = y[1];
        
    double r = 0.1;
    // x = [x1 x2]'; //'
    // y = [y1 y2]'; //'
    
    F[0] = r*(x1*x1 + x2*x2) - 3*y1 - 4*y2 + 0.5*(y1*y1 + y2*y2);
    
    //////////////////////////////////////////////////////////////////////////
    // constraints here
    //G = [];
    //Lower level constraints included at upper level
    G[0] = -0.333*y1 + y2 - 2.0;
    G[1] = y1 - 0.333*y2 -2;
    //////////////////////////////////////////////////////////////////////////
}
void TP6_leader(int nx, int ny, double *x, double *y, double *F, double *G){

    //At optima x1=1.888,y1=0.888,y2=0,fu=1.2098,fl=-7.61
    double x1 = x[0];
    double y1 = y[0];
    double y2 = y[1];

    F[0] = pow((x1-1), 2) + 2*y1 - 2*x1;
        
    
    //////////////////////////////////////////////////////////////////////////
    // constraints here
    //G = [];
    //Lower level constraints included at upper level
    G[0] = 4*x1+5*y1+4*y2-12;
    G[1] = 4*y2-4*x1-5*y1+4;
    G[2] = 4*x1-4*y1+5*y2-4;
    G[3] = 4*y1-4*x1+5*y2-4;
    //////////////////////////////////////////////////////////////////////////
}
void TP7_leader(int nx, int ny, double *x, double *y, double *F, double *G){

    //At optima x1=5.24,x2=5.63,y1=5.24,y2=0,y3=0.27,fu=2.0714,fl=-2.0714
    double x1 = x[0];
    double x2 = x[1];
    double y1 = y[0];
    double y2 = y[1];

    F[0] = -(x1+y1) * (x2+y2) / (1+x1 * y1+x2 * y2);
        
    
    //////////////////////////////////////////////////////////////////////////
    // constraints here
    G[0] = pow(x1, 2)+pow(x2, 2) - 100.0;
    G[1] = x1-x2; //New constraint added at upper level to correct the bilevel problem
    //Lower level constraints included at upper level
    G[2] = y1-x1;
    G[3] = y2-x2;
    //////////////////////////////////////////////////////////////////////////
}
void TP8_leader(int nx, int ny, double *x, double *y, double *F, double *G){

    //At optima x1=0,x2=30,y1=-10,y2=10,fu=0,fl=-100
    double x1 = x[0];
    double x2 = x[1];
    double y1 = y[0];
    double y2 = y[1];
    
    F[0] = abs(2*x1+2*x2-3*y1-3*y2-60);
        
    
    //////////////////////////////////////////////////////////////////////////
    // constraints here
    G[0] = x1+x2+y1-2*y2-40.0;
    //Lower level constraints included at upper level
    G[1] = 2*y1-x1 + 10.0;
    G[2] = 2*y2-x2 + 10.0;
    //////////////////////////////////////////////////////////////////////////
}
void TP9_leader(int nx, int ny, double *x, double *y, double *F){

    //At optima x1=0,x2=30,y1=-10,y2=10,fu=0,fl=-100
    
    F[0] = 0.0;

    int i;
    for (i = 0; i < nx; ++i) {
        F[0] += fabs(x[i] - 1.0);
    }

    for (i = 0; i < ny; ++i) {
        F[0] += fabs(y[i]);
    }  
}
void TP10_leader(int nx, int ny, double *x, double *y, double *F){
    F[0] = 0.0;

    int i;
    for (i = 0; i < nx; ++i) {
        F[0] += fabs(x[i] - 1.0);
    }

    for (i = 0; i < ny; ++i) {
        F[0] += fabs(y[i]);
    }
}


void TP1_follower(int nx, int ny, double *x, double *y, double *f){
    f[0] = pow(x[1] - y[1], 2) + pow(x[2] - y[2], 2);
}
void TP2_follower(int nx, int ny, double *x, double *y, double *f, double *g){

    double x1 = x[0];
    double x2 = x[1];
    double y1 = y[0];
    double y2 = y[1];
    f[0] = pow(y1-x1+20, 2) +  pow(y2-x2+20, 2);


    //////////////////////////
    // constraints here
    g[0] = 10-x1+2*y1;
    g[1] = 10-x2+2*y2;
    //////////////////////////
}
void TP3_follower(int nx, int ny, double *x, double *y, double *f, double *g){
    
    double x1 = x[0];
    double x2 = x[1];
    double y1 = y[0];
    double y2 = y[1];
    f[0] = 2*pow(x1, 2)+pow(y1, 2)-5*y2;


    //////////////////////////
    // constraints here
    g[0] = -3-pow(x1, 2)+2*x1-pow(x2, 2)+2*y1-y2;
    g[1] = 4-x2-3*y1+4*y2;
    //////////////////////////
}
void TP4_follower(int nx, int ny, double *x, double *y, double *f, double *g){
    
    double x1 = x[0];
    double x2 = x[1];
    double y1 = y[0];
    double y2 = y[1];
    double y3 = y[2];
    f[0] = x1+2*x2+y1+y2+2*y3;


    //////////////////////////
    // constraints here
    g[0] = y2+y3-y1-1;
    g[1] = 2*x1-y1+2*y2-0.5*y3-1;
    g[2] = 2*x2+2*y1-y2-0.5*y3-1;
    //////////////////////////
}
void TP5_follower(int nx, int ny, double *x, double *y, double *f, double *g){
    
    double x1 = x[0];
    double x2 = x[1];
    double y1 = y[0];
    double y2 = y[1];
    
    // H = [1 3; 3 10];
    // b = [-1 2; 3 -3];
    // x = [x1 x2]';
    // y = [y1 y2]';
    
    // f[0] = 0.5*y'*H*y + (b*x)'*y;


    //////////////////////////
    // constraints here
    g[0] = -0.333*y1 + y2 - 2.0;
    g[1] = y1 - 0.333*y2 -2;
    //////////////////////////
}
void TP6_follower(int nx, int ny, double *x, double *y, double *f, double *g){
    
    double x1 = x[0];
    double y1 = y[0];
    double y2 = y[1];
    
    f[0] = pow(2*y1-4, 2) + pow(2*y2-1, 2) + x1*y1;


    //////////////////////////
    // constraints here
    g[0] = 4*x1+5*y1+4*y2-12;
    g[1] = 4*y2-4*x1-5*y1+4;
    g[2] = 4*x1-4*y1+5*y2-4;
    g[3] = 4*y1-4*x1+5*y2-4;
    //////////////////////////
}
void TP7_follower(int nx, int ny, double *x, double *y, double *f, double *g){
    
    double x1 = x[0];
    double x2 = x[1];
    double y1 = y[0];
    double y2 = y[1];

    f[0] = (x1+y1) * (x2+y2) / (1+x1 * y1+x2 * y2);


    //////////////////////////
    // constraints here
    g[0] = y1-x1;
    g[1] = y2-x2;
    //////////////////////////
}
void TP8_follower(int nx, int ny, double *x, double *y, double *f, double *g){
    
    double x1 = x[0];
    double x2 = x[1];
    double y1 = y[0];
    double y2 = y[1];

    f[0] = pow(y1-x1+20, 2)+ pow(y2-x2+20, 2);


    //////////////////////////
    // constraints here
    g[0] = 2*y1-x1 + 10.0;
    g[1] = 2*y2-x2 + 10.0;
    //////////////////////////
}
void TP9_follower(int nx, int ny, double *x, double *y, double *f){
    
    int i;
    double sum_y2 = 0.0;
    double prod_cosy = 1.0;
    for (i = 0; i < ny; ++i){
        sum_y2    += y[i]*y[i];
        prod_cosy *= cos(y[i] / sqrt( (double) (i + 1.0)) );
    }

    double exponent = 1.0 + 1.0/4000.0 * sum_y2 - prod_cosy;
    
    sum_y2 = 0.0;
    for (i = 0; i < nx; ++i){
        sum_y2  += x[i]*x[i];
    }
    exponent = exponent*sum_y2;
    f[0] = exp(exponent);
}
void TP10_follower(int nx, int ny, double *x, double *y, double *f){

    int i;
    double sum_y2 = 0.0;
    double prod_cosy = 1.0;
    for (i = 0; i < ny; ++i){
        sum_y2    += y[i]*y[i] * x[i]*x[i];
        prod_cosy *= cos(y[i]*y[i] * x[i]*x[i] / sqrt( (double) (i + 1.0)) );
    }

    double exponent = 1.0 + 1.0/4000.0 * sum_y2 - prod_cosy;
    
    
    f[0] = exp(exponent);
}
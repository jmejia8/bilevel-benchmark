#include "constants.h"

void TP1_leader(int nx, int ny, double *x, double *y, double *F, double *G){
    double x1 = x[1];
    double x2 = x[2];
    double y1 = y[1];
    double y2 = y[2];    
    //At optima x1=20,x2=5,y1=10,y2=5,fu=-225
    F[0] = pow(x[1]-30.0, 2) + pow(x[2]-20.0, 2) - 20.0*y[1] + 20.0*y[2];
        
    
    //////////////////////////////////////////////////////////////////////////
    //Write the constraints here
    G[1] = 30.0 - x[1] - 2*x[2];
    G[2] = x[1] + x[2] - 25.0;
    //////////////////////////////////////////////////////////////////////////
}
void TP2_leader(int nx, int ny, double *x, double *y, double *F, double *G){

    //At optima x1=0,x2=30,y1=-10,y2=10,fu=0,fl=-100
    double x1 = x[1];
    double x2 = x[2];
    double y1 = y[1];
    double y2 = y[2];
    F[0] = 2*x1 + 2*x2 - 3*y1 - 3*y2 - 60.0;
        
    
    //////////////////////////////////////////////////////////////////////////
    //Write the constraints here
    G[1] = x1 + x2 + y1 - 2*y2 - 40.0;
    //Lower level constraints included at upper level
    G[2] = 10.0 - x1 + 2*y1;
    G[3] = 10.0 - x2 + 2*y2;
    //////////////////////////////////////////////////////////////////////////
}
void TP3_leader(int nx, int ny, double *x, double *y, double *F, double *G){

    //At optima x1=0,x2=2,y1=1.8750,y2=0.9062,fu=18.6787,fl=1.0156
    double x1 = x[1];
    double x2 = x[2];
    double y1 = y[1];
    double y2 = y[2];
    F[0] = -pow(x1, 2) - 3*pow(x2, 2) - 4*y1 + pow(y2, 2);
        
    
    //////////////////////////////////////////////////////////////////////////
    //Write the constraints here
    G[1] = pow(x1, 2) + 2*x2 - 4;
    //Lower level constraints included at upper level
    G[2] = -3-pow(x1, 2)+2*x1 - pow(x2, 2)+2*y1-y2;
    G[3] = 4-x2-3*y1+4*y2;
    //////////////////////////////////////////////////////////////////////////
}
void TP4_leader(int nx, int ny, double *x, double *y, double *F, double *G){

    //At optima x1=0.29,x2=0.70,y1=0,y2=0.27,y3=0.27,fu=29.2,fl=-3.2
    double x1 = x[1];
    double x2 = x[2];
    double y1 = y[1];
    double y2 = y[2];
    double y3 = y[3];
    F[0] = -8*x1 - 4*x2 + 4*y1 - 40*y2 - 4*y3;
        
    
    //////////////////////////////////////////////////////////////////////////
    //Write the constraints here
    //G = [];
    //Lower level constraints included at upper level
    G[1] = y2+y3-y1-1;
    G[2] = 2*x1-y1+2*y2-0.5*y3-1;
    G[3] = 2*x2+2*y1-y2-0.5*y3-1;
    //////////////////////////////////////////////////////////////////////////
}
void TP5_leader(int nx, int ny, double *x, double *y, double *F, double *G){

    //At optima x1=2,x2=0,y1=2,y2=0,fu=3.6,fl=2
    double x1 = x[1];
    double x2 = x[2];
    double y1 = y[1];
    double y2 = y[2];
        
    double r = 0.1;
    // x = [x1 x2]'; //'
    // y = [y1 y2]'; //'
    
    F[0] = r*(x1*x1 + x2*x2) - 3*y1 - 4*y2 + 0.5*(y1*y1 + y2*y2);
    
    //////////////////////////////////////////////////////////////////////////
    //Write the constraints here
    //G = [];
    //Lower level constraints included at upper level
    G[1] = -0.333*y1 + y2 - 2;
    G[2] = y1 - 0.333*y2 -2;
    //////////////////////////////////////////////////////////////////////////
}
 void TP6_leader(int nx, int ny, double *x, double *y, double *F, double *G){

    //At optima x1=1.888,y1=0.888,y2=0,fu=1.2098,fl=-7.61
    double x1 = x[1];
    double y1 = y[1];
    double y2 = y[2];

    F[0] = pow((x1-1), 2) + 2*y1 - 2*x1;
        
    
    //////////////////////////////////////////////////////////////////////////
    //Write the constraints here
    //G = [];
    //Lower level constraints included at upper level
    G[1] = 4*x1+5*y1+4*y2-12;
    G[2] = 4*y2-4*x1-5*y1+4;
    G[3] = 4*x1-4*y1+5*y2-4;
    G[4] = 4*y1-4*x1+5*y2-4;
    //////////////////////////////////////////////////////////////////////////
}
void TP7_leader(int nx, int ny, double *x, double *y, double *F, double *G){

    //At optima x1=5.24,x2=5.63,y1=5.24,y2=0,y3=0.27,fu=2.0714,fl=-2.0714
    double x1 = x[1];
    double x2 = x[2];
    double y1 = y[1];
    double y2 = y[2];

    F[0] = -(x1+y1) * (x2+y2) / (1+x1 * y1+x2 * y2);
        
    
    //////////////////////////////////////////////////////////////////////////
    //Write the constraints here
    G[1] = pow(x1, 2)+pow(x2, 2) - 100;
    G[2] = x1-x2; //New constraint added at upper level to correct the bilevel problem
    //Lower level constraints included at upper level
    G[3] = y1-x1;
    G[4] = y2-x2;
    //////////////////////////////////////////////////////////////////////////
}
void TP8_leader(int nx, int ny, double *x, double *y, double *F, double *G){

    //At optima x1=0,x2=30,y1=-10,y2=10,fu=0,fl=-100
    double x1 = x[1];
    double x2 = x[2];
    double y1 = y[1];
    double y2 = y[2];
    
    F[0] = abs(2*x1+2*x2-3*y1-3*y2-60);
        
    
    //////////////////////////////////////////////////////////////////////////
    //Write the constraints here
    G[1] = x1+x2+y1-2*y2-40;
    //Lower level constraints included at upper level
    G[2] = 2*y1-x1+10;
    G[3] = 2*y2-x2+10;
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

void TP10_leader(int nx, int ny, double *x, double *y, double *F, double *G){
    F[0] = 0.0;

    int i;
    for (i = 0; i < nx; ++i) {
        F[0] += fabs(x[i] - 1.0);
    }

    for (i = 0; i < ny; ++i) {
        F[0] += fabs(y[i]);
    }

}
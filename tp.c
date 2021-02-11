#include <stdio.h>
#include "constants.h"

//  D_ul, D_ll, lenG, leng = settings
void TP_config(int *settings, int fnum){
    if (fnum <= 8) {
        settings[0] = 2;
    }else if (fnum == 9) {
        settings[0] = 5;
    }else if (fnum == 10) {
        settings[0] = 10;
    }else if (fnum == 6) {
        settings[0] = 1;
    }else{
        settings[0] = 0;
    }

    if (fnum == 4) {
        settings[1] = 3;
    }else if (fnum <= 8) {
        settings[1] = 2;
    }else if (fnum == 9) {
        settings[1] = 5;
    }else if (fnum == 10) {
        settings[1] = 10;
    }else{
        settings[1] = 0;
    }

    if (fnum == 1 || fnum == 5) {
        settings[2] = 2;
    }else if (fnum == 2 || fnum == 3 || fnum == 4 || fnum == 8) {
        settings[2] = 3;
    }else if (fnum == 7 || fnum == 6) {
        settings[2] = 4;
    }else{
        settings[2] = 0;
    }

    if (fnum == 2 || fnum == 3 || fnum == 5 || fnum == 7 || fnum == 8) {
        settings[3] = 2;
    }else if (fnum == 4) {
        settings[3] = 3;
    }else if (fnum == 6) {
        settings[3] = 4;
    }else{
        settings[3] = 0;
    }
}

void TP_solutions(int nx, int ny, double *x, double *y, int fnum){
    if (fnum == 1){
        x[0]=20; x[1]=5; y[0]=10; y[1]=5; // fu=-225
    } else if (fnum == 2){
        x[0]=0; x[1]=30; y[0]=-10; y[1]=10; // fu=0; // fl=-100
    } else if (fnum == 3){
        x[0]=0; x[1]=2; y[0]=1.8750; y[1]=0.9062; // fu=18.6787; // fl=1.0156
    } else if (fnum == 4){
        x[0]=0.29; x[1]=0.70; y[0]=0; y[1]=0.27; y[2]=0.27; // fu=29.2; // fl=-3.2
    } else if (fnum == 5){
        x[0]=2; x[1]=0; y[0]=2; y[1]=0; // fu=3.6; // fl=2
    } else if (fnum == 6){
        x[0]=1.888; y[0]=0.888; y[1]=0; // fu=1.2098; // fl=-7.61
    } else if (fnum == 7) {
        x[0]=5.24; x[1]=5.63; y[0]=5.24; y[1]=0; //y[2]=0.27; // fu=2.0714; // fl=-2.0714
    } else if (fnum == 8){
        x[0]=0; x[1]=30; y[0]=-10; y[1]=10; // fu=0; // fl=-100
    } else if (fnum == 9 || fnum == 10) {
        x[0]=0; x[1]=30; y[0]=-10; y[1]=10; // fu=0; // fl=-100
    }
    
}


void TP_optimum(double *F, double *f, int fnum){
    if (fnum == 1){
        F[0] = 225.0; f[0] = 100.0;
    } else if (fnum == 2){
        F[0] = 0.0; f[0] = 100.0;
    } else if (fnum == 3){
        F[0] = -18.6787; f[0] = -1.0156;
    } else if (fnum == 4){
        F[0] = -29.2; f[0] = 3.2;
    } else if (fnum == 5){
        F[0] = -3.6; f[0] = -2.0;
    } else if (fnum == 6){
        F[0] = -1.2098; f[0] = 7.6172;
    } else if (fnum == 7) {
        F[0] = -1.961; f[0] = 1.961;
    } else if (fnum == 8){
        F[0] = 0.0; f[0] = 100.0;
    } else if (fnum == 9 || fnum == 10) {
        F[0] = 0.0; f[0] = 1.0;
    }
    
}

void TP_ranges(int D_ul, int D_ll, double *bounds_ul, double *bounds_ll, int fnum){
    if (fnum == 1){
        bounds_ul[0] = -30.0; bounds_ul[2] = 30.0;
        bounds_ul[1] = -30.0; bounds_ul[3] = 15.0;

        bounds_ll[0] = 0.0;  bounds_ll[2] = 10.0;
        bounds_ll[1] = 0.0;  bounds_ll[3] = 10.0;
    }else if (fnum == 2) {
        bounds_ul[0] = 0.0; bounds_ul[2] = 50.0;
        bounds_ul[1] = 0.0; bounds_ul[3] = 50.0;

        bounds_ll[0] = -10.0;  bounds_ll[2] = 20.0;
        bounds_ll[1] = -10.0;  bounds_ll[3] = 20.0;
    }else if (fnum == 3 || fnum == 5) {
        bounds_ul[0] = 0.0; bounds_ul[2] = 10.0;
        bounds_ul[1] = 0.0; bounds_ul[3] = 10.0;

        bounds_ll[0] = 0.0;  bounds_ll[2] = 10.0;
        bounds_ll[1] = 0.0;  bounds_ll[3] = 10.0;
    }else if (fnum == 4) {
        bounds_ul[0] = 0.0; bounds_ul[2] = 1.0;
        bounds_ul[1] = 0.0; bounds_ul[3] = 1.0;
        int i = 0;
        for (i = 0; i < 3; ++i) {
            bounds_ll[i] = 0.0; bounds_ll[3+i] = 1.0;
        }

    }else if (fnum == 6) {
        bounds_ul[0] = 0.0; bounds_ul[2] = 2.0;
        bounds_ul[1] = 0.0; bounds_ul[3] = 2.0;

        bounds_ll[0] = 0.0;  bounds_ll[2] = 2.0;
        bounds_ll[1] = 0.0;  bounds_ll[3] = 2.0;
    }else if (fnum == 7) {
        bounds_ul[0] = 0.0; bounds_ul[2] = 10.0;
        bounds_ul[1] = 0.0; bounds_ul[3] = 10.0;

        bounds_ll[0] = 0.0;  bounds_ll[2] = 1.0;
        bounds_ll[1] = 0.0;  bounds_ll[3] = 10.0;
    }else if (fnum == 8) {
        bounds_ul[0] = 0.0; bounds_ul[2] = 50.0;
        bounds_ul[1] = 0.0; bounds_ul[3] = 50.0;

        bounds_ll[0] = -10.0;  bounds_ll[2] = 20.0;
        bounds_ll[1] = -10.0;  bounds_ll[3] = 20.0;
    }else if (fnum == 9 ) {
        int i;
        for (i = 0; i < 5; ++i) {
            bounds_ul[i] = -1.0; bounds_ul[5 + i] = 1.0;
            bounds_ll[i] = -PI;  bounds_ll[5 + i] = PI;
        }
    } else if (fnum == 10) {
        int i;
        for (i = 0; i < 10; ++i) {
            bounds_ul[i] = -1.0; bounds_ul[10 + i] = 1.0;
            bounds_ll[i] = -PI;  bounds_ll[10 + i] = PI;
        }
    }
}


void TP1_leader(int nx, int ny, double *x, double *y, double *F, double *G){
    // double x1 = x[0];
    // double x2 = x[1];
    // double y1 = y[0];
    // double y2 = y[1];    
    //At optima x1=20,x2=5,y1=10,y2=5,fu=-225
    F[0] = pow(x[0]-30.0, 2) + pow(x[1]-20.0, 2) - 20.0*y[0] + 20.0*y[1];
        
    
    //////////////////////////////////////////////////////////////////////////
    // constraints here
    G[0] = 30.0 - x[0] - 2.0*x[1];
    G[1] = x[0] + x[1] - 25.0;
    //////////////////////////////////////////////////////////////////////////
}
void TP2_leader(int nx, int ny, double *x, double *y, double *F, double *G){

    //At optima x1=0,x2=30,y1=-10,y2=10,fu=0,fl=-100
    double x1 = x[0];
    double x2 = x[1];
    double y1 = y[0];
    double y2 = y[1];
    F[0] = 2.0*x1 + 2.0*x2 - 3*y1 - 3*y2 - 60.0;
        
    
    //////////////////////////////////////////////////////////////////////////
    // constraints here
    G[0] = x1 + x2 + y1 - 2.0*y2 - 40.0;
    //Lower level constraints included at upper level
    G[1] = 10.0 - x1 + 2.0*y1;
    G[2] = 10.0 - x2 + 2.0*y2;
    //////////////////////////////////////////////////////////////////////////
}
void TP3_leader(int nx, int ny, double *x, double *y, double *F, double *G){

    //At optima x1=0,x2=2,y1=1.8750,y2=0.9062,fu=18.6787,fl=1.0156
    double x1 = x[0];
    double x2 = x[1];
    double y1 = y[0];
    double y2 = y[1];
    F[0] = -pow(x1, 2) - 3*pow(x2, 2) - 4.0*y1 + pow(y2, 2);
        
    
    //////////////////////////////////////////////////////////////////////////
    // constraints here
    G[0] = pow(x1, 2) + 2.0*x2 - 4;
    //Lower level constraints included at upper level
    G[1] = -3-pow(x1, 2)+2.0*x1 - pow(x2, 2)+2.0*y1-y2;
    G[2] = 4-x2-3*y1+4.0*y2;
    //////////////////////////////////////////////////////////////////////////
}
void TP4_leader(int nx, int ny, double *x, double *y, double *F, double *G){

    //At optima x1=0.29,x2=0.70,y1=0,y2=0.27,y3=0.27,fu=29.2,fl=-3.2
    double x1 = x[0];
    double x2 = x[1];
    double y1 = y[0];
    double y2 = y[1];
    double y3 = y[2];
    F[0] = -8.0*x1 - 4.0*x2 + 4.0*y1 - 40.0*y2 - 4.0*y3;
        
    
    //////////////////////////////////////////////////////////////////////////
    // constraints here
    //G = [];
    //Lower level constraints included at upper level
    G[0] = y2+y3-y1-1.0;
    G[1] = 2.0*x1-y1+2.0*y2-0.5*y3-1.0;
    G[2] = 2.0*x2+2.0*y1-y2-0.5*y3-1.0;
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
    
    F[0] = r*(x1*x1 + x2*x2) - 3*y1 - 4.0*y2 + 0.5*(y1*y1 + y2*y2);
    
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

    F[0] = pow((x1-1.0), 2) + 2.0*y1 - 2.0*x1;
        
    
    //////////////////////////////////////////////////////////////////////////
    // constraints here
    //G = [];
    //Lower level constraints included at upper level
    G[0] = 4.0*x1+5.0*y1+4.0*y2-12.0;
    G[1] = 4.0*y2-4.0*x1-5.0*y1+4.0;
    G[2] = 4.0*x1-4.0*y1+5.0*y2-4.0;
    G[3] = 4.0*y1-4.0*x1+5.0*y2-4.0;
    //////////////////////////////////////////////////////////////////////////
}
void TP7_leader(int nx, int ny, double *x, double *y, double *F, double *G){

    //At optima x1=5.24,x2=5.63,y1=5.24,y2=0,y3=0.27,fu=2.0714,fl=-2.0714
    double x1 = x[0];
    double x2 = x[1];
    double y1 = y[0];
    double y2 = y[1];

    F[0] = -(x1+y1) * (x2+y2) / (1.0+x1 * y1+x2 * y2);
        
    
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
    
    F[0] = fabs(2.0*x1+2.0*x2-3*y1-3*y2-60);
        
    
    //////////////////////////////////////////////////////////////////////////
    // constraints here
    G[0] = x1+x2+y1-2.0*y2-40.0;
    //Lower level constraints included at upper level
    G[1] = 2.0*y1-x1 + 10.0;
    G[2] = 2.0*y2-x2 + 10.0;
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
    f[0] = pow(x[1] - y[1], 2) + pow(x[0] - y[0], 2);
}
void TP2_follower(int nx, int ny, double *x, double *y, double *f, double *g){

    double x1 = x[0];
    double x2 = x[1];
    double y1 = y[0];
    double y2 = y[1];
    f[0] = pow(y1-x1+20.0, 2) +  pow(y2-x2+20.0, 2);


    //////////////////////////
    // constraints here
    g[0] = 10-x1+2.0*y1;
    g[1] = 10-x2+2.0*y2;
    //////////////////////////
}
void TP3_follower(int nx, int ny, double *x, double *y, double *f, double *g){
    
    double x1 = x[0];
    double x2 = x[1];
    double y1 = y[0];
    double y2 = y[1];
    f[0] = 2.0*pow(x1, 2)+pow(y1, 2)-5.0*y2;


    //////////////////////////
    // constraints here
    g[0] = -3-pow(x1, 2)+2.0*x1-pow(x2, 2)+2.0*y1-y2;
    g[1] = 4-x2-3*y1+4.0*y2;
    //////////////////////////
}
void TP4_follower(int nx, int ny, double *x, double *y, double *f, double *g){
    
    double x1 = x[0];
    double x2 = x[1];
    double y1 = y[0];
    double y2 = y[1];
    double y3 = y[2];
    f[0] = x1+2.0*x2+y1+y2+2.0*y3;


    //////////////////////////
    // constraints here
    g[0] = y2+y3-y1-1.0;
    g[1] = 2.0*x1-y1+2.0*y2-0.5*y3-1;
    g[2] = 2.0*x2+2.0*y1-y2-0.5*y3-1;
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
    
    f[0] =  y1 * (- x1 + 2.0*x2) + y1 * (0.5*y1 + 1.5*y2) + y2 * (3.0*x1 - 3.0*x2) + y2 * (1.5*y1 + 5.0*y2);


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
    
    f[0] = pow(2.0*y1-4, 2) + pow(2.0*y2-1, 2) + x1*y1;


    //////////////////////////
    // constraints here
    g[0] = 4.0*x1+5.0*y1+4.0*y2-12;
    g[1] = 4.0*y2-4.0*x1-5.0*y1+4;
    g[2] = 4.0*x1-4.0*y1+5.0*y2-4;
    g[3] = 4.0*y1-4.0*x1+5.0*y2-4;
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

    f[0] = pow(y1-x1+20.0, 2)+ pow(y2-x2+20.0, 2);


    //////////////////////////
    // constraints here
    g[0] = 2.0*y1-x1 + 10.0;
    g[1] = 2.0*y2-x2 + 10.0;
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

void TP_test(int fnum)
{
    printf("========================\n");
    printf("Test TP%d\n", fnum);
    printf("========================\n");
    int settings[4];
    TP_config(settings, fnum);
    int D_ul = settings[0], D_ll = settings[1];
    int lenG = settings[2], leng = settings[3];

    double x[D_ul], y[D_ll], F[1], f[1], F_true[1], f_true[1];

    TP_optimum(F, f, fnum);

    // printf("F[0] = %g\n", F[0]);
    // printf("f[0] = %g\n", f[0]);
    // printf("_-_\n");
    double G[lenG], g[leng];
    TP_solutions(0,0, x, y, fnum);
    TP_leader(0,0,x,y,F_true,G, fnum);
    TP_follower(0,0,x,y,f_true,g, fnum);

    if (fabs(F_true[0] - F[0]) > 1e-1) {
        printf("Error TP%d_leader\n", fnum);
        printf("F[0] = %g\n", F[0]);
        printf("F[0] = %g\n", F_true[0]);
    }

    int i;
    printf("G = ");
    for (i = 0; i < lenG; ++i) {
        printf("%g ", G[i]);
    }
    printf("\ng = ");

    for (i = 0; i < leng; ++i) {
        printf("%g ", g[i]);
    }
    printf("\n");
    
    if (fabs(f_true[0] - f[0]) > 1e-1) {
        printf("Error TP%d_follwer\n", fnum);
        printf("f[0] = %g\n", f[0]);
        printf("f[0] = %g\n", f_true[0]);
    }



}

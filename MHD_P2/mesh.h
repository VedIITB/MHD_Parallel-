
#ifndef MESH_H_
#define MESH_H_

#define  GAMMA                1.67
#define  L                    1.0
#define  B                    1.0

#define  mu                   1
#define  ghost                2

#define  xcells                2000
#define  xpoints               xcells + ghost+1
#define  xsize                 xpoints-1

#define  ycells                2000
#define  ypoints               ycells+ ghost +1
#define  ysize                 ypoints-1




#define    elipson    1e-06

 
using namespace std;


     int i,j;

     double  lamdaxMAX;
     double  lamdayMAX;
     double  lamdaMAX;
 


  struct primitive
 {
    double **rho;
    double **vx;
    double **vy;
    double **p;
    double **h; 
    double **vz;   
    double **bx;   
    double **by;   
    double **bz;   
    double **b;   
    double **a;  

 };

 struct waves
 {
    double **fwx;   
    double **fwy;   
    double **afx;   
    double **afy;   
  

 };

 struct conservative
 {
    double  **U1;
    double  **U2;
    double  **U3;
    double  **U4;
    double  **U5;
    double  **U6;
    double  **U7;
    double  **U8;
   
     };

struct flux
 {
    double  **F1;
    double  **F2;
    double  **F3;
    double  **F4;
    double  **F5;
    double  **F6;
    double  **F7;
    double  **F8;


 
    double  **G1;
    double  **G2;
    double  **G3;
    double  **G4;
    double  **G5;
    double  **G6;
    double  **G7;
    double  **G8;
   
     };



struct cellinfo
 {
     double  **lx1;
     double  **lx2;
     double  **ly1;
     double  **ly2;
     double  **lMaxx;
     double  **lMaxy;
     double  **lx;
     double  **ly;
  
     };


struct RungeKutta
 {
    double  **u1;
    double  **u2;
    double  **u3;
    double  **u4;
    double  **u5;
    double  **u6;
    double  **u7;
    double  **u8;
  
    
   };


 struct Powell
 {
    double  **s1;
    double  **s2;
    double  **s3;
    double  **s4;
    double  **s5;
    double  **s6;
    double  **s7;
    double  **s8;
  
    
   };


primitive    q;
conservative u;
flux         f;
cellinfo     c;
RungeKutta   k;
RungeKutta   k2;
Powell       S;
waves        w;
   
 



 #endif


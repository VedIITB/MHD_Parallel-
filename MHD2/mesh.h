
#ifndef MESH_H_
#define MESH_H_

#define  GAMMA                1.67
#define  L                    1.0
#define  B                    1.0

#define  mu                   1
#define  ghost                6

#define  xcells              200
#define  xpoints               xcells + ghost+1
#define  xsize                 xpoints-1

#define  ycells               200
#define  ypoints               ycells+ ghost +1
#define  ysize                 ypoints-1



#define max_no_variable         8


 
using namespace std;


     int i,j;

     double  lamdaxMAX;
     double  lamdayMAX;
     double  lamdaMAX;
 


  struct primitive
 {
    double rho[ypoints][xpoints];
    double vx[ypoints][xpoints];
    double vy[ypoints][xpoints];
    double p[ypoints][xpoints];
    double h[ypoints][xpoints]; 
    double vz[ypoints][xpoints];   
    double bx[ypoints][xpoints];   
    double by[ypoints][xpoints];   
    double bz[ypoints][xpoints];   
    double b[ypoints][xpoints];   
    double a[ypoints][xpoints];  

 };

 struct waves
 {
    double fwx[ypoints][xpoints];   
    double fwy[ypoints][xpoints];   
    double afx[ypoints][xpoints];   
    double afy[ypoints][xpoints];   
  

 };

 struct conservative
 {
    double  U1[ypoints][xpoints];
    double  U2[ypoints][xpoints];
    double  U3[ypoints][xpoints];
    double  U4[ypoints][xpoints];
    double  U5[ypoints][xpoints];
    double  U6[ypoints][xpoints];
    double  U7[ypoints][xpoints];
    double  U8[ypoints][xpoints];
   
     };

struct flux
 {
    double  F1[ypoints][xpoints];
    double  F2[ypoints][xpoints];
    double  F3[ypoints][xpoints];
    double  F4[ypoints][xpoints];
    double  F5[ypoints][xpoints];
    double  F6[ypoints][xpoints];
    double  F7[ypoints][xpoints];
    double  F8[ypoints][xpoints];


 
    double  G1[ypoints][xpoints];
    double  G2[ypoints][xpoints];
    double  G3[ypoints][xpoints];
    double  G4[ypoints][xpoints];
    double  G5[ypoints][xpoints];
    double  G6[ypoints][xpoints];
    double  G7[ypoints][xpoints];
    double  G8[ypoints][xpoints];
   
     };



struct cellinfo
 {
    double  lx1[ypoints][xpoints];
    double  lx2[ypoints][xpoints];
    double  ly1[ypoints][xpoints];
    double  ly2[ypoints][xpoints];
    double  lMaxx[ypoints][xpoints];
    double  lMaxy[ypoints][xpoints];
    double  lx[ypoints][xpoints];
    double  ly[ypoints][xpoints];
  
     };


struct RungeKutta
 {
    double  u1[ypoints][xpoints];
    double  u2[ypoints][xpoints];
    double  u3[ypoints][xpoints];
    double  u4[ypoints][xpoints];
   double  u5[ypoints][xpoints];
   double  u6[ypoints][xpoints];
   double  u7[ypoints][xpoints];
   double  u8[ypoints][xpoints];
  
    
   };

struct Powell
 {
    double  s1[ypoints][xpoints];
    double  s2[ypoints][xpoints];
    double  s3[ypoints][xpoints];
    double  s4[ypoints][xpoints];
   double  s5[ypoints][xpoints];
   double  s6[ypoints][xpoints];
   double  s7[ypoints][xpoints];
   double  s8[ypoints][xpoints];
  
    
   };



primitive    q;
conservative u;
flux         f,Fr;
cellinfo     c;
RungeKutta   k;
RungeKutta   k2;
Powell       S;
waves        w;
   
 



 #endif


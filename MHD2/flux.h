#include "mesh.h"
#include "reconstruction.h"

#define  W    0.02

  struct Reconstruct
 {

  double UR1[ypoints][xpoints];
  double UR2[ypoints][xpoints];
  double UR3[ypoints][xpoints];
  double UR4[ypoints][xpoints];
  double UR5[ypoints][xpoints];
  double UR6[ypoints][xpoints];
  double UR7[ypoints][xpoints];
  double UR8[ypoints][xpoints];
  double hR[ypoints][xpoints];
  double pR[ypoints][xpoints];

  double UL1[ypoints][xpoints];
  double UL2[ypoints][xpoints];
  double UL3[ypoints][xpoints];
  double UL4[ypoints][xpoints];
  double UL5[ypoints][xpoints];
  double UL6[ypoints][xpoints];
  double UL7[ypoints][xpoints];
  double UL8[ypoints][xpoints];
  double hL[ypoints][xpoints];
  double pL[ypoints][xpoints];
 };



Reconstruct *X= new Reconstruct ;
Reconstruct *Y =new Reconstruct;




void flux1D( double U2[][xpoints],  double U1[][xpoints],
           double U3[][xpoints],
           double U4[][xpoints],
           double U5[][xpoints], 
           double U6[][xpoints], 
           double U7[][xpoints] , 
           double U8[][xpoints] , double p[][xpoints],double h[][xpoints]);

void flux_WENO( double U2[][xpoints],  double U1[][xpoints],
           double U3[][xpoints],
           double U4[][xpoints],
           double U5[][xpoints], 
           double U6[][xpoints], 
           double U7[][xpoints] , 
           double U8[][xpoints] , double p[][xpoints],double h[][xpoints]);


void flux_WENO3( double U2[][xpoints],  double U1[][xpoints],
           double U3[][xpoints],
           double U4[][xpoints],
           double U5[][xpoints], 
           double U6[][xpoints], 
           double U7[][xpoints] , 
           double U8[][xpoints] , double p[][xpoints],double h[][xpoints]);




void flux1D( double U2[][xpoints],
           double U1[][xpoints],
           double U3[][xpoints],
           double U4[][xpoints],
           double U5[][xpoints], 
           double U6[][xpoints], 
           double U7[][xpoints] , 
           double U8[][xpoints] ,double  p[][xpoints], double h[][xpoints]) {

for(i=2; i<=ysize-3; i++)
        { 
 	
          for(j=2; j<=xsize-3; j++)
            {


  f.F1[i][j] = (( U2[i][j]+ U2[i][j+1]) - ( U1[i][j+1]-U1[i][j])*c.lx[i][j])*0.5;
  f.F2[i][j] = 0.5*(U2[i][j]*U2[i][j]/U1[i][j]+ U2[i][j+1]*U2[i][j+1]/U1[i][j+1])  + 0.5*(p[i][j]+p[i][j+1])+ 

  0.25*(((U6[i][j]*U6[i][j])+(U7[i][j]*U7[i][j]))+(U6[i][j+1]*U6[i][j+1])+(U7[i][j+1]*U7[i][j+1])+(U8[i][j]*U8[i][j])+(U8[i][j+1]*U8[i][j+1])) -((U7[i][j]*U7[i][j])+(U7[i][j+1]*U7[i][j+1]))/2 - (U2[i][j+1]-U2[i][j])*c.lx[i][j]*0.5;
                  
 f.F3[i][j] = 0.5*(U2[i][j]*U3[i][j]/U1[i][j] + U2[i][j+1]*U3[i][j+1]/U1[i][j+1])  -((U6[i][j]*U7[i][j])+U6[i][j+1]*U7[i][j+1])*0.5 - (U3[i][j+1]-U3[i][j])*c.lx[i][j]*0.5; 

 f.F4[i][j] = 0.5*(U2[i][j]*U4[i][j]/U1[i][j]+ U2[i][j+1]*U4[i][j+1]/U1[i][j+1])  -((U8[i][j]*U7[i][j])+(U8[i][j+1]*U7[i][j+1]))*0.5- (U4[i][j+1]-U4[i][j])*c.lx[i][j]*0.5; 

f.F5[i][j] = 0.5*(U2[i][j]*h[i][j]+ U2[i][j+1]*h[i][j+1])- (U7[i][j]*(U7[i][j]*U2[i][j]/U1[i][j]+U6[i][j]*U3[i][j]/U1[i][j]+ U8[i][j]*U4[i][j]/U1[i][j])+ (U7[i][j+1]*(U7[i][j+1]*U2[i][j+1]/U1[i][j+1]+U6[i][j+1]*U3[i][j+1]/U1[i][j+1]+ U8[i][j+1]*U4[i][j+1]/U1[i][j+1])))*0.5 - (U5[i][j+1]-U5[i][j])*c.lx[i][j]*0.5;
	
f.F6[i][j]=((U6[i][j]*U2[i][j]/U1[i][j])+(U6[i][j+1]*U2[i][j+1]/U1[i][j+1]))*.5 - ((U7[i][j]*U3[i][j]/U1[i][j])+(U7[i][j+1]*U3[i][j+1]/U1[i][j+1]))*.5-(U6[i][j+1]-U6[i][j])*c.lx[i][j]*0.5; 

f.F7[i][j]= -(U7[i][j+1]-U7[i][j])*c.lx[i][j]*0.5; 
		   
f.F8[i][j]= ((U8[i][j]*U2[i][j]/U1[i][j])+(U8[i][j+1]*U2[i][j+1]/U1[i][j+1]))*.5 - ((U7[i][j]*U4[i][j]/U1[i][j])+(U7[i][j+1]*U4[i][j+1]/U1[i][j+1]))*.5-(U8[i][j+1]-U8[i][j])*c.lx[i][j]*0.5;


//------------------------------------------------------------------------------------------------------------------------------------//
//                                                     FOR Y-DIRECTION                                                    
//------------------------------------------------------------------------------------------------------------------------------------//

f.G1[i][j] = ((U1[i][j]*U3[i][j]/U1[i][j] + U1[i+1][j]*U3[i+1][j]/U1[i+1][j]) - (U1[i+1][j]-U1[i][j])*c.ly[i][j])*0.5;       
		  
f.G2[i][j] = 0.5*(U3[i][j]*U2[i][j]/U1[i][j]+ U3[i+1][j]*U2[i+1][j]/U1[i+1][j]) -((U6[i][j]*U7[i][j])+U6[i+1]
[j]*U7[i+1][j])*0.5 - (U2[i+1][j]-U2[i][j])*c.ly[i][j]*0.5; 

 f.G3[i][j] = 0.5*(U3[i][j]*U3[i][j]/U1[i][j]+ U3[i+1][j]*U3[i+1][j]/U1[i+1][j]+p[i][j]+p[i+1][j]) - (U3[i+1][j]-U3[i][j])*c.ly[i][j]*0.5 +0.25*(((U6[i][j]*U6[i][j])+(U7[i][j]*U7[i][j])+(U8[i][j]*U8[i][j]))+(U6[i+1][j]*U6[i+1][j])+(U7[i+1][j]*U7[i+1][j])+(U8[i+1][j]*U8[i+1][j])) - ((U6[i+1][j]*U6[i+1][j])+(U6[i][j]*U6[i][j]))/2;
 
f.G4[i][j] = 0.5*(U3[i][j]*U4[i][j]/U1[i][j]+ u.U3[i+1][j]*u.U4[i+1][j]/U1[i+1][j]) -((U6[i][j]*U8[i][j])+U6[i+1][j]*U8[i+1][j])*0.5 - (U4[i+1][j]-U4[i][j])*c.ly[i][j]*0.5; 

 f.G5[i][j] = 0.5*(U3[i][j]*h[i][j]+ U3[i+1][j]*h[i+1][j]) - (U5[i+1][j]-U5[i][j])*c.ly[i][j]*0.5 - (.5*(U6[i][j]*(U7[i][j]*U2[i][j]/U1[i][j]+(U6[i][j]*U3[i][j]/U1[i][j])+(U8[i][j]*U4[i][j]/U1[i][j]))+(U6[i+1][j]*(U7[i+1][j]*U2[i+1][j]/U1[i+1][j]+(U6[i+1][j]*U3[i+1][j]/U1[i+1][j])+(U8[i+1][j]*U4[i+1][j]/U1[i+1][j]))))) ;
	   
f.G6[i][j]=-(U6[i+1][j]-U6[i][j])*c.ly[i][j]*0.5; 
 
f.G7[i][j]= ((U7[i][j]*U3[i][j]/U1[i][j])+(U7[i+1][j]*U3[i+1][j]/U1[i+1][j]))*.5 - ((U6[i][j]*U2[i][j]/U1[i][j])+(U6[i+1][j]*U2[i+1][j]/U1[i+1][j]))*.5- (U7[i+1][j]-U7[i][j])*c.ly[i][j]*0.5; 

f.G8[i][j]= ((U8[i][i]*U3[i][j]/U1[i][j])+(U8[i+1][j]*U3[i+1][j]/U1[i+1][j]))*.5 - ((U6[i][j]*U4[i][j]/U1[i][j])+(U6[i+1][j]*U4[i+1][j]/U1[i+1][j]))*.5-(U8[i+1][j]-U8[i][j])*c.ly[i][j]*0.5;
    
 }

 }

}
  
   void flux_WENO( double U2[][xpoints],
           double U1[][xpoints],
           double U3[][xpoints],
           double U4[][xpoints],
           double U5[][xpoints], 
           double U6[][xpoints], 
           double U7[][xpoints] , 
           double U8[][xpoints] , 
           double  p[][xpoints], double h[][xpoints]) {


   
 WENO_X(U1, X->UR1,X->UL1);  //rho
 WENO_X(U2, X->UR2,X->UL2);  // vx
 WENO_X(U3, X->UR3,X->UL3);  // vy
 WENO_X(U4, X->UR4,X->UL4);  // vz
 WENO_X(U5, X->UR5,X->UL5);  // bx
 WENO_X(U6, X->UR6,X->UL6);  // by
 WENO_X(U7, X->UR7,X->UL7);  
 WENO_X(U8, X->UR8,X->UL8);  

 WENO_X(p, X->pR,X->pL);
 WENO_X(h, X->hR,X->hL);

 WENO_Y(U1, Y->UR1,Y->UL1);
 WENO_Y(U2, Y->UR2,Y->UL2);
 WENO_Y(U3, Y->UR3,Y->UL3);
 WENO_Y(U4, Y->UR4,Y->UL4);
 WENO_Y(U5, Y->UR5,Y->UL5);
 WENO_Y(U6, Y->UR6,Y->UL6);
 WENO_Y(U7, Y->UR7,Y->UL7);
 WENO_Y(p, Y->pR,Y->pL);
 WENO_Y(h, Y->hR,Y->hL);
 WENO_Y(U8, Y->UR8,Y->UL8);  


 //-----------------------------------------------------------------------------------------------------------------------------------//

//---------------------------------------------- RECONSTRUCT FLUX --------------------------------------------------------------------//
 
//--------------------------------------------------------------------------------------------------------------------------------------//

for(i=2; i<=ysize-2; i++)
        { 
 	
          for(j=2; j<=xsize-2; j++)
            {




  f.F1[i][j] = (( X->UL2[i][j]+ X->UR2[i][j+1]) - ( X->UR1[i][j+1]-X->UL1[i][j])*c.lx[i][j])*0.5;

  f.F2[i][j] = 0.5*(X->UL2[i][j]*X->UL2[i][j]/X->UL1[i][j]+ X->UR2[i][j+1]*X->UR2[i][j+1]/X->UR1[i][j+1])  + 0.5*(X->pL[i][j]+X->pR[i][j+1])+ 0.25*(((X->UL6[i][j]*X->UL6[i][j])+(X->UL7[i][j]*X->UL7[i][j]))+(X->UR6[i][j+1]*X->UR6[i][j+1])+(X->UR7[i][j+1]*X->UR7[i][j+1])+(X->UL8[i][j]*X->UL8[i][j])+(X->UR8[i][j+1]*X->UR8[i][j+1])) -((X->UL7[i][j]*X->UL7[i][j])+(X->UR7[i][j+1]*X->UR7[i][j+1]))/2 - (X->UR2[i][j+1]-X->UL2[i][j])*c.lx[i][j]*0.5;
                  
 f.F3[i][j] = 0.5*(X->UL2[i][j]*X->UL3[i][j]/X->UL1[i][j] + X->UR2[i][j+1]*X->UR3[i][j+1]/X->UR1[i][j+1])  -((X->UL6[i][j]*X->UL7[i][j])+X->UR6[i][j+1]*X->UR7[i][j+1])*0.5 - (X->UR3[i][j+1]-X->UL3[i][j])*c.lx[i][j]*0.5; 

 f.F4[i][j] = 0.5*(X->UL2[i][j]*X->UL4[i][j]/X->UL1[i][j]+ X->UR2[i][j+1]*X->UR4[i][j+1]/X->UR1[i][j+1])  -((X->UL8[i][j]*X->UL7[i][j])+(X->UR8[i][j+1]*X->UR7[i][j+1]))*0.5- (X->UR4[i][j+1]-X->UL4[i][j])*c.lx[i][j]*0.5; 

f.F5[i][j] = 0.5*(X->UL2[i][j]*X->hL[i][j]+ X->UR2[i][j+1]*X->hR[i][j+1])- (X->UL7[i][j]*(X->UL7[i][j]*X->UL2[i][j]/X->UL1[i][j]+X->UL6[i][j]*X->UL3[i][j]/X->UL1[i][j]+ X->UL8[i][j]*X->UL4[i][j]/X->UL1[i][j])+ (X->UR7[i][j+1]*(X->UR7[i][j+1]*X->UR2[i][j+1]/X->UR1[i][j+1]+X->UR6[i][j+1]*X->UR3[i][j+1]/X->UR1[i][j+1]+ X->UR8[i][j+1]*X->UR4[i][j+1]/X->UR1[i][j+1])))*0.5 - (X->UR5[i][j+1]-X->UL5[i][j])*c.lx[i][j]*0.5;
	
f.F6[i][j]=((X->UL6[i][j]*X->UL2[i][j]/X->UL1[i][j])+(X->UR6[i][j+1]*X->UR2[i][j+1]/X->UR1[i][j+1]))*.5 - ((X->UL7[i][j]*X->UL3[i][j]/X->UL1[i][j])+(X->UR7[i][j+1]*X->UR3[i][j+1]/X->UR1[i][j+1]))*.5-(X->UR6[i][j+1]-X->UL6[i][j])*c.lx[i][j]*0.5; 

f.F7[i][j]= -(X->UR7[i][j+1]-X->UL7[i][j])*c.lx[i][j]*0.5; 
		   
f.F8[i][j]= ((X->UL8[i][j]*X->UL2[i][j]/X->UL1[i][j])+(X->UR8[i][j+1]*X->UR2[i][j+1]/X->UR1[i][j+1]))*.5 - ((X->UL7[i][j]*X->UL4[i][j]/X->UL1[i][j])+(X->UR7[i][j+1]*X->UR4[i][j+1]/X->UR1[i][j+1]))*.5-(X->UR8[i][j+1]-X->UL8[i][j])*c.lx[i][j]*0.5;

//------------------------------------------------------------------------------------------------------------------------------------//
//                                                     FOR Y-DIRECTION                                                    
//------------------------------------------------------------------------------------------------------------------------------------//

f.G1[i][j] = ((Y->UL1[i][j]*Y->UL3[i][j]/Y->UL1[i][j] + Y->UR1[i+1][j]*Y->UR3[i+1][j]/Y->UR1[i+1][j]) - (Y->UR1[i+1][j]-Y->UL1[i][j])*c.ly[i][j])*0.5;       
		  
f.G2[i][j] = 0.5*(Y->UL3[i][j]*Y->UL2[i][j]/Y->UL1[i][j]+ Y->UR3[i+1][j]*Y->UR2[i+1][j]/Y->UR1[i+1][j]) -((Y->UL6[i][j]*Y->UL7[i][j])+Y->UR6[i+1][j]*Y->UR7[i+1][j])*0.5 - (Y->UR2[i+1][j]-Y->UL2[i][j])*c.ly[i][j]*0.5; 

 f.G3[i][j] = 0.5*(Y->UL3[i][j]*Y->UL3[i][j]/Y->UL1[i][j]+ Y->UR3[i+1][j]*Y->UR3[i+1][j]/Y->UR1[i+1][j]+ Y->pL[i][j]+ Y->pR[i+1][j]) - (Y->UR3[i+1][j]-Y->UL3[i][j])*c.ly[i][j]*0.5 +0.25*(((Y->UL6[i][j]*Y->UL6[i][j])+(Y->UL7[i][j]*Y->UL7[i][j])+(Y->UL8[i][j]*Y->UL8[i][j]))+(Y->UR6[i+1][j]*Y->UR6[i+1][j])+(Y->UR7[i+1][j]*Y->UR7[i+1][j])+(Y->UR8[i+1][j]*Y->UR8[i+1][j])) - ((Y->UR6[i+1][j]*Y->UR6[i+1][j])+(Y->UL6[i][j]*Y->UL6[i][j]))/2;
 
f.G4[i][j] = 0.5*(Y->UL3[i][j]*Y->UL4[i][j]/Y->UL1[i][j]+ Y->UR3[i+1][j]*Y->UR4[i+1][j]/Y->UR1[i+1][j]) -((Y->UL6[i][j]*Y->UL8[i][j])+Y->UR6[i+1][j]*Y->UL8[i+1][j])*0.5 - (Y->UR4[i+1][j]-Y->UL4[i][j])*c.ly[i][j]*0.5; 

 f.G5[i][j] = 0.5*(Y->UL3[i][j]*Y->hL[i][j]+ Y->UR3[i+1][j]*Y->hR[i+1][j]) - (Y->UR5[i+1][j]-Y->UL5[i][j])*c.ly[i][j]*0.5 - (.5*(Y->UL6[i][j]*(Y->UL7[i][j]*Y->UL2[i][j]/Y->UL1[i][j]+(Y->UL6[i][j]*Y->UL3[i][j]/Y->UL1[i][j])+(Y->UL8[i][j]*Y->UL4[i][j]/Y->UL1[i][j]))+(Y->UR6[i+1][j]*(Y->UR7[i+1][j]*Y->UR2[i+1][j]/Y->UR1[i+1][j]+(Y->UR6[i+1][j]*Y->UR3[i+1][j]/Y->UR1[i+1][j])+(Y->UR8[i+1][j]*Y->UR4[i+1][j]/Y->UR1[i+1][j]))))) ;
	   
f.G6[i][j]=-(Y->UR6[i+1][j]-Y->UL6[i][j])*c.ly[i][j]*0.5; 
 
f.G7[i][j]= ((Y->UL7[i][j]*Y->UL3[i][j]/Y->UL1[i][j])+(Y->UR7[i+1][j]*Y->UR3[i+1][j]/Y->UR1[i+1][j]))*.5 - ((Y->UL6[i][j]*Y->UL2[i][j]/Y->UL1[i][j])+(Y->UR6[i+1][j]*Y->UR2[i+1][j]/Y->UR1[i+1][j]))*.5- (Y->UR7[i+1][j]-Y->UL7[i][j])*c.ly[i][j]*0.5; 

f.G8[i][j]= ((Y->UL8[i][i]*Y->UL3[i][j]/Y->UL1[i][j])+(Y->UR8[i+1][j]*Y->UR3[i+1][j]/Y->UR1[i+1][j]))*.5 - ((Y->UL6[i][j]*Y->UL4[i][j]/Y->UL1[i][j])+(Y->UR6[i+1][j]*Y->UR4[i+1][j]/Y->UR1[i+1][j]))*.5-(Y->UR8[i+1][j]-Y->UL8[i][j])*c.ly[i][j]*0.5;
    


   }
 }

}

void flux_WENO3( double U2[][xpoints],
           double U1[][xpoints],
           double U3[][xpoints],
           double U4[][xpoints],
           double U5[][xpoints], 
           double U6[][xpoints], 
           double U7[][xpoints] , 
           double U8[][xpoints] , 
           double  p[][xpoints], double h[][xpoints]) {


   
 WENO_3X(U1, X->UR1,X->UL1);  //rho
 WENO_3X(U2, X->UR2,X->UL2);  // vx
 WENO_3X(U3, X->UR3,X->UL3);  // vy
 WENO_3X(U4, X->UR4,X->UL4);  // vz
 WENO_3X(U5, X->UR5,X->UL5);  // bx
 WENO_3X(U6, X->UR6,X->UL6);  // by
 WENO_3X(U7, X->UR7,X->UL7);  // bz
 WENO_3X(U8, X->UR8,X->UL8);  // bz
 WENO_3X(p, X->pR,X->pL);
 WENO_3X(h, X->hR,X->hL);

 WENO_3Y(U1, Y->UR1,Y->UL1);
 WENO_3Y(U2, Y->UR2,Y->UL2);
 WENO_3Y(U3, Y->UR3,Y->UL3);
 WENO_3Y(U4, Y->UR4,Y->UL4);
 WENO_3Y(U5, Y->UR5,Y->UL5);
 WENO_3Y(U6, Y->UR6,Y->UL6);
 WENO_3Y(U7, Y->UR7,Y->UL7);
 WENO_3Y(U8, Y->UR8,Y->UL8);
 WENO_3Y(p, Y->pR,Y->pL);
 WENO_3Y(h, Y->hR,Y->hL);


 //-----------------------------------------------------------------------------------------------------------------------------------//

//---------------------------------------------- RECONSTRUCT FLUX --------------------------------------------------------------------//
 
//--------------------------------------------------------------------------------------------------------------------------------------//

for(i=2; i<=ysize-2; i++)
        { 
 	
          for(j=2; j<=xsize-2; j++)
            {




  Fr.F1[i][j] = (( X->UL2[i][j]+ X->UR2[i][j+1]) - ( X->UR1[i][j+1]-X->UL1[i][j])*c.lx[i][j])*0.5;

  Fr.F2[i][j] = 0.5*(X->UL2[i][j]*X->UL2[i][j]/X->UL1[i][j]+ X->UR2[i][j+1]*X->UR2[i][j+1]/X->UR1[i][j+1])  + 0.5*(X->pL[i][j]+X->pR[i][j+1])+ 0.25*(((X->UL6[i][j]*X->UL6[i][j])+(X->UL7[i][j]*X->UL7[i][j]))+(X->UR6[i][j+1]*X->UR6[i][j+1])+(X->UR7[i][j+1]*X->UR7[i][j+1])+(X->UL8[i][j]*X->UL8[i][j])+(X->UR8[i][j+1]*X->UR8[i][j+1])) -((X->UL7[i][j]*X->UL7[i][j])+(X->UR7[i][j+1]*X->UR7[i][j+1]))/2 - (X->UR2[i][j+1]-X->UL2[i][j])*c.lx[i][j]*0.5;
                  
 Fr.F3[i][j] = 0.5*(X->UL2[i][j]*X->UL3[i][j]/X->UL1[i][j] + X->UR2[i][j+1]*X->UR3[i][j+1]/X->UR1[i][j+1])  -((X->UL6[i][j]*X->UL7[i][j])+X->UR6[i][j+1]*X->UR7[i][j+1])*0.5 - (X->UR3[i][j+1]-X->UL3[i][j])*c.lx[i][j]*0.5; 

 Fr.F4[i][j] = 0.5*(X->UL2[i][j]*X->UL4[i][j]/X->UL1[i][j]+ X->UR2[i][j+1]*X->UR4[i][j+1]/X->UR1[i][j+1])  -((X->UL8[i][j]*X->UL7[i][j])+(X->UR8[i][j+1]*X->UR7[i][j+1]))*0.5- (X->UR4[i][j+1]-X->UL4[i][j])*c.lx[i][j]*0.5; 

Fr.F5[i][j] = 0.5*(X->UL2[i][j]*X->hL[i][j]+ X->UR2[i][j+1]*X->hR[i][j+1])- (X->UL7[i][j]*(X->UL7[i][j]*X->UL2[i][j]/X->UL1[i][j]+X->UL6[i][j]*X->UL3[i][j]/X->UL1[i][j]+ X->UL8[i][j]*X->UL4[i][j]/X->UL1[i][j])+ (X->UR7[i][j+1]*(X->UR7[i][j+1]*X->UR2[i][j+1]/X->UR1[i][j+1]+X->UR6[i][j+1]*X->UR3[i][j+1]/X->UR1[i][j+1]+ X->UR8[i][j+1]*X->UR4[i][j+1]/X->UR1[i][j+1])))*0.5 - (X->UR5[i][j+1]-X->UL5[i][j])*c.lx[i][j]*0.5;
	
Fr.F6[i][j]=((X->UL6[i][j]*X->UL2[i][j]/X->UL1[i][j])+(X->UR6[i][j+1]*X->UR2[i][j+1]/X->UR1[i][j+1]))*.5 - ((X->UL7[i][j]*X->UL3[i][j]/X->UL1[i][j])+(X->UR7[i][j+1]*X->UR3[i][j+1]/X->UR1[i][j+1]))*.5-(X->UR6[i][j+1]-X->UL6[i][j])*c.lx[i][j]*0.5; 

Fr.F7[i][j]= -(X->UR7[i][j+1]-X->UL7[i][j])*c.lx[i][j]*0.5; 
		   
Fr.F8[i][j]= ((X->UL8[i][j]*X->UL2[i][j]/X->UL1[i][j])+(X->UR8[i][j+1]*X->UR2[i][j+1]/X->UR1[i][j+1]))*.5 - ((X->UL7[i][j]*X->UL4[i][j]/X->UL1[i][j])+(X->UR7[i][j+1]*X->UR4[i][j+1]/X->UR1[i][j+1]))*.5-(X->UR8[i][j+1]-X->UL8[i][j])*c.lx[i][j]*0.5;

//------------------------------------------------------------------------------------------------------------------------------------//
//                                                     FOR Y-DIRECTION                                                    
//------------------------------------------------------------------------------------------------------------------------------------//

Fr.G1[i][j] = ((Y->UL1[i][j]*Y->UL3[i][j]/Y->UL1[i][j] + Y->UR1[i+1][j]*Y->UR3[i+1][j]/Y->UR1[i+1][j]) - (Y->UR1[i+1][j]-Y->UL1[i][j])*c.ly[i][j])*0.5;       
		  
Fr.G2[i][j] = 0.5*(Y->UL3[i][j]*Y->UL2[i][j]/Y->UL1[i][j]+ Y->UR3[i+1][j]*Y->UR2[i+1][j]/Y->UR1[i+1][j]) -((Y->UL6[i][j]*Y->UL7[i][j])+Y->UR6[i+1][j]*Y->UR7[i+1][j])*0.5 - (Y->UR2[i+1][j]-Y->UL2[i][j])*c.ly[i][j]*0.5; 

 Fr.G3[i][j] = 0.5*(Y->UL3[i][j]*Y->UL3[i][j]/Y->UL1[i][j]+ Y->UR3[i+1][j]*Y->UR3[i+1][j]/Y->UR1[i+1][j]+ Y->pL[i][j]+ Y->pR[i+1][j]) - (Y->UR3[i+1][j]-Y->UL3[i][j])*c.ly[i][j]*0.5 +0.25*(((Y->UL6[i][j]*Y->UL6[i][j])+(Y->UL7[i][j]*Y->UL7[i][j])+(Y->UL8[i][j]*Y->UL8[i][j]))+(Y->UR6[i+1][j]*Y->UR6[i+1][j])+(Y->UR7[i+1][j]*Y->UR7[i+1][j])+(Y->UR8[i+1][j]*Y->UR8[i+1][j])) - ((Y->UR6[i+1][j]*Y->UR6[i+1][j])+(Y->UL6[i][j]*Y->UL6[i][j]))/2;
 
Fr.G4[i][j] = 0.5*(Y->UL3[i][j]*Y->UL4[i][j]/Y->UL1[i][j]+ Y->UR3[i+1][j]*Y->UR4[i+1][j]/Y->UR1[i+1][j]) -((Y->UL6[i][j]*Y->UL8[i][j])+Y->UR6[i+1][j]*Y->UL8[i+1][j])*0.5 - (Y->UR4[i+1][j]-Y->UL4[i][j])*c.ly[i][j]*0.5; 

 Fr.G5[i][j] = 0.5*(Y->UL3[i][j]*Y->hL[i][j]+ Y->UR3[i+1][j]*Y->hR[i+1][j]) - (Y->UR5[i+1][j]-Y->UL5[i][j])*c.ly[i][j]*0.5 - (.5*(Y->UL6[i][j]*(Y->UL7[i][j]*Y->UL2[i][j]/Y->UL1[i][j]+(Y->UL6[i][j]*Y->UL3[i][j]/Y->UL1[i][j])+(Y->UL8[i][j]*Y->UL4[i][j]/Y->UL1[i][j]))+(Y->UR6[i+1][j]*(Y->UR7[i+1][j]*Y->UR2[i+1][j]/Y->UR1[i+1][j]+(Y->UR6[i+1][j]*Y->UR3[i+1][j]/Y->UR1[i+1][j])+(Y->UR8[i+1][j]*Y->UR4[i+1][j]/Y->UR1[i+1][j]))))) ;
	   
Fr.G6[i][j]=-(Y->UR6[i+1][j]-Y->UL6[i][j])*c.ly[i][j]*0.5; 
 
Fr.G7[i][j]= ((Y->UL7[i][j]*Y->UL3[i][j]/Y->UL1[i][j])+(Y->UR7[i+1][j]*Y->UR3[i+1][j]/Y->UR1[i+1][j]))*.5 - ((Y->UL6[i][j]*Y->UL2[i][j]/Y->UL1[i][j])+(Y->UR6[i+1][j]*Y->UR2[i+1][j]/Y->UR1[i+1][j]))*.5- (Y->UR7[i+1][j]-Y->UL7[i][j])*c.ly[i][j]*0.5; 

Fr.G8[i][j]= ((Y->UL8[i][i]*Y->UL3[i][j]/Y->UL1[i][j])+(Y->UR8[i+1][j]*Y->UR3[i+1][j]/Y->UR1[i+1][j]))*.5 - ((Y->UL6[i][j]*Y->UL4[i][j]/Y->UL1[i][j])+(Y->UR6[i+1][j]*Y->UR4[i+1][j]/Y->UR1[i+1][j]))*.5-(Y->UR8[i+1][j]-Y->UL8[i][j])*c.ly[i][j]*0.5;
    


   }
 }


 for(i=2; i<=ysize-2; i++)
        { 
 	
          for(j=2; j<=xsize-2; j++)
            {

     f.F1[i][j] = f.F1[i][j]*W + (1-W)*Fr.F1[i][j]; 
     f.F2[i][j] = f.F2[i][j]*W + (1-W)*Fr.F2[i][j];
     f.F3[i][j] = f.F3[i][j]*W + (1-W)*Fr.F3[i][j];
     f.F4[i][j] = f.F4[i][j]*W + (1-W)*Fr.F4[i][j]; 
     f.F5[i][j] = f.F5[i][j]*W + (1-W)*Fr.F5[i][j];
     f.F6[i][j] = f.F6[i][j]*W + (1-W)*Fr.F6[i][j]; 
 
     f.F7[i][j] = f.F7[i][j]*W + (1-W)*Fr.F7[i][j]; 
     f.F8[i][j] = f.F8[i][j]*W + (1-W)*Fr.F8[i][j]; 

    f.G1[i][j] = f.G1[i][j]*W + (1-W)*Fr.G1[i][j]; 
    f.G2[i][j] = f.G2[i][j]*W + (1-W)*Fr.G2[i][j]; 
    f.G3[i][j] = f.G3[i][j]*W + (1-W)*Fr.G3[i][j]; 
    f.G4[i][j] = f.G4[i][j]*W + (1-W)*Fr.G4[i][j]; 
    f.G5[i][j] = f.G5[i][j]*W + (1-W)*Fr.G5[i][j]; 
    f.G6[i][j] = f.G6[i][j]*W + (1-W)*Fr.G6[i][j]; 
    f.G7[i][j] = f.G7[i][j]*W + (1-W)*Fr.G7[i][j]; 
    f.G8[i][j] = f.G8[i][j]*W + (1-W)*Fr.G8[i][j]; 
  

  }


}



}

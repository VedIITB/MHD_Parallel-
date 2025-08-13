#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fstream>
#include <iostream>
#include "mesh.h"
#include <time.h>
#include "initialize.h"
#include "reconstruction.h"
#include "flux.h"


bool WENO=true;

bool Third=true;

void initilization();
void time_step();
void iteration();
void eigen_value();
void maxEigen();
void Rusanov();
void RK3();    
void tracer();
void boundaryConditions();
void output();
void input();
void conserve();


double dt=1e-3 , dx, dy , t=0.0 , CFL , tf ,  bc1 ,bc2 ,bc3 , bc4 , maxD=0;

double x[xpoints],y[ypoints], start, stop;



using namespace std;

   int main()
  {


    dx = L/xcells;
    dy = B/ycells;
   
start = clock();
    input();
    initilization();
	

         do{

        
	     RK3();
	     time_step();
	    
	
          }
	   while(t<=tf);
stop = clock();   	  

   delete X,Y;//Fr;
   output();
	  
   cout<<"THE END OF THE PROGRAMME IN \t"<<(stop-start)/CLOCKS_PER_SEC<<"\tsecs"<<endl;

   


    }

     void initilization()
     {

      

        int l=0,m=0;
   
     for(j=3; j<=xsize-3; j++)
          {
 	
          
	     x[j]=  -0.0+ dx*(l);
           l++;
 
   
          } 
  
  
       for(i=3; i<=ysize-3; i++)
         {
  
 
             y[i]= -0.0 + dy*(m);
           m++;
            
         }  
  
	Primitive( q.p,  q.rho, q.vx, q.vy, q.vz, q.bx,q.by,q.bz, x,y );


 
     for(i=3; i<=ysize-3; i++)
   {
 	
    for(j=3; j<=xsize-3; j++)
     {
          
    
           q.a[i][j] = sqrt(GAMMA*q.p[i][j]/q.rho[i][j]) ;
	 
	   q.h[i][j] = 0.5*(q.vx[i][j]*q.vx[i][j]) + 0.5*(q.vy[i][j]*q.vy[i][j])+0.5*(q.vz[i][j]*q.vz[i][j]) + (q.a[i][j]*q.a[i][j])/(GAMMA-1)+ ((q.bz[i][j]*q.bz[i][j] + q.by[i][j]*q.by[i][j]+q.bx[i][j]*q.bx[i][j])/q.rho[i][j]);
         
//.........................................................................................................//

           u.U1[i][j] = q.rho[i][j];
      
	   u.U2[i][j] = q.rho[i][j]*q.vx[i][j];
      
	   u.U3[i][j] = q.rho[i][j]*q.vy[i][j];

           u.U4[i][j] = q.rho[i][j]*q.vz[i][j];
       
           u.U5[i][j] = q.p[i][j]/(GAMMA-1) + 0.5*q.rho[i][j]*(q.vy[i][j]*q.vy[i][j]+ q.vx[i][j]*q.vx[i][j]+q.vz[i][j]*q.vz[i][j])         +0.5*(q.bz[i][j]*q.bz[i][j] + q.by[i][j]*q.by[i][j]+q.bx[i][j]*q.bx[i][j]);

           u.U6[i][j] = q.by[i][j];

           u.U7[i][j] = q.bx[i][j];

           u.U8[i][j] = q.bz[i][j];

       }
 
    }

boundaryConditions();

  
  }
void wave_system()
{

double temp1,alfven;
  
     for(i=3; i<=ysize-3; i++)
   {
 	
    for(j=3; j<=xsize-3; j++)
     {
 

    w.afy[i][j]  =  u.U6[i][j]/sqrt(u.U1[i][j]);

   w.afx[i][j]  =  u.U7[i][j]/sqrt(u.U1[i][j]);

   alfven     = sqrt((u.U7[i][j]*u.U7[i][j])+(u.U6[i][j]*u.U6[i][j])+u.U8[i][j]*u.U8[i][j])/sqrt(u.U1[i][j]);

   temp1 = (q.a[i][j]*q.a[i][j])+(alfven*alfven);
  
    
   w.fwx[i][j] = sqrt(0.5*((temp1)+sqrt((temp1*temp1)-(4*q.a[i][j]*q.a[i][j]*w.afx[i][j]*w.afx[i][j]))));
	   
  	

   w.fwy[i][j] = sqrt(0.5*((temp1)+sqrt((temp1*temp1)-(4*q.a[i][j]*q.a[i][j]*w.afy[i][j]*w.afy[i][j]))));
	

 	

 	//cout<<w.fwy[i][j]<<endl; 

     
     }
  
   }
     
 }


void eigen_value()
{
      
   
 for(i=2; i<=ysize-2; i++)
 {
 	
    for(j=2; j<=xsize-2; j++)
     {
 

  
     c.lx1[i][j]  = q.vx[i][j] + w.fwx[i][j];
 
     c.lx2[i][j]  = q.vx[i][j] - w.fwx[i][j];

     c.ly1[i][j]  = q.vy[i][j] + w.fwy[i][j];

     c.ly2[i][j] =  q.vy[i][j] - w.fwy[i][j];
   
 


     }
  
   }
     
 }
  


  
 

 void maxEigen()
   {


   lamdaxMAX=0;
   lamdayMAX=0;
  
      

    for(i=2; i<=ysize-2; i++)
   { 
 	
    for(j=2; j<=xsize-2; j++)
        {


    if(lamdaxMAX<fabs(c.lx1[i][j]))
       lamdaxMAX=fabs(c.lx1[i][j]);

      if(lamdayMAX<fabs(c.ly1[i][j]))
        lamdayMAX=fabs(c.ly1[i][j]);
     
 
       }
  

    }

   
   for(i=2; i<=ysize-2; i++)
     { 
 	
    for(j=2; j<=xsize-2; j++)
        {


     if(lamdaxMAX<fabs(c.lx2[i][j]))
       lamdaxMAX=fabs(c.lx2[i][j]);

     if(lamdayMAX<fabs(c.ly2[i][j]))
        lamdayMAX=fabs(c.ly2[i][j]);
 
        }
  

     } 


  if(lamdaxMAX>lamdayMAX)
  lamdaMAX=lamdaxMAX;

  else 
  lamdaMAX=lamdayMAX;

   }
  	

  void Rusanov()
  {


       for(i=2; i<=ysize-2; i++)
        { 
 	
          for(j=2; j<=xsize-2; j++)
            {


               if(fabs(c.lx2[i][j])>fabs(c.lx1[i][j]))
                 c.lMaxx[i][j] = fabs(c.lx2[i][j]);

               else 
                 c.lMaxx[i][j] = fabs(c.lx1[i][j]);

            }   
  

        }

     

      for(i=2; i<=ysize-2; i++)
        { 
 	
          for(j=2; j<=xsize-2; j++)
            {


                 if(fabs(c.ly2[i][j])>fabs(c.ly1[i][j]))
                 c.lMaxy[i][j] = fabs(c.ly2[i][j]);

               else 
                 c.lMaxy[i][j]= fabs(c.ly1[i][j]);
  
            
 
             }   
  

         }


  

        for(i=2; i<=ysize-3; i++)
        { 
 	
          for(j=2; j<=xsize-3; j++)
            {


              if(c.lMaxx[i][j+1] > c.lMaxx[i][j] )
                 c.lx[i][j] = c.lMaxx[i][j+1];

               else
                 c.lx[i][j] = c.lMaxx[i][j];
 
             }   
  

         }

       for(i=2; i<=ysize-3; i++)
        { 
 	
          for(j=2; j<=xsize-3; j++)
            {



                if(c.lMaxy[i+1][j] > c.lMaxy[i][j] )
                 c.ly[i][j] = c.lMaxy[i+1][j];

               else
                 c.ly[i][j]= c.lMaxy[i][j];
     
 
             }   
  

         }


     
if(WENO==false){
      flux1D( u.U2,
           u.U1,
           u.U3,
           u.U4,
           u.U5, 
           u.U6, 
           u.U7 , 
           u.U8, 
           q.p,q.h);
      }


	 
else{

      flux1D( u.U2,
           u.U1,
           u.U3,
           u.U4,
           u.U5, 
           u.U6, 
           u.U7 , 
           u.U8, 
           q.p,q.h);
      


	flux_WENO3( u.U2,
           u.U1,
           u.U3,
           u.U4,
           u.U5, 
           u.U6, 
           u.U7 , 
           u.U8, 
           q.p,q.h);
      }
  


      
   
   
         
double NetdivBx, NetdivBy;

   for(i=2; i<=ysize-3; i++)
        { 
 	
          for(j=2; j<=xsize-3; j++)
            {

  

    

                NetdivBx= (u.U7[i][j+1]-u.U7[i][j])/dx; 
   
                NetdivBy =  (u.U6[i+1][j]-u.U6[i][j])/dy;


                S.s1[i][j]=0;
        
		S.s1[i][j] =(NetdivBx+ NetdivBy)*u.U7[i][j]; 
		
		S.s3[i][j] =(NetdivBx+ NetdivBy)*u.U6[i][j]; 
		 
		
		S.s4[i][j] =(NetdivBx+ NetdivBy)*u.U8[i][j]; 
		 
		
		S.s5[i][j] =((u.U7[i][j]*u.U2[i][j]/u.U1[i][j])+(u.U6[i][j]*u.U3[i][j]/u.U1[i][j])+(u.U8[i][j]*u.U4[i][j]/u.U1[i][j]))*(NetdivBx +NetdivBx); 
		
		S.s8[i][j] = u.U4[i][j]/u.U1[i][j]*(NetdivBx+ NetdivBy);  
		
		S.s7[i][j] = u.U2[i][j]/u.U1[i][j]*(NetdivBx+ NetdivBy);
        
		S.s6[i][j] =u.U3[i][j]/u.U1[i][j]*(NetdivBx+ NetdivBy);
           
        if(maxD< fabs(NetdivBx+NetdivBy))
          maxD=fabs(NetdivBx+NetdivBy);
  
          }

       }
    
   }
          
   
     

 
void RK3()
{


   iteration();

    for(i=3; i<=ysize-3; i++)
        { 
 	
          for(j=3; j<=xsize-3; j++)
            {



  
       
         k.u1[i][j]  = u.U1[i][j]  + (dt/dx)*(f.F1[i][j-1] -f.F1[i][j])+(dt/dy)*(f.G1[i-1][j] -f.G1[i][j])-dt*S.s1[i][j]; 

         k.u2[i][j]  = u.U2[i][j]  + (dt/dx)*(f.F2[i][j-1] -f.F2[i][j])+(dt/dy)*(f.G2[i-1][j] -f.G2[i][j])-dt*S.s2[i][j];

         k.u3[i][j]  = u.U3[i][j]  + (dt/dx)*(f.F3[i][j-1] -f.F3[i][j])+(dt/dy)*(f.G3[i-1][j] -f.G3[i][j])-dt*S.s3[i][j];

         k.u4[i][j]  = u.U4[i][j]  + (dt/dx)*(f.F4[i][j-1] -f.F4[i][j])+(dt/dy)*(f.G4[i-1][j] -f.G4[i][j])-dt*S.s4[i][j];
   
         k.u5[i][j]  = u.U5[i][j]  + (dt/dx)*(f.F5[i][j-1] -f.F5[i][j])+(dt/dy)*(f.G5[i-1][j] -f.G5[i][j])-dt*S.s5[i][j];
    
         k.u6[i][j]  = u.U6[i][j]  + (dt/dx)*(f.F6[i][j-1] -f.F6[i][j])+(dt/dy)*(f.G6[i-1][j] -f.G6[i][j])-dt*S.s6[i][j];
   
         k.u7[i][j]  = u.U7[i][j]  + (dt/dx)*(f.F7[i][j-1] -f.F7[i][j])+(dt/dy)*(f.G7[i-1][j] -f.G7[i][j])-dt*S.s7[i][j];
   
         k.u8[i][j]  = u.U8[i][j]  + (dt/dx)*(f.F8[i][j-1] -f.F8[i][j])+(dt/dy)*(f.G8[i-1][j] -f.G8[i][j])-dt*S.s8[i][j];


         q.rho[i][j]   =  k.u1[i][j];
     
         q.vx[i][j]    =  k.u2[i][j]/k.u1[i][j];

         q.vy[i][j]    =  k.u3[i][j]/k.u1[i][j];

         q.vz[i][j]    =  k.u4[i][j]/k.u1[i][j];

         q.bz[i][j]=    k.u8[i][j];

         q.by[i][j]=    k.u6[i][j];
 
         q.bx[i][j]=    k.u7[i][j];

      
         q.p[i][j]     = (GAMMA-1)*(k.u5[i][j] - (0.5*k.u2[i][j]*k.u2[i][j])/k.u1[i][j] - 0.5*(k.u3[i][j]*k.u3[i][j])/k.u1[i][j]-0.5*(k.u4[i][j]*k.u4[i][j])/k.u1[i][j]-0.5*(k.u8[i][j]*k.u8[i][j]+k.u6[i][j]*k.u6[i][j] + k.u7[i][j]*k.u7[i][j])) ;
   
         q.a[i][j]    =  sqrt(GAMMA*q.p[i][j]/q.rho[i][j]);
	
	 q.h[i][j] =  0.5*(q.vx[i][j]*q.vx[i][j]) + 0.5*(q.vy[i][j]*q.vy[i][j])+0.5*(q.vz[i][j]*q.vz[i][j])+ (q.a[i][j]*q.a[i][j])/(GAMMA-1)+ ((q.bz[i][j]*q.bz[i][j] + q.by[i][j]*q.by[i][j]+q.bx[i][j]*q.bx[i][j])/q.rho[i][j]);
 
    
 
        }


  }

boundaryConditions();

iteration();
 
   for(i=3; i<=ysize-3; i++)
        { 
 	
          for(j=3; j<=xsize-3; j++)
            {

 
                

   k2.u1[i][j]  = 0.75*u.U1[i][j]  +0.25*k.u1[i][j]  + (0.25*dt/dx)*(f.F1[i][j-1] -f.F1[i][j])+(0.25*dt/dy)*(f.G1[i-1][j] -f.G1[i][j])-0.25*dt*S.s1[i][j] ;
 
     k2.u2[i][j]  = 0.75*u.U2[i][j]  +0.25*k.u2[i][j]  + (0.25*dt/dx)*(f.F2[i][j-1] -f.F2[i][j])+(0.25*dt/dy)*(f.G2[i-1][j] -f.G2[i][j])-dt*S.s2[i][j]*0.25 ;
 
     k2.u3[i][j]  = 0.75*u.U3[i][j]  +0.25*k.u3[i][j]  + (0.25*dt/dx)*(f.F3[i][j-1] -f.F3[i][j])+(0.25*dt/dy)*(f.G3[i-1][j] -f.G3[i][j])-dt*S.s3[i][j]*0.25 ;
  
     k2.u4[i][j]  = 0.75*u.U4[i][j]  +0.25*k.u4[i][j]  + (0.25*dt/dx)*(f.F4[i][j-1] -f.F4[i][j])+(0.25*dt/dy)*(f.G4[i-1][j] -f.G4[i][j])-dt*S.s4[i][j]*0.25;
 
 k2.u5[i][j]  = 0.75*u.U5[i][j]  +0.25*k.u5[i][j]  + (0.25*dt/dx)*(f.F5[i][j-1] -f.F5[i][j])+(0.25*dt/dy)*(f.G5[i-1][j] -f.G5[i][j])-dt*S.s5[i][j]*0.25;
 
  k2.u6[i][j]  = 0.75*u.U6[i][j]  +0.25*k.u6[i][j]  + (0.25*dt/dx)*(f.F6[i][j-1] -f.F6[i][j])+(0.25*dt/dy)*(f.G6[i-1][j] -f.G6[i][j])-dt*S.s6[i][j]*0.25;
       
  k2.u7[i][j]  = 0.75*u.U7[i][j]  +0.25*k.u7[i][j]  + (0.25*dt/dx)*(f.F7[i][j-1] -f.F7[i][j])+(0.25*dt/dy)*(f.G7[i-1][j] -f.G7[i][j])-dt*S.s7[i][j]*0.25;
 
 k2.u8[i][j]  = 0.75*u.U8[i][j]  +0.25*k.u8[i][j]  + (0.25*dt/dx)*(f.F8[i][j-1] -f.F8[i][j])+(0.25*dt/dy)*(f.G8[i-1][j] -f.G8[i][j])-dt*S.s8[i][j]*0.25;
 
         q.rho[i][j]   =  k2.u1[i][j];
     
         q.vx[i][j]    =  k2.u2[i][j]/k2.u1[i][j];

         q.vy[i][j]    =  k2.u3[i][j]/k2.u1[i][j];

         q.vz[i][j]    =  k2.u4[i][j]/k2.u1[i][j];

         q.bz[i][j]=    k2.u8[i][j];

         q.by[i][j]=    k2.u6[i][j];
 
         q.bx[i][j]=    k2.u7[i][j];

      
         q.p[i][j]     = (GAMMA-1)*(k2.u5[i][j] - (0.5*k2.u2[i][j]*k2.u2[i][j])/k2.u1[i][j] - 0.5*(k2.u3[i][j]*k2.u3[i][j])/k2.u1[i][j]-0.5*(k2.u4[i][j]*k2.u4[i][j])/k2.u1[i][j]-0.5*(k2.u8[i][j]*k2.u8[i][j]+k2.u6[i][j]*k2.u6[i][j] + k2.u7[i][j]*k2.u7[i][j])) ;
   
         q.a[i][j]    =  sqrt(GAMMA*q.p[i][j]/q.rho[i][j]);
	
	 q.h[i][j] =  0.5*(q.vx[i][j]*q.vx[i][j]) + 0.5*(q.vy[i][j]*q.vy[i][j])+ 0.5*(q.vz[i][j]*q.vz[i][j]) +(q.a[i][j]*q.a[i][j])/(GAMMA-1)+
	  
               ((q.bz[i][j]*q.bz[i][j] + q.by[i][j]*q.by[i][j]+q.bx[i][j]*q.bx[i][j])/q.rho[i][j]);
 
   
   

           }

  
      }
   

boundaryConditions();
 iteration();
 

   for(i=3; i<=ysize-3; i++)
        { 
 	
          for(j=3; j<=xsize-3; j++)
            {


    u.U1[i][j]  = 0.33*u.U1[i][j]  +0.67*k2.u1[i][j]  + (0.67*dt/dx)*(f.F1[i][j-1] -f.F1[i][j])+(0.67*dt/dy)*(f.G1[i-1][j] -f.G1[i][j]) 
-dt*S.s1[i][j]*0.67;

    u.U2[i][j]  = 0.33*u.U2[i][j]  +0.67*k2.u2[i][j]  + (0.67*dt/dx)*(f.F2[i][j-1] -f.F2[i][j])+(0.67*dt/dy)*(f.G2[i-1][j] -f.G2[i][j])-dt*S.s2[i][j]*0.67; 

    u.U3[i][j]  = 0.33*u.U3[i][j]  +0.67*k2.u3[i][j]  + (0.67*dt/dx)*(f.F3[i][j-1] -f.F3[i][j])+(0.67*dt/dy)*(f.G3[i-1][j] -f.G3[i][j])-
dt*S.s3[i][j]*0.67; 

    u.U4[i][j]  = 0.33*u.U4[i][j]  +0.67*k2.u4[i][j]  + (0.67*dt/dx)*(f.F4[i][j-1] -f.F4[i][j])+(0.67*dt/dy)*(f.G4[i-1][j] -f.G4[i][j])-dt*S.s4[i][j]*0.67; 

    u.U5[i][j]  = 0.33*u.U5[i][j]  +0.67*k2.u5[i][j]  + (0.67*dt/dx)*(f.F5[i][j-1] -f.F5[i][j])+(0.67*dt/dy)*(f.G5[i-1][j] -f.G5[i][j])-dt*S.s5[i][j]*0.67; 

 u.U6[i][j]  = 0.33*u.U6[i][j]  +0.67*k2.u6[i][j]  + (0.67*dt/dx)*(f.F6[i][j-1] -f.F6[i][j])+(0.67*dt/dy)*(f.G6[i-1][j] -f.G6[i][j])-dt*S.s6[i][j]*0.67; 
     
 u.U7[i][j]  = 0.33*u.U7[i][j]  +0.67*k2.u7[i][j]  + (0.67*dt/dx)*(f.F7[i][j-1] -f.F7[i][j])+(0.67*dt/dy)*(f.G7[i-1][j] -f.G7[i][j])-dt*S.s7[i][j]*0.67; 


 u.U8[i][j]  = 0.33*u.U8[i][j]  +0.67*k2.u8[i][j]  + (0.67*dt/dx)*(f.F8[i][j-1] -f.F8[i][j])+(0.67*dt/dy)*(f.G8[i-1][j] -f.G8[i][j])-dt*S.s8[i][j]*0.67; 

         q.rho[i][j]   =  u.U1[i][j];
     
         q.vx[i][j]    =  u.U2[i][j]/u.U1[i][j];

         q.vy[i][j]    =  u.U3[i][j]/u.U1[i][j];

         q.vz[i][j]    =  u.U4[i][j]/u.U1[i][j];

         q.bz[i][j]=    u.U8[i][j];

         q.by[i][j]=    u.U6[i][j];
 
         q.bx[i][j]=    u.U7[i][j];

      
         q.p[i][j]     = (GAMMA-1)*(u.U5[i][j] - (0.5*u.U2[i][j]*u.U2[i][j])/u.U1[i][j] - 0.5*(u.U3[i][j]*u.U3[i][j])/u.U1[i][j]-0.5*(u.U4[i][j]*u.U4[i][j])/u.U1[i][j]-0.5*(u.U8[i][j]*u.U8[i][j]+u.U6[i][j]*u.U6[i][j] + u.U7[i][j]*u.U7[i][j])) ;
   
         q.a[i][j]    =  sqrt(GAMMA*q.p[i][j]/q.rho[i][j]);
	
	 q.h[i][j] =  0.5*(q.vx[i][j]*q.vx[i][j]) + 0.5*(q.vy[i][j]*q.vy[i][j])+0.5*(q.vz[i][j]*q.vz[i][j])+ (q.a[i][j]*q.a[i][j])/(GAMMA-1)+
	  
               ((q.bz[i][j]*q.bz[i][j] + q.by[i][j]*q.by[i][j]+q.bx[i][j]*q.bx[i][j])/q.rho[i][j]);
 
         q.b[i][j] = sqrt(q.bx[i][j]*q.bx[i][j]+ q.by[i][j]*q.by[i][j] + q.bz[i][j]*q.bz[i][j]);
   
     }


  }

boundaryConditions();


   tracer();

}



void tracer()
   {

  for(i=3; i<=ysize-3; i++)
        { 
 	
          for(j=3; j<=xsize-3; j++)
            {


                if(i==50) {
                 cout<<"==============================================================================="<<endl;
                 cout<<"EVOLVING FOR TIME :: "<<t<<endl;
                 cout<<"DENSITY = "<<q.rho[i][j]<<"\n";
                 cout<<"PRESSURE = "<<q.p[i][j]<<"\n";
                 cout<<"VELOCITY = "<<q.vx[i][j]<<"\n";
                 cout<<"MAGNETIC-Y = "<<q.by[i][j]<<"\n";
                 cout<<"MAGNETIC-x = "<<q.bx[i][j]<<"\n";
                 cout<<"LAMDA_MAX ="<<lamdaxMAX<<"\n";
                 cout<<"DIV.B= "<<maxD<<"\n";
                 cout<<"DELTA_T ="<<dt<<"\n";
         
                  }    
   
             }

        }

 }

void time_step()
{


	dt= 0.5*CFL*(dx/lamdaMAX);
         t += dt;


   }
		
 void iteration()
	{
        
	wave_system();
       	eigen_value();
        maxEigen();
        Rusanov();
         
       }	
 
  void output()
  {


     ofstream vx,p,vy,rho,by,Bt;
     vx.open ("v_x.dat", ios_base::out);
     vy.open ("v_y.dat", ios_base::out);
      p.open ("p.dat", ios_base::out);
     rho.open ("rho.dat", ios_base::out);
      by.open ("by.dat", ios_base::out);
    Bt.open ("B.dat", ios_base::out);

 rho<<"TITLE= RHO"<<endl;
 rho<<"Variables = X , Y , Rho"<<endl;
 rho<<"Zone I = "<<xcells+1<<" ,J= "<<ycells+1<<endl;

 p<<"TITLE= PRESSURE"<<endl;
 p<<"Variables = X , Y , p"<<endl;
 p<<"Zone I = "<<xcells+1<<" ,J= "<<ycells+1<<endl;


vx<<"TITLE= velocity-x"<<endl;
 vx<<"Variables = X , Y , vx"<<endl;
 vx<<"Zone I = "<<xcells+1<<" ,J= "<<ycells+1<<endl;

vy<<"TITLE= velocity-y"<<endl;
 vy<<"Variables = X , Y , vy"<<endl;
 vy<<"Zone I = "<<xcells+1<<" ,J= "<<ycells+1<<endl;

by<<"TITLE= q.by-y"<<endl;
 by<<"Variables = X , Y , By"<<endl;
 by<<"Zone I = "<<xcells+1<<" ,J= "<<ycells+1<<endl;

Bt<<"TITLE= TOTAL-B"<<endl;
 Bt<<"Variables = X , Y , B"<<endl;
 Bt<<"Zone I = "<<xcells+1<<" ,J= "<<ycells+1<<endl;

  for(i=3; i<=ysize-3; i++)
        { 
 	
          for(j=3; j<=xsize-3; j++)
            {


             rho<<x[j]<<"\t"<<y[i]<<"\t"<<q.rho[i][j]<<endl;  
              vy<<x[j]<<"\t"<<y[i]<<"\t"<<q.vy[i][j]<<endl;  
              vx<<x[j]<<"\t"<<y[i]<<"\t"<<q.vx[i][j]<<endl;  
              p<<x[j]<<"\t"<<y[i]<<"\t"<<q.p[i][j]<<endl;  
              by<<x[j]<<"\t"<<y[i]<<"\t"<<q.by[i][j]<<endl;  
              Bt<<x[j]<<"\t"<<y[i]<<"\t"<<q.b[i][j]<<endl;  
 

            }

         
       
             vx<<endl;  
             vy<<endl;  
             rho<<endl;  
             p<<endl;  
              
       }
   
   std::string filename = "rho.dat";
    std::string command = "paraview ";
   command += filename;   

   system(command.c_str());

 
   }

 void input()
 {

  string s;  
  ifstream f;
  f.open ("input.dat", ios_base::in);
  f>>s>>tf;
  f>>s>>CFL;
  f>>s>>bc1; 
  f>>s>>bc2;
  f>>s>>bc3;
  f>>s>>bc4;

  f.close();
  
  cout<<"tf = "<<tf<<endl;
  cout<<"CFL = "<<CFL<<endl;

 if(bc1==0)
  cout<<"BC FOR TOP = ZERO GRADIENT"<<endl;
 else if(bc1==1)
  cout<<"BC FOR TOP = PERIODIC"<<endl;
 else if(bc1==2)
 cout<<"BC FOR TOP = REFLECTIVE"<<endl;
      

if(bc2==0)
  cout<<"BC FOR BOTTOM = ZERO GRADIENT"<<endl;
else if(bc2==1)
  cout<<"BC FOR BOTTOM = PERIODIC"<<endl;
else if(bc2==2)
cout<<"BC FOR BOTTOM = REFLECTIVE"<<endl;



if(bc3==0)
  cout<<"BC FOR LEFT = ZERO GRADIENT"<<endl;
else if(bc3==1)
  cout<<"BC FOR LEFT = PERIODIC"<<endl;
else if(bc3==2)
cout<<"BC FOR LEFT= REFLECTIVE"<<endl;


if(bc4==0)
  cout<<"BC FOR RIGHT = ZERO GRADIENT"<<endl;
else if(bc4==1)
  cout<<"BC FOR RIGHT= PERIODIC"<<endl;
else if(bc4==2)
cout<<"BC FOR RIGHT= REFLECTIVE"<<endl;



 }

void boundaryConditions()
{


 for(i=0; i<=ysize; i++)
        { 
 	
          for(j=0; j<=xsize; j++)
            {

               
 
              if(bc1==0)
                {
  
                    q.rho[xsize][j] = q.rho[xsize-3][j] ;           
           
                    q.p[xsize][j] = q.p[xsize-3][j] ;          
          
                    q.vx[xsize][j] = q.vx[xsize-3][j];            
          
                    q.vy[xsize][j] = q.vy[xsize-3][j] ;           
          
                    q.h[xsize][j] = q.h[xsize-3][j]   ;         
          
                        
          
                    q.a[xsize][j] = q.a[xsize-3][j]  ;  
                  
                    q.vz[xsize][j] = q.vz[xsize-3][j]  ;  
           
                    q.bx[xsize][j] = q.bx[xsize-3][j]  ;  
          
                    q.by[xsize][j] = q.by[xsize-3][j]  ;  
                   
                    q.bz[xsize][j] = q.bz[xsize-3][j]  ;  
        
          
                 
                    u.U1[xsize][j] = u.U1[xsize-3][j]  ;          

                    u.U2[xsize][j] = u.U2[xsize-3][j]  ;          
              
                    u.U3[xsize][j] = u.U3[xsize-3][j]  ;          

                    u.U4[xsize][j] = u.U4[xsize-3][j]  ;          
                   
                     u.U5[xsize][j] = u.U5[xsize-3][j]  ;          
                   
                     u.U6[xsize][j] = u.U6[xsize-3][j]  ;          
               
                     u.U7[xsize][j] = u.U7[xsize-3][j]  ;          
                   
                      u.U8[xsize][j] = u.U8[xsize-3][j]  ;          
                   

             //------------------------------------------//


                  q.rho[xsize-2][j] = q.rho[xsize-3][j] ;           
           
                    q.p[xsize-2][j] = q.p[xsize-3][j] ;          
          
                    q.vx[xsize-2][j] = q.vx[xsize-3][j];            
          
                    q.vy[xsize-2][j] = q.vy[xsize-3][j] ;           
          
                    q.h[xsize-2][j] = q.h[xsize-3][j]   ;         
          
                        
          
                    q.a[xsize-2][j] = q.a[xsize-3][j]  ;  
                  
                    q.vz[xsize-2][j] = q.vz[xsize-3][j]  ;  
           
                  q.bx[xsize-2][j] = q.bx[xsize-3][j]  ;  
          
                   q.by[xsize-2][j] = q.by[xsize-3][j]  ;  
                   
                    q.bz[xsize-2][j] = q.bz[xsize-3][j]  ;  
        
          
                 
                    u.U1[xsize-2][j] = u.U1[xsize-3][j]  ;          

                    u.U2[xsize-2][j] = u.U2[xsize-3][j]  ;          
              
                    u.U3[xsize-2][j] = u.U3[xsize-3][j]  ;          

                    u.U4[xsize-2][j] = u.U4[xsize-3][j]  ;          
                   
                     u.U5[xsize-2][j] = u.U5[xsize-3][j]  ;          
                   
                     u.U6[xsize-2][j] = u.U6[xsize-3][j]  ;          
               
                    u.U7[xsize-2][j] = u.U7[xsize-3][j]  ;          
                   
                   u.U8[xsize-2][j] = u.U8[xsize-3][j]  ;          
                             
               //----------------------------------------//
       
                   q.rho[xsize-1][j] = q.rho[xsize-3][j] ;           
           
                    q.p[xsize-1][j] = q.p[xsize-3][j] ;          
          
                    q.vx[xsize-1][j] = q.vx[xsize-3][j];            
          
                    q.vy[xsize-1][j] = q.vy[xsize-3][j] ;           
          
                    q.h[xsize-1][j] = q.h[xsize-3][j]   ;         
          
                        
          
                    q.a[xsize-1][j] = q.a[xsize-3][j]  ;  
                  
                    q.vz[xsize-1][j] = q.vz[xsize-3][j]  ;  
           
                  q.bx[xsize-1][j] = q.bx[xsize-3][j]  ;  
          
                   q.by[xsize-1][j] = q.by[xsize-3][j]  ;  
                   
                    q.bz[xsize-1][j] = q.bz[xsize-3][j]  ;  
        
          
                 
                    u.U1[xsize-1][j] = u.U1[xsize-3][j]  ;          

                    u.U2[xsize-1][j] = u.U2[xsize-3][j]  ;          
              
                    u.U3[xsize-1][j] = u.U3[xsize-3][j]  ;          

                    u.U4[xsize-1][j] = u.U4[xsize-3][j]  ;          
                   
                     u.U5[xsize-1][j] = u.U5[xsize-3][j]  ;          
                   
                     u.U6[xsize-1][j] = u.U6[xsize-3][j]  ;          
               
                    u.U7[xsize-1][j] = u.U7[xsize-3][j]  ;          
                   
                     u.U8[xsize-1][j] = u.U8[xsize-3][j]  ;          
                   
     






                  }
               
             else if(bc1==1)
                {
                    q.p[xsize-1][j] = q.p[5][j]  ;         
                    q.vx[xsize-1][j] =q.vx[5][j]  ;         
                    q.vy[xsize-1][j] =q.vy[5][j]  ;         
                    q.vz[xsize-1][j] =q.vz[5][j]  ;         
                    q.bx[xsize-1][j] =q.bx[5][j]  ;         
                    q.by[xsize-1][j] =q.by[5][j]  ;         
                    q.bz[xsize-1][j] =q.bz[5][j]  ;         
                    q.rho[xsize-1][j] = q.p[5][j]  ;
         
                
                    q.h[xsize-1][j] = q.h[5][j] ;           
                  
                      
                    q.a[xsize-1][j] = q.a[5][j]  ;          
           
                 
                    u.U1[xsize-1][j] = u.U1[5][j]  ;          

                    u.U2[xsize-1][j] = u.U2[5][j]  ;          
              
                    u.U3[xsize-1][j] = u.U3[5][j]  ;          

                    u.U4[xsize-1][j] = u.U4[5][j]  ;          
                   
                     u.U5[xsize-1][j] = u.U5[5][j]  ;          
                   
                     u.U6[xsize-1][j] = u.U6[5][j]  ;          
               
                     u.U7[xsize-1][j] = u.U7[5][j]  ;          
                   
                     u.U8[xsize-1][j] = u.U8[5][j]  ;
                 



                    q.h[xsize-2][j] = q.h[4][j] ;           
                    q.p[xsize-2][j] = q.p[4][j]  ;         
                    q.vx[xsize-2][j] =q.vx[4][j]  ;         
                    q.vy[xsize-2][j] =q.vy[4][j]  ;         
                    q.vz[xsize-2][j] =q.vz[4][j]  ;         
                    q.bx[xsize-2][j] =q.bx[4][j]  ;         
                    q.by[xsize-2][j] =q.by[4][j]  ;         
                    q.bz[xsize-2][j] =q.bz[4][j]  ;  
                    q.a[xsize-2][j] = q.a[4][j]  ;          
                   q.rho[xsize-2][j] = q.rho[4][j]  ;          
           
                 
                    u.U1[xsize-2][j] = u.U1[4][j]  ;          

                    u.U2[xsize-2][j] = u.U2[4][j]  ;          
              
                    u.U3[xsize-2][j] = u.U3[4][j]  ;          

                    u.U4[xsize-2][j] = u.U4[4][j]  ;          
                   
                     u.U5[xsize-2][j] = u.U5[4][j]  ;          
                   
                     u.U6[xsize-2][j] = u.U6[4][j]  ;          
               
                     u.U7[xsize-2][j] = u.U7[4][j]  ;          
                   
                     u.U8[xsize-2][j] = u.U8[4][j]  ;



                    q.p[xsize][j] = q.p[6][j]  ;
                    q.rho[xsize][j] = q.rho[6][j]  ;
         
                     q.vx[xsize][j] =q.vx[6][j]  ;         
                    q.vy[xsize][j] =q.vy[6][j]  ;         
                    q.vz[xsize][j] =q.vz[6][j]  ;         
                    q.bx[xsize][j] =q.bx[6][j]  ;         
                    q.by[xsize][j] =q.by[6][j]  ;         
                    q.bz[xsize][j] =q.bz[6][j]  ;       
                
                
                    q.h[xsize][j] = q.h[6][j] ;           
          
                      
                    q.a[xsize][j] = q.a[6][j]  ;          
           
                 
                    u.U1[xsize][j] = u.U1[6][j]  ;          

                    u.U2[xsize][j] = u.U2[6][j]  ;          
              
                    u.U3[xsize][j] = u.U3[6][j]  ;          

                    u.U4[xsize][j] = u.U4[6][j]  ;          
                   
                    u.U5[xsize][j] = u.U5[6][j]  ;          
                   
                    u.U6[xsize][j] = u.U6[6][j]  ;          
               
                    u.U7[xsize][j] = u.U7[6][j]  ;          
                   
                    u.U8[xsize][j] = u.U8[6][j]  ;          
                   
                  
 
                  }
               /*
             else if (bc1==2)
                {
  
                    q.rho[xsize][j] = q.rho[xsize-1][j]  ;          
           
                    q.p[xsize][j] = q.p[xsize-1][j] ;          
          
                    q.vx[xsize][j] = -q.vx[xsize-1][j];            
          
                    q.vy[xsize][j] = -q.vy[xsize-1][j];            
          
                    q.h[xsize][j] = q.h[xsize-1][j] ;           
          
                    q.M[xsize][j] = q.M[xsize-1][j]  ;          
          
                    q.a[xsize][j] = q.a[xsize-1][j] ;    


                 }
             */
              if(bc4==0)
                {
  
                    q.rho[i][xsize-1] = q.rho[i][xsize-3] ;           
           
                    q.p[i][xsize-1] = q.p[i][xsize-3]  ;         
          
                    q.vx[i][xsize-1] = q.vx[i][xsize-3]  ;          
          
                    q.vy[i][xsize-1] = q.vy[i][xsize-3]  ;          
          
                    q.h[i][xsize-1] = q.h[i][xsize-3]   ;         
          
                       
          
                    q.a[i][xsize-1] = q.a[i][xsize-3]  ;          
          
                    q.vz[i][xsize-1] = q.vz[i][xsize-3]  ;          
                    q.by[i][xsize-1] = q.by[i][xsize-3]  ;          
                    q.bx[i][xsize-1] = q.bx[i][xsize-3]  ;          
                    q.bz[i][xsize-1] = q.bz[i][xsize-3]  ;          
             
                 
                    u.U1[i][xsize-1] = u.U1[i][xsize-3]  ;          
           
                    u.U2[i][xsize-1] = u.U2[i][xsize-3]  ;  
  
                    u.U3[i][xsize-1] = u.U3[i][xsize-3]  ;  

                    u.U4[i][xsize-1] = u.U4[i][xsize-3]  ;  
 
                   u.U5[i][xsize-1] = u.U5[i][xsize-3]  ;  
                   u.U6[i][xsize-1] = u.U6[i][xsize-3]  ;  
                   u.U7[i][xsize-1] = u.U7[i][xsize-3]  ;  
                   u.U8[i][xsize-1] = u.U8[i][xsize-3]  ;  
                //---------------------------------------//
  
                   q.rho[i][xsize-2] = q.rho[i][xsize-3] ;           
           
                    q.p[i][xsize-2] = q.p[i][xsize-3]  ;         
          
                    q.vx[i][xsize-2] = q.vx[i][xsize-3]  ;          
          
                    q.vy[i][xsize-2] = q.vy[i][xsize-3]  ;          
          
                    q.h[i][xsize-2] = q.h[i][xsize-3]   ;         
          
                  //  q.M[i][xsize-2] = q.M[i][xsize-3]  ;          
          
                    q.a[i][xsize-2] = q.a[i][xsize-3]  ;          
          
                     q.vz[i][xsize-2] = q.vz[i][xsize-3]  ;          
                    q.by[i][xsize-2] = q.by[i][xsize-3]  ;          
                    q.bx[i][xsize-2] = q.bx[i][xsize-3]  ;          
                    q.bz[i][xsize-2] = q.bz[i][xsize-3]  ;          
             
                 
                    u.U1[i][xsize-2] = u.U1[i][xsize-3]  ;          
           
                    u.U2[i][xsize-2] = u.U2[i][xsize-3]  ;  
  
                    u.U3[i][xsize-2] = u.U3[i][xsize-3]  ;  

                    u.U4[i][xsize-2] = u.U4[i][xsize-3]  ;  
 
                   u.U5[i][xsize-2] = u.U5[i][xsize-3]  ;  
                   u.U6[i][xsize-2] = u.U6[i][xsize-3]  ;  
                   u.U7[i][xsize-2] = u.U7[i][xsize-3]  ;  
                   u.U8[i][xsize-2] = u.U8[i][xsize-3]  ;  
 

                   //----------------------------------//
 
                   q.rho[i][xsize] = q.rho[i][xsize-3] ;           
           
                    q.p[i][xsize] = q.p[i][xsize-3]  ;         
          
                    q.vx[i][xsize] = q.vx[i][xsize-3]  ;          
          
                    q.vy[i][xsize] = q.vy[i][xsize-3]  ;          
          
                    q.h[i][xsize] = q.h[i][xsize-3]   ;         
          
               
                    q.a[i][xsize] = q.a[i][xsize-3]  ;          
          
                     q.vz[i][xsize] = q.vz[i][xsize-3]  ;          
                    q.by[i][xsize] = q.by[i][xsize-3]  ;          
                    q.bx[i][xsize] = q.bx[i][xsize-3]  ;          
                    q.bz[i][xsize] = q.bz[i][xsize-3]  ;          
             
                 
                    u.U1[i][xsize] = u.U1[i][xsize-3]  ;          
           
                    u.U2[i][xsize] = u.U2[i][xsize-3]  ;  
  
                    u.U3[i][xsize] = u.U3[i][xsize-3]  ;  

                    u.U4[i][xsize] = u.U4[i][xsize-3]  ;  
 
                   u.U5[i][xsize] = u.U5[i][xsize-3]  ;  
                   u.U6[i][xsize] = u.U6[i][xsize-3]  ;  
                   u.U7[i][xsize] = u.U7[i][xsize-3]  ;  
                   u.U8[i][xsize] = u.U8[i][xsize-3]  ;  
 
 
                  }
               
             else if(bc4==1)
                {

            

                    q.p[i][xsize] = q.p[i][6]     ;      
                    q.rho[i][xsize] = q.rho[i][6]     ;      
            q.vx[i][xsize] = q.vx[i][6]     ;      
            q.vy[i][xsize] = q.vy[i][6]     ;      
            q.vz[i][xsize] = q.vz[i][6]     ;      
            q.bx[i][xsize] = q.bx[i][6]     ;      
            q.by[i][xsize] = q.by[i][6]     ;      
             q.bz[i][xsize] = q.bz[i][6]     ;      
            
           
                    q.h[i][xsize] = q.h[i][6]     ;       
          
             
                    q.a[i][xsize] = q.a[i][6]   ;         
           
                     u.U1[i][xsize] = u.U1[i][6]   ;         
                      u.U2[i][xsize] = u.U2[i][6]   ;         
                     u.U3[i][xsize] = u.U3[i][6]   ;         
                    u.U4[i][xsize] = u.U4[i][6]   ;         
                   u.U5[i][xsize] = u.U5[i][6]   ;         
                   u.U6[i][xsize] = u.U6[i][6]   ;         
                   u.U7[i][xsize] = u.U7[i][6]   ;         
                   u.U8[i][xsize] = u.U8[i][6]   ;        

 
           
                     q.p[i][xsize-1] = q.p[i][5]     ;      
                     q.rho[i][xsize-1] = q.rho[i][5] ;
                     q.vx[i][xsize-1] = q.vx[i][5] ;
                   q.vy[i][xsize-1] = q.vy[i][5] ;
          q.bx[i][xsize-1] = q.bx[i][5] ;
          q.vz[i][xsize-1] = q.vz[i][5] ;
          q.by[i][xsize-1] = q.by[i][5] ;
          q.bz[i][xsize-1] = q.bz[i][5] ;
          
                    q.h[i][xsize-1] = q.h[i][5]     ;       
          
                   // q.M[xsize][j] = q.M[1][j]   ;         
          
                    q.a[i][xsize-1] = q.a[i][5]   ;         
           
                     u.U1[i][xsize-1] = u.U1[i][5]   ;         
                      u.U2[i][xsize-1] = u.U2[i][5]   ;         
                     u.U3[i][xsize-1] = u.U3[i][5]   ;         
                    u.U4[i][xsize-1] = u.U4[i][5]   ;         
                   u.U5[i][xsize-1] = u.U5[i][5]   ;         
                   u.U6[i][xsize-1] = u.U6[i][5]   ;         
                   u.U7[i][xsize-1] = u.U7[i][5]   ;         
                   u.U8[i][xsize-1] = u.U8[i][5]   ;         
           


                 
          
                 
                     q.p[i][xsize-2] = q.p[i][4]     ;      
                     q.h[i][xsize-2] = q.h[i][4]     ;       
                     q.a[i][xsize-2] = q.a[i][4]   ;         
                    q.rho[i][xsize-2] = q.rho[i][4]  ;
                     q.vx[i][xsize-2] = q.vx[i][4]  ;
                    q.vy[i][xsize-2] = q.vy[i][4]  ;
                    q.vz[i][xsize-2] = q.vz[i][4]  ;
                    q.bx[i][xsize-2] = q.bx[i][4]  ;
                    q.by[i][xsize-2] = q.by[i][4]  ;
                    q.bz[i][xsize-2] = q.bz[i][4]  ;
                  



 
                      u.U1[i][xsize-2] = u.U1[i][4]   ;         
                      u.U2[i][xsize-2] = u.U2[i][4]   ;         
                     u.U3[i][xsize-2] = u.U3[i][4]   ;         
                    u.U4[i][xsize-2] = u.U4[i][4]   ;         
                   u.U5[i][xsize-2] = u.U5[i][4]   ;         
                   u.U6[i][xsize-2] = u.U6[i][4]   ;         
                   u.U7[i][xsize-2] = u.U7[i][4]   ;         
                   u.U8[i][xsize-2] = u.U8[i][4]   ;         
           
                  }
              /*
             else if (bc4==2)
                {
  
                    q.rho[xsize][j] = q.rho[xsize-1][j] ;           
           
                    q.p[xsize][j] = q.p[xsize-1][j]   ;        
          
                    q.vx[xsize][j] = -q.vx[xsize-1][j] ;           
          
                    q.vy[xsize][j] = -q.vy[xsize-1][j]  ;          
          
                    q.h[xsize][j] = q.h[xsize-1][j]    ;        
          
                    q.M[xsize][j] = q.M[xsize-1][j]    ;        
          
                    q.a[xsize][j] = q.a[xsize-1][j]    ; 


                 }
              */
           if(bc2==0)
                {
  
                    q.rho[0][j] = q.rho[3][j]   ;         
           
                    q.p[0][j] = q.p[3][j]   ;        
          
                    q.vx[0][j] = q.vx[3][j]   ;         
          
                    q.vy[0][j] = q.vy[3][j]   ;         
          
                    q.h[0][j] = q.h[3][j]  ;          
          
             
                    q.a[0][j] = q.a[3][j]   ;    

                     q.vz[0][j] = q.vz[3][j]   ;    
                    q.by[0][j] = q.by[3][j]   ;    
                   q.bx[0][j] = q.bx[3][j]   ;    
                   q.bz[0][j] = q.bz[3][j]   ;    

                    u.U1[0][j] = u.U1[3][j]   ;    
     
                    u.U2[0][j] = u.U2[3][j]   ;    
     
                    u.U3[0][j] = u.U3[3][j]   ;    
     
                    u.U4[0][j] = u.U4[3][j]   ;    
     
                     u.U5[0][j] = u.U5[3][j]   ;    
                     u.U6[0][j] = u.U6[3][j]   ;    
                     u.U7[0][j] = u.U7[3][j]   ;    
                     u.U8[0][j] = u.U8[3][j]   ;    
               //--------------------------------//
                 
                    q.rho[1][j] = q.rho[3][j]   ;         
           
                    q.p[1][j] = q.p[3][j]   ;        
          
                    q.vx[1][j] = q.vx[3][j]   ;         
          
                    q.vy[1][j] = q.vy[3][j]   ;         
          
                    q.h[1][j] = q.h[3][j]  ;          
          
             
                    q.a[1][j] = q.a[3][j]   ;    

                     q.vz[1][j] = q.vz[3][j]   ;    
                    q.by[1][j] = q.by[3][j]   ;    
                   q.bx[1][j] = q.bx[3][j]   ;    
                   q.bz[1][j] = q.bz[3][j]   ;    

                    u.U1[1][j] = u.U1[3][j]   ;    
     
                    u.U2[1][j] = u.U2[3][j]   ;    
     
                    u.U3[1][j] = u.U3[3][j]   ;    
     
                    u.U4[1][j] = u.U4[3][j]   ;    
     
                     u.U5[1][j] = u.U5[3][j]   ;    
                     u.U6[1][j] = u.U6[3][j]   ;    
                     u.U7[1][j] = u.U7[3][j]   ;    
                     u.U8[1][j] = u.U8[3][j]   ;    
     

                  //----------------------------//
                    
                        q.rho[2][j] = q.rho[3][j]   ;         
           
                    q.p[2][j] = q.p[3][j]   ;        
          
                    q.vx[2][j] = q.vx[3][j]   ;         
          
                    q.vy[2][j] = q.vy[3][j]   ;         
          
                    q.h[2][j] = q.h[3][j]  ;          
          
             
                    q.a[2][j] = q.a[3][j]   ;    

                     q.vz[2][j] = q.vz[3][j]   ;    
                    q.by[2][j] = q.by[3][j]   ;    
                   q.bx[2][j] = q.bx[3][j]   ;    
                   q.bz[2][j] = q.bz[3][j]   ;    

                    u.U1[2][j] = u.U1[3][j]   ;    
     
                    u.U2[2][j] = u.U2[3][j]   ;    
     
                    u.U3[2][j] = u.U3[3][j]   ;    
     
                    u.U4[2][j] = u.U4[3][j]   ;    
     
                     u.U5[2][j] = u.U5[3][j]   ;    
                     u.U6[2][j] = u.U6[3][j]   ;    
                     u.U7[2][j] = u.U7[3][j]   ;    
                     u.U8[2][j] = u.U8[3][j]   ;    
     
            }
               
             else if(bc2==1)
                {
  
                
              q.p[0][j] = q.p[xsize-6][j]  ;         
              q.h[0][j] = q.h[xsize-6][j]  ;          
              q.a[0][j] = q.a[xsize-6][j]  ; 
             q.rho[0][j] = q.rho[xsize-6][j]  ;         
              q.vx[0][j] = q.vx[xsize-6][j]  ;          
              q.vy[0][j] = q.vy[xsize-6][j]  ; 
              q.vz[0][j] = q.vz[xsize-6][j]  ;         
              q.bx[0][j] = q.bx[xsize-6][j]  ;          
              q.by[0][j] = q.by[xsize-6][j]  ; 
              q.bz[0][j] = q.bz[xsize-6][j]  ; 
         
             u.U1[0][j] = u.U1[xsize-6][j]  ;          
             u.U2[0][j] = u.U2[xsize-6][j]  ;          
             u.U3[0][j] = u.U3[xsize-6][j]  ;          
             u.U4[0][j] = u.U4[xsize-6][j]  ;          
             u.U5[0][j] = u.U5[xsize-6][j]  ;          
             u.U6[0][j] = u.U6[xsize-6][j]  ;          
             u.U7[0][j] = u.U7[xsize-6][j]  ;          
             u.U8[0][j] = u.U8[xsize-6][j]  ;                         
              
                    

             q.p[1][j] = q.p[xsize-5][j]  ;         
             q.h[1][j] = q.h[xsize-5][j]  ;          
             q.a[1][j] = q.a[xsize-5][j]  ;          
             q.vx[1][j] =q.vx[xsize-5][j]  ;          
             q.rho[1][j] = q.rho[xsize-5][j]  ;          
             q.vy[1][j] = q.vy[xsize-5][j]  ;          
             q.vz[1][j] = q.vz[xsize-5][j]  ;          
             q.bx[1][j] = q.bx[xsize-5][j]  ;          
             q.by[1][j] = q.by[xsize-5][j]  ;          
             q.bz[1][j] = q.bz[xsize-5][j]  ;          
            
 
            u.U1[1][j] = u.U1[xsize-5][j]  ;          
            u.U2[1][j] = u.U2[xsize-5][j]  ;          
            u.U3[1][j] = u.U3[xsize-5][j]  ;          
            u.U4[1][j] = u.U4[xsize-5][j]  ;          
            u.U5[1][j] = u.U5[xsize-5][j]  ;          
            u.U6[1][j] = u.U6[xsize-5][j]  ;          
            u.U7[1][j] = u.U7[xsize-5][j]  ;          
            u.U8[1][j] = u.U8[xsize-5][j]  ;          
                 
          
             q.p[2][j] = q.p[xsize-4][j]  ;         
             q.h[2][j] = q.h[xsize-4][j]  ;          
             q.a[2][j] = q.a[xsize-4][j]  ;          
                
             q.rho[2][j] = q.rho[xsize-4][j]  ;          
             q.vx[2][j] = q.vx[xsize-4][j]  ;          
             q.vy[2][j] = q.vy[xsize-4][j]  ;          
             q.vz[2][j] = q.vz[xsize-4][j]  ;          
             q.bx[2][j] = q.bx[xsize-4][j]  ;          
             q.by[2][j] = q.by[xsize-4][j]  ;          
             q.bz[2][j] = q.bz[xsize-4][j]  ;          
           
            u.U1[2][j] = u.U1[xsize-4][j]  ;          
            u.U2[2][j] = u.U2[xsize-4][j]  ;          
            u.U3[2][j] = u.U3[xsize-4][j]  ;          
            u.U4[2][j] = u.U4[xsize-4][j]  ;          
            u.U5[2][j] = u.U5[xsize-4][j]  ;          
            u.U6[2][j] = u.U6[xsize-4][j]  ;          
            u.U7[2][j] = u.U7[xsize-4][j]  ;          
            u.U8[2][j] = u.U8[xsize-4][j]  ;          
            
                  
                  
                  }
              /*
             else if (bc2==2)
                {
  
                    q[0][j].rho = q.rho[1][j] ;           
           
                    q[0][j].p = q.p[1][j]   ;       
          
                    q[0][j].vx = -q.vx[1][j] ;           
          
                    q[0][j].vy = -q.vy[1][j] ;           
          
                    q[0][j].h = q.h[1][j] ;           
          
                    q.M[0][j] = q.M[1][j]   ;         
          
                    q.a[0][j] = q.a[1][j];    

                 }
     */
               if(bc3==0)
                {
  
                 
                          
          
                 
                          
          
                     q.p[i][0] = q.p[i][3];   
                     q.h[i][0]= q.h[i][3] ;   
                     q.a[i][0] = q.a[i][3]  ;          
                     u.U1[i][0] = u.U1[i][3]   ;    
     
                     u.U2[i][0] = u.U2[i][3]   ;    
     
                     u.U3[i][0] = u.U3[i][3]  ;    
       
                     u.U4[i][0] = u.U4[i][3]   ;    
     
                     u.U5[i][0] = u.U5[i][3]   ;    
                     u.U6[i][0] = u.U6[i][3]   ;    
                     u.U7[i][0] = u.U7[i][3]   ;    
                     u.U8[i][0] = u.U8[i][3]   ;
                 
  
               
                     q.p[i][1] = q.p[i][3];   
                     q.h[i][1]= q.h[i][3] ;   
                     q.a[i][1] = q.a[i][3]  ;          
                     u.U1[i][1] = u.U1[i][3]   ;    
     
                     u.U2[i][1] = u.U2[i][3]   ;    
     
                     u.U3[i][1] = u.U3[i][3]  ;    
       
                     u.U4[i][1] = u.U4[i][3]   ;    
     
                     u.U5[i][1] = u.U5[i][3]   ;    
                     u.U6[i][1] = u.U6[i][3]   ;    
                     u.U7[i][1] = u.U7[i][3]   ;    
                     u.U8[i][1] = u.U8[i][3]   ;


                     q.p[i][2] = q.p[i][3];   
                     q.h[i][2]= q.h[i][3] ;   
                     q.a[i][2] = q.a[i][3]  ;          
                     u.U1[i][2] = u.U1[i][3]   ;    
     
                     u.U2[i][2] = u.U2[i][3]   ;    
     
                     u.U3[i][2] = u.U3[i][3]  ;    
       
                     u.U4[i][2] = u.U4[i][3]   ;    
     
                     u.U5[i][2] = u.U5[i][3]   ;    
                     u.U6[i][2] = u.U6[i][3]   ;    
                     u.U7[i][2] = u.U7[i][3]   ;    
                     u.U8[i][2] = u.U8[i][3]   ;
               }
               
             else if(bc3==1)
                {
  
                  q.p[i][0] = q.p[i][xsize-6];           
                  q.h[i][0]= q.h[i][xsize-6] ;           
                  q.a[i][0] = q.a[i][xsize-6]  ;          
                  q.rho[i][0] = q.rho[i][xsize-6]  ;          
                  q.vx[i][0] = q.vx[i][xsize-6]  ;          
                  q.vy[i][0] = q.vy[i][xsize-6]  ;          
                  q.vz[i][0] = q.vz[i][xsize-6]  ;          
                  q.bx[i][0] = q.bx[i][xsize-6]  ;          
                  q.by[i][0] = q.by[i][xsize-6]  ;          
                  q.bz[i][0] = q.bz[i][xsize-6]  ;          
                
                 
                  u.U1[i][0] = u.U1[i][xsize-6]   ;    
                  u.U2[i][0] = u.U2[i][xsize-6]   ;    
                  u.U3[i][0] = u.U3[i][xsize-6]  ;    
                  u.U4[i][0] = u.U4[i][xsize-6]   ;    
                  u.U5[i][0] = u.U5[i][xsize-6]   ;    
                  u.U6[i][0] = u.U6[i][xsize-6]   ;    
                  u.U7[i][0] = u.U7[i][xsize-6]   ;    
                  u.U8[i][0] = u.U8[i][xsize-6]   ;  
                                 


                  q.p[i][1] = q.p[i][xsize-5];           
                  q.h[i][1]= q.h[i][xsize-5] ;           
                  q.a[i][1] = q.a[i][xsize-5]  ;          
                  q.vx[i][1] = q.vx[i][xsize-5]  ;          
                  q.rho[i][1] = q.rho[i][xsize-5]  ;          
                  q.vy[i][1] = q.vy[i][xsize-5]  ;          
                  q.vz[i][1] = q.vz[i][xsize-5]  ;          
                  q.bx[i][1] = q.bx[i][xsize-5]  ;          
                  q.by[i][1] = q.by[i][xsize-5]  ;          
                  q.bz[i][1] = q.bz[i][xsize-5]  ;          
                  
 
                  u.U1[i][1] = u.U1[i][xsize-5]   ;    
                  u.U2[i][1] = u.U2[i][xsize-5]   ;    
                  u.U3[i][1] = u.U3[i][xsize-5]  ;    
                  u.U4[i][1] = u.U4[i][xsize-5]   ;    
                  u.U5[i][1] = u.U5[i][xsize-5]   ;    
                  u.U6[i][1] = u.U6[i][xsize-5]   ;    
                  u.U7[i][1] = u.U7[i][xsize-5]   ;    
                  u.U8[i][1] = u.U8[i][xsize-5]   ;  
                  
                  q.p[i][2] = q.p[i][xsize-4];           
                  q.h[i][2]= q.h[i][xsize-4] ;           
                  q.a[i][2] = q.a[i][xsize-4]  ;          
                  q.rho[i][2] = q.rho[i][xsize-4]  ;          
                  q.vx[i][2] = q.vx[i][xsize-4]  ;          
                  q.vy[i][2] = q.vy[i][xsize-4]  ;          
                  q.vz[i][2] = q.vz[i][xsize-4]  ;          
                  q.bx[i][2] = q.bx[i][xsize-4]  ;          
                  q.by[i][2] = q.by[i][xsize-4]  ;          
                  q.bz[i][2] = q.bz[i][xsize-4]  ;          
                
 
                  u.U1[i][2] = u.U1[i][xsize-4]   ;    
                  u.U2[i][2] = u.U2[i][xsize-4]   ;    
                  u.U3[i][2] = u.U3[i][xsize-4]  ;    
                  u.U4[i][2] = u.U4[i][xsize-4]   ;    
                  u.U5[i][2] = u.U5[i][xsize-4]   ;    
                  u.U6[i][2] = u.U6[i][xsize-4]   ;    
                  u.U7[i][2] = u.U7[i][xsize-4]   ;    
                  u.U8[i][2] = u.U8[i][xsize-4]   ;  
                 
                 
                  
  
                  }
              
           
            }

 
        }
  
  
   }
  



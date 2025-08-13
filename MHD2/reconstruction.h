#ifndef RECONSTRUCTION_H_
#define RECONSTRUCTION_H_


#include <cmath>

#define  cr1         0.3333
#define  cr2         1.16666
#define  cr3         1.8333
#define  cr4          .83333
#define   cr5         .166666
 
#define   gamma1     0.1
#define   gamma2     0.6
#define   gamma3     0.3
 
#define    elipson    1e-06
 
using namespace std;


 
double  SI1,SI2,SI3,SI4,omega1,omega2,omega3,omega_1,omega_2,omega_3,w1,w2,w3,w_1,w_2,w_3,S1,S2,S3,S_1,S_2,S_3;

double VANLEER1(double l1 ,double l2 , double l3 );

double VANLEER2(double l1 ,double l2 , double l3 );


void WENO_3X(double u[][xpoints],double UR[][xpoints],double UL[][xpoints]);
void WENO_3Y(double u[][xpoints],double UR[][xpoints],double UL[][xpoints]);


void WENO_X(double u[][xpoints],double UR[][xpoints],double UL[][xpoints]);
void WENO_Y(double u[][xpoints],double UR[][xpoints],double UL[][xpoints]);

/*
void WENO_X(double u[][xpoints],double UR[][xpoints],double UL[][xpoints])
 {

 	
for(i=2; i<=ysize-2; i++)
        { 
 	
          for(j=2; j<=xsize-2; j++)
            {


  
 SI1 = (1.0833)*(u[i][j-2]-2*u[i][j-1]+u[i][j])*(u[i][j-2]-2*u[i][j-1]+u[i][j]) + 0.25*(u[i][j-2]-4*u[i][j-1]+3*u[i][j])*(u[i][j-2]-4*u[i][j-1]+3*u[i][j]);
 
 SI2 = (1.08333)*(u[i][j-1]-2*u[i][j]+u[i][j+1])*(u[i][j-1]-2*u[i][j]+u[i][j+1]) + 0.25*(u[i][j-1]-u[i][j+1])*(u[i][j-1]-u[i][j+1]);
 
 SI3 = (1.08333)*(u[i][j]-2*u[i][j+1]+u[i][j+2])*(u[i][j]-2*u[i][j+1]+u[i][j+2]) + 0.25*(3*u[i][j]-4*u[i][j+1]+u[i][j+2])*(3*u[i][j]-4*u[i][j+1]+u[i][j+2]);
 
 
 omega1= gamma1/((elipson + SI1)*(elipson + SI1));
 omega2= gamma2/((elipson + SI2)*(elipson + SI2));
 omega3= gamma3/((elipson + SI3)*(elipson + SI3));
 
 
 
 omega_1= gamma3/((elipson+SI1)*(elipson+SI1));
 omega_2= gamma2/((elipson+SI2)*(elipson +SI2));
 omega_3= gamma1/((elipson+SI3)*(elipson+SI3));
 
 
 w1=omega1/(omega1+omega2+omega3);
 w2=omega2/(omega1+omega2+omega3);
 w3=omega3/(omega1+omega2+omega3);
 
 w_1=omega_1/(omega_1+omega_2+omega_3);
 w_2=omega_2/(omega_1+omega_2+omega_3);
 w_3=omega_3/(omega_1+omega_2+omega_3);
 
 
 
 S1 = cr1*u[i][j-2]- cr2*u[i][j-1] + cr3*u[i][j]; 
 	
 S2 = - cr5*u[i][j-1] + cr4*u[i][j] + cr1*u[i][j+1]; 
 
 S3 = - cr5*u[i][j+2] + cr4*u[i][j+1] + cr1*u[i][j]; 
 	
 	
 S_3 = cr1*u[i][j+2]- cr2*u[i][j+1] + cr3*u[i][j];
 
 S_2 = - cr5*u[i][j+1] + cr4*u[i][j] + cr1*u[i][j-1]; 
 
 S_1 = - cr5*u[i][j-2] + cr4*u[i][j-1] + cr1*u[i][j]; 
 
 
  UL[i][j]=(S1*w1)+(S2*w2)+(S3*w3);	
 	
  UR[i][j]=(S_1*w_1)+(S_2*w_2)+(S_3*w_3);	
 	
 
 }

}


}

//..................................................................... Y-AXIS.......................................................// 


 //-------------------------------------------------FOLLOWING JIANG AND SHU METHOD.--------------------------------------------------//

//------------------------------------------------------------------------------------------------------------------------------------//

void WENO_Y(double u[][xpoints],double UR[][xpoints],double UL[][xpoints])
{

 	
for(i=2; i<=ysize-2; i++)
        { 
 	
          for(j=2; j<=xsize-2; j++)
            {


  SI1 = (1.08)*(u[i-2][j]-2*u[i-1][j]+u[i][j])*(u[i-2][j]-2*u[i-1][j]+u[i][j]) + 0.25*(u[i-2][j]-4*u[i-1][j]+3*u[i][j])*(u[i-2][j]-4*u[i-1][j]+3*u[i][j]);
 
 SI2 = (1.08)*(u[i-1][j]-2*u[i][j]+u[i+1][j])*(u[i-1][j]-2*u[i][j]+u[i+1][j]) + 0.25*(u[i-1][j]-u[i+1][j])*(u[i-1][j]-u[i+1][j]);
 
 SI3 = (1.08)*(u[i][j]-2*u[i+1][j]+u[i+2][j])*(u[i][j]-2*u[i+1][j]+u[i+2][j]) + 0.25*(3*u[i][j]-4*u[i+1][j]+u[i+2][j])*(3*u[i][j]-4*u[i+1][j]+u[i+2][j]);
 
 
 omega1= gamma1/((elipson + SI1)*(elipson + SI1));
 omega2= gamma2/((elipson + SI2)*(elipson + SI2));
 omega3= gamma3/((elipson + SI3)*(elipson + SI3));
 
 
 
 omega_1= gamma3/((elipson+SI1)*(elipson+SI1));
 omega_2= gamma2/((elipson+SI2)*(elipson +SI2));
 omega_3= gamma1/((elipson+SI3)*(elipson+SI3));
 
 
 w1=omega1/(omega1+omega2+omega3);
 w2=omega2/(omega1+omega2+omega3);
 w3=omega3/(omega1+omega2+omega3);
 
 w_1=omega_1/(omega_1+omega_2+omega_3);
 w_2=omega_2/(omega_1+omega_2+omega_3);
 w_3=omega_3/(omega_1+omega_2+omega_3);
 
 
 
 S1 = cr1*u[i-2][j]- cr2*u[i-1][j] + cr3*u[i][j]; 
 	
 S2 = - cr5*u[i-1][j] + cr4*u[i][j] + cr1*u[i+1][j]; 
 
 S3 = - cr5*u[i+2][j] + cr4*u[i+1][j] + cr1*u[i][j]; 
 	
 	
 S_3 = cr1*u[i+2][j]- cr2*u[i+1][j] + cr3*u[i][j];
 
 S_2 = - cr5*u[i+1][j] + cr4*u[i][j] + cr1*u[i-1][j]; 
 
 S_1 = - cr5*u[i-2][j] + cr4*u[i-1][j] + cr1*u[i][j]; 
 
  
  
  UL[i][j]=(S1*w1)+(S2*w2)+(S3*w3);	
 	
  UR[i][j]=(S_1*w_1)+(S_2*w_2)+(S_3*w_3);	
  
 	 }

   }

return;
}


*/

void WENO_X(double u[][xpoints],double UR[][xpoints],double UL[][xpoints])
{

 	
for(i=2; i<=ysize-2; i++)
        { 
 	
          for(j=2; j<=xsize-2; j++)
            {



   UL[i][j]  =    u[i][j] + (u[i][j+1]-u[i][j])*0.5*VANLEER1(u[i][j],u[i][j+1],u[i][j-1]);	
 	
   UR[i][j]  =    u[i][j] - (u[i][j+1]-u[i][j])*0.5*VANLEER1(u[i][j],u[i][j+1],u[i][j-1]);

            }
   
        }
    
    }

void WENO_Y(double u[][xpoints],double UR[][xpoints],double UL[][xpoints])
{

 	
for(i=2; i<=ysize-2; i++)
        { 
 	
          for(j=2; j<=xsize-2; j++)
            {

   UL[i][j]  =    u[i][j] + (u[i+1][j]-u[i][j])*0.5*VANLEER1(u[i][j],u[i+1][j],u[i-1][j]);	
 	
   UR[i][j]  =    u[i][j] - (u[i+1][j]-u[i][j])*0.5*VANLEER1(u[i][j],u[i+1][j],u[i-1][j]);


 	 
	     }
       
          }

   }

double VANLEER1(double l1 ,double l2 , double l3 )

{

    double phi,r;

   r= (l1-l3)/(l2-l1);

    if(r<=0) 
     phi=0;
 
    else if(r>0)
    phi = (r+abs(r))/(1+abs(r));
 
   

    return phi;

         }

double VANLEER2(double l1 ,double l2 , double l3 )

{

 double phi,r;

    r= (l1-l3)/(l2-l1);
  
  /*
   if(r<=0) 
    phi=0;

   else if(r > 0 && r < 1)
   phi=r;
 
   else 
   phi=1;
*/

    phi = std::max(0.0, std::min(1.0, r));

   return phi;

}


void WENO_3X(double u[][xpoints],double UR[][xpoints],double UL[][xpoints])
 {

 	
for(i=2; i<=ysize-2; i++)
        { 
 	
          for(j=2; j<=xsize-2; j++)
            {


  
   SI1 =  (u[i][j]-u[i][j-1])*(u[i][j]-u[i][j-1]);
 
   SI2 =  (u[i][j+1]-u[i][j])*(u[i][j+1]-u[i][j]);
  
  omega1= (0.333)/((elipson + SI1)*(elipson + SI1));
  omega2= (0.666)/((elipson + SI2)*(elipson + SI2));
 
 
 
  omega_1= (0.666)/((elipson+SI1)*(elipson+ SI1));
  omega_2= (0.333)/((elipson+SI2)*(elipson +SI2));
 
 
  w1=omega1/(omega1+omega2);
  w2=omega2/(omega1+omega2);
 
  w_1=omega_1/(omega_1+omega_2);
  w_2=omega_2/(omega_1+omega_2);
  

 S1= -0.5*u[i][j-1] + 1.5*u[i][j];

 S2=  0.5*u[i][j] + 0.5*u[i][j+1];

 S_1 =  0.5*u[i][j] + 0.5*u[i][j-1];

 S_2 =  -0.5*u[i][j+1] + 1.5*u[i][j];

 
  UL[i][j]=(S1*w1)+(S2*w2);	
 	
  UR[i][j]=(S_1*w_1)+(S_2*w_2);	
  
 
       }

    }
}

//..................................................................... Y-AXIS.......................................................// 


 //-------------------------------------------------FOLLOWING JIANG AND SHU METHOD.--------------------------------------------------//

//------------------------------------------------------------------------------------------------------------------------------------//

void WENO_3Y(double u[][xpoints],double UR[][xpoints],double UL[][xpoints])
{

 	
for(i=2; i<=ysize-2; i++)
        { 
 	
          for(j=2; j<=xsize-2; j++)
            {



   SI1 =  (u[i][j]-u[i-1][j])*(u[i][j]-u[i-1][j]);
 
   SI2 =  (u[i+1][j]-u[i][j])*(u[i+1][j]-u[i][j]);
  
   omega1= (0.333)/((elipson + SI1)*(elipson + SI1));
   omega2= (0.666)/((elipson + SI2)*(elipson + SI2));
 
 
 
  omega_1= (0.666)/((elipson+SI1)*(elipson+ SI1));
  omega_2= (0.333)/((elipson+SI2)*(elipson +SI2));


  w1=omega1/(omega1+omega2);
  w2=omega2/(omega1+omega2);
 
  w_1=omega_1/(omega_1+omega_2);
  w_2=omega_2/(omega_1+omega_2);
  

 S1= -0.5*u[i-1][j] + 1.5*u[i][j];

 S2=  0.5*u[i][j] + 0.5*u[i+1][j];

 S_1 =  0.5*u[i][j] + 0.5*u[i-1][j];

 S_2 =  -0.5*u[i+1][j] + 1.5*u[i][j];

 
  UL[i][j]=(S1*w1)+(S2*w2);	
 	
  UR[i][j]=(S_1*w_1)+(S_2*w_2);	
  
  
 	 }

   }

}

 #endif


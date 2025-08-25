#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fstream>
#include <iostream>
#include "mesh.h"
#include <mpi.h>
#include <time.h>
#include <sstream>
#include "initiate.h"
//#include "deallocate.h"


#define  DOWN                101
#define  UP                  100
#define  EIGEN               200
#define  SENDY               210

#define  q1              120
#define  q2              121
#define  q3              123
#define  q4              124
#define  q5              125
#define  q6              656


#define  S1              127
#define  S2              128
#define  S3              129
#define  S4              130

#define  S8              548
#define  S5              549
#define  S6              550
#define  S7              551


#define  q7              131
#define  q8              132
#define  q9              133
#define  q10             134
#define  q11             553
#define  q12             554
#define  q13             555
#define  q14             556


#define  R               135
#define  P               136
#define  VY              137
#define  VX              138
#define  Vz              139
#define  BY              140
#define  tb               141

#define  BC2             4002
#define  BC3             4003
#define  BC4             4004
#define  BC5             4005
#define  BC6             4006
#define  BC7             4007
#define  BC8             4008
#define  BC9             4011
#define  BC91             4025
#define  BC92             4026
#define  BC93             4027

#define  PC2             4012
#define  PC3             4013
#define  PC4             4014
#define  PC5             4015
#define  PC6             4016
#define  PC7             4017
#define  PC8             4018
#define  PC9             4019
#define  PC91             4020
#define  PC92             4021
#define  PC93             4022

#define  PC1               4001

#define  ROOT                 0
#define  Y23             40225

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
void grid();
void wave_system();


int offset2(int prcs);
int offset1(int prcs);
int offset3(int prcs);
int offset4(int prcs);

double dt=1e-3 , dx, dy , t=0.0 , CFL , tf ,  bc1 ,bc2 ,bc3 , bc4 ;

double x[xpoints], t1, t2;

double *y;     

int prcs,np, width , st, stp, ymax, iter, PartitionSize;

 MPI_Status status;

 using namespace std;  

   int main()
  {


    dx = L/xcells;
    dy = B/ycells;
    ymax=ypoints;
   
   
 
   MPI_Init(NULL, NULL);
  
   MPI_Comm_size(MPI_COMM_WORLD, &np);

   MPI_Comm_rank(MPI_COMM_WORLD, &prcs);

 
   t1=clock();

   if(prcs==0){
   input();
 
   width=ycells/np;
 
   st=1;

    PartitionSize = width+3;
    stp =  PartitionSize-2;
   

   for(i=1; i<=np-1; i++){
  
   MPI_Send( &width,  1, MPI_INT, prcs+i, DOWN,   MPI_COMM_WORLD);
   


          }

    }

  //grid();
  MPI_Bcast(&tf, 1, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);
  MPI_Bcast(&CFL, 1, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);
  MPI_Bcast(&bc1, 1, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);
  MPI_Bcast(&bc2, 1, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);
  MPI_Bcast(&bc3, 1, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);
  MPI_Bcast(&bc4, 1, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);
 
  
if(prcs!=0){


   MPI_Recv( &width, 1 , MPI_INT, ROOT, DOWN  ,  MPI_COMM_WORLD, &status);

   if(prcs!=np-1) 
    { 
                                    // will be in a loop separately ./
   st=1;// width*prcs + 1;

   //stp=width*(prcs+1);
  
   PartitionSize = width + 2;
   
   stp =  PartitionSize-2; 


       }

  else if (prcs==np-1){

   st=1;//width*prcs + 1;

 //  stp=width*(prcs+1)+ (ymax%np-1);
 
    PartitionSize = width + 2;
    stp =  PartitionSize-2;

     }

  }


      
     distribution(PartitionSize , st , stp);

     grid();   

     initilization();

 //cout<<"WIDTH\t"<<y[st]<<"\t"<<y[stp]<<"\tPRCS\t"<<prcs<<"\tPS\t"<<PartitionSize<<endl;
     boundaryConditions();
         
     cout<<"WIDTH\t"<<y[st]<<"\t"<<y[stp]<<"\tPRCS\t"<<prcs<<"\tPS\t"<<PartitionSize<<endl;


     iter=0;
     
           t1=clock();
	   do{

             RK3();
	    
          if(prcs==0)
            time_step();
     
           MPI_Bcast(&dt, 1, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);
           MPI_Bcast(&t, 1, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);  
            iter++;
       
	   if(prcs==0)
             tracer();
 	   
          }	   while(t<=tf);
 

   t2 = clock();   	 
 
  MPI_Barrier( MPI_COMM_WORLD);
  

     

    MPI_Send( y, PartitionSize  , MPI_DOUBLE, 0, Y23,   MPI_COMM_WORLD);

  
  
     for(i=st; i<=stp; i++)
    {

	     

     MPI_Send( q.rho[i], xpoints  , MPI_DOUBLE, 0, R,   MPI_COMM_WORLD);
     MPI_Send( q.p[i], xpoints  , MPI_DOUBLE, 0, P,   MPI_COMM_WORLD);
     MPI_Send( q.vx[i], xpoints  , MPI_DOUBLE, 0, VX,   MPI_COMM_WORLD);
     MPI_Send( q.vy[i], xpoints  , MPI_DOUBLE, 0, VY,   MPI_COMM_WORLD);
     MPI_Send( q.by[i], xpoints  , MPI_DOUBLE, 0, BY,   MPI_COMM_WORLD);
     MPI_Send( q.b[i], xpoints  , MPI_DOUBLE, 0, tb,   MPI_COMM_WORLD);
  
    }


 
  if(prcs==0){
     cout<<"WRITING FILES....\n"<<endl;
     output();
     cout<<"THE END OF THE PROGRAMME IN =\t"<<(t2-t1)/CLOCKS_PER_SEC<<"\tsecs"<<endl;
	  
    }

    MPI_Finalize(); 
 return 0;

 }
   


    

     void initilization()
     {

      
    

     double r,f,r_0=0.1, r_1=0.115;; 
    for(i=st; i<=stp; i++)
    {
 	
    for(j=0; j<=xsize; j++)
          
        {
	   q.p[i][j] = 1/GAMMA;
	
	    q.rho[i][j] = 1; 	         
    
	    q.by[i][j] = sin(4*3.142*x[j])/GAMMA; 
  	
	    q.vx[i][j] = -sin(2*3.142*y[i]);
      
	    q.vy[i][j] = sin(2*3.142*x[j]);
      
	  	q.vz[i][j] = 0;
	  
	  	q.bz[i][j]         = 0 ;
	  
	    q.bx[i][j]        = -sin(2*3.142*y[i])/GAMMA; 
   

        }

      }
   
 // need to be checked ..//	  
	 

     for(i=st; i<=stp; i++)
   {
 	
    for(j=0; j<=xsize; j++)
     {
          
    

           q.a[i][j] = sqrt(GAMMA*q.p[i][j]/q.rho[i][j]) ;
	 
	   q.h[i][j] = 0.5*(q.vx[i][j]*q.vx[i][j]) + 0.5*(q.vy[i][j]*q.vy[i][j])+ (q.a[i][j]*q.a[i][j])/(GAMMA-1)+
	  
               ((q.bz[i][j]*q.bz[i][j] + q.by[i][j]*q.by[i][j]+q.bx[i][j]*q.bx[i][j])/q.rho[i][j]);
         
           u.U1[i][j] = q.rho[i][j];
      
	   u.U2[i][j] = q.rho[i][j]*q.vx[i][j];
      
	   u.U3[i][j] = q.rho[i][j]*q.vy[i][j];

           u.U4[i][j] = q.rho[i][j]*q.vz[i][j];
       
           u.U5[i][j] = q.p[i][j]/(GAMMA-1) + 0.5*q.rho[i][j]*(q.vy[i][j]*q.vy[i][j]+ q.vx[i][j]*q.vx[i][j])+0.5*(q.bz[i][j]*q.bz[i][j] + q.by[i][j]*q.by[i][j]+q.bx[i][j]*q.bx[i][j]);

           u.U6[i][j] = q.by[i][j];

           u.U7[i][j] = q.bx[i][j];

           u.U8[i][j] = q.bz[i][j];


       }
 
    }

    //  cout<<"yes"<<endl;
//boundaryConditions();  // need to be checked ...//
 //cout<<"yes"<<endl;
  
  }

void wave_system()
{

double temp1,alfven;

int start , stop;

if(prcs==0)
  start = 0;

else 
   start = st;

	
for(i=start; i<=stp; i++)
 {
 	
    for(j=0; j<=xsize; j++)
     {
 

   w.afy[i][j]  =  u.U6[i][j]/sqrt(u.U1[i][j]);

   w.afx[i][j]  =  u.U7[i][j]/sqrt(u.U1[i][j]);

   alfven     = sqrt((u.U7[i][j]*u.U7[i][j])+(u.U6[i][j]*u.U6[i][j])+u.U8[i][j]*u.U8[i][j])/sqrt(u.U1[i][j]);

   temp1 = (q.a[i][j]*q.a[i][j])+(alfven*alfven);
  
    
   w.fwx[i][j] = sqrt(0.5*((temp1)+sqrt((temp1*temp1)-(4*q.a[i][j]*q.a[i][j]*w.afx[i][j]*w.afx[i][j]))));
	   
  	

   w.fwy[i][j] = sqrt(0.5*((temp1)+sqrt((temp1*temp1)-(4*q.a[i][j]*q.a[i][j]*w.afy[i][j]*w.afy[i][j]))));
	

 	

     
     }
  
   }
     
//cout<<"yes\t"<<prcs<<endl;
 

}



void eigen_value()
{

  int start , stop;

if(prcs==0)
  start = 0;

else
   start = st;
   
 for(i=start; i<=stp; i++)
 {
 	
    for(j=0; j<=xsize; j++)
     {
 

  

     c.lx1[i][j]  =  u.U2[i][j]/u.U1[i][j] + w.fwx[i][j];

     c.lx2[i][j]  = u.U2[i][j]/u.U1[i][j] - w.fwx[i][j];

     c.ly1[i][j]  = u.U3[i][j]/u.U1[i][j] + w.fwy[i][j];

     c.ly2[i][j] =  u.U3[i][j]/u.U1[i][j] - w.fwy[i][j];
   
 


     }
  
   }
     
 }
  
 void maxEigen()
   {


   lamdaxMAX=0;
   lamdayMAX=0;
  
      

    for(i=st; i<=stp; i++)
   { 
 	
    for(j=0; j<=xsize; j++)
        {


      if(lamdaxMAX<fabs(c.lx1[i][j]))
       lamdaxMAX=fabs(c.lx1[i][j]);

      if(lamdayMAX<fabs(c.ly1[i][j]))
        lamdayMAX=fabs(c.ly1[i][j]);
     
 
       }
  

    }

   
   for(i=st; i<=stp; i++)
     { 
 	
    for(j=0; j<=xsize; j++)
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



if( prcs != 0 )    
         MPI_Send( &lamdaMAX,  1, MPI_DOUBLE, 0, EIGEN,   MPI_COMM_WORLD);  




   }
  
	

  void Rusanov()
  {

 int start , stop;

if(prcs==0)
  start = 0;

else
   start = st;


       for(i=st; i<=stp; i++)
        { 
 	
          for(j=0; j<=xsize; j++)
            {


               if(fabs(c.lx2[i][j])>fabs(c.lx1[i][j]))
                 c.lMaxx[i][j] = fabs(c.lx2[i][j]);

               else 
                 c.lMaxx[i][j] = fabs(c.lx1[i][j]);

            
 
            }   
  

        }

     

      for(i=start; i<=stp; i++)
        { 
 	
          for(j=0; j<=xsize; j++)
            {


               if(fabs(c.ly2[i][j])>fabs(c.ly1[i][j]))
                 c.lMaxy[i][j] = fabs(c.ly2[i][j]);

               else 
                 c.lMaxy [i][j]= fabs(c.ly1[i][j]);
  
            
 
             }   
  

         }


    if( prcs != 0 )    
         MPI_Send( c.lMaxy[st],  xpoints, MPI_DOUBLE, prcs-1, SENDY,   MPI_COMM_WORLD);  




   
  if( prcs != np-1 )       
        MPI_Recv( c.lMaxy[stp+1], xpoints, MPI_DOUBLE, prcs+1,SENDY, MPI_COMM_WORLD, &status);
     

      for(i=st; i<=stp; i++)
        { 
 	
          for(j=0; j<=xsize-1; j++)
            {


               if(c.lMaxx[i][j+1] > c.lMaxx[i][j] )
                 c.lx[i][j] = c.lMaxx[i][j+1];

               else
                 c.lx[i][j] = c.lMaxx[i][j];
     
 
             }   
  

         }

      for(i=start; i<=stp; i++)
        { 
 	
          for(j=0; j<=xsize-1; j++)
            {


               if(c.lMaxy[i+1][j] > c.lMaxy[i][j] )
                 c.ly[i][j] = c.lMaxy[i+1][j];

               else
                 c.ly [i][j]= c.lMaxy[i][j];
     
 
             }   
  

         }


     
 if( prcs != 0 ){     
                         MPI_Send( u.U1[st],  xpoints, MPI_DOUBLE, prcs-1, S1,   MPI_COMM_WORLD);  
                         MPI_Send( u.U2[st],  xpoints, MPI_DOUBLE, prcs-1, S2,   MPI_COMM_WORLD);   
                         MPI_Send( u.U3[st],  xpoints, MPI_DOUBLE, prcs-1, S3,   MPI_COMM_WORLD);   
                         MPI_Send( u.U4[st],  xpoints, MPI_DOUBLE, prcs-1, S4,   MPI_COMM_WORLD);   
                        
                         MPI_Send( u.U5[st],  xpoints, MPI_DOUBLE, prcs-1, S5,   MPI_COMM_WORLD);  
                         MPI_Send( u.U6[st],  xpoints, MPI_DOUBLE, prcs-1, S6,   MPI_COMM_WORLD);  
                         MPI_Send( u.U7[st],  xpoints, MPI_DOUBLE, prcs-1, S7,   MPI_COMM_WORLD);  
                         MPI_Send( u.U8[st],  xpoints, MPI_DOUBLE, prcs-1, S8,   MPI_COMM_WORLD);  
      
        
 
 
              }
      

  
   
    if( prcs != np-1 ){       
                      MPI_Recv( u.U1[stp+1], xpoints, MPI_DOUBLE, prcs+1,S1, MPI_COMM_WORLD, &status);
                      MPI_Recv( u.U2[stp+1], xpoints, MPI_DOUBLE, prcs+1,S2, MPI_COMM_WORLD, &status);
                      MPI_Recv( u.U3[stp+1], xpoints, MPI_DOUBLE, prcs+1,S3, MPI_COMM_WORLD, &status);
                      MPI_Recv( u.U4[stp+1], xpoints, MPI_DOUBLE, prcs+1,S4, MPI_COMM_WORLD, &status);
                     
                      MPI_Recv( u.U5[stp+1], xpoints, MPI_DOUBLE, prcs+1,S5, MPI_COMM_WORLD, &status);   
                      MPI_Recv( u.U6[stp+1], xpoints, MPI_DOUBLE, prcs+1,S6, MPI_COMM_WORLD, &status);   
                      MPI_Recv( u.U7[stp+1], xpoints, MPI_DOUBLE, prcs+1,S7, MPI_COMM_WORLD, &status);   
                      MPI_Recv( u.U8[stp+1], xpoints, MPI_DOUBLE, prcs+1,S8, MPI_COMM_WORLD, &status);   
               
 

 
                     }


    

 if( prcs != 0 ){     
         
                         MPI_Send( q.h[st],  xpoints, MPI_DOUBLE, prcs-1, q4,   MPI_COMM_WORLD);  
                         MPI_Send( q.p[st],  xpoints, MPI_DOUBLE, prcs-1, q5,   MPI_COMM_WORLD);                                      
                        
                     

 
              }

  
   
    if( prcs != np-1 ){       
          
                      MPI_Recv( q.h[stp+1], xpoints, MPI_DOUBLE,prcs+1,q4, MPI_COMM_WORLD, &status);    
                      MPI_Recv( q.p[stp+1], xpoints, MPI_DOUBLE,prcs+1,q5, MPI_COMM_WORLD, &status);   
                   
                     }
      	  

  for(i=start; i<=stp; i++)
        { 
 	
          for(j=0; j<=xsize-1; j++)
            {


    f.F1[i][j] = ((u.U2[i][j]+ u.U2[i][j+1]) - (u.U1[i][j+1]-u.U1[i][j])*c.lx[i][j])*0.5;
	

    f.F2[i][j] = 0.5*(u.U2[i][j]*u.U2[i][j]/u.U1[i][j]+ u.U2[i][j+1]*u.U2[i][j+1]/u.U1[i][j+1])  + 0.5*(q.p[i][j]+q.p[i][j+1])+ 

0.25*(((u.U6[i][j]*u.U6[i][j])+(u.U7[i][j]*u.U7[i][j]))+(u.U6[i][j+1]*u.U6[i][j+1])+(u.U7[i][j+1]*u.U7[i][j+1])+(u.U8[i][j]*u.U8[i][j])+(u.U8[i][j+1]*u.U8[i][j+1]))

    -((u.U7[i][j]*u.U7[i][j])+(u.U7[i][j+1]*u.U7[i][j+1]))/2 - (u.U2[i][j+1]-u.U2[i][j])*c.lx[i][j]*0.5;
                  

		   f.F3[i][j] = 0.5*(u.U2[i][j]*u.U3[i][j]/u.U1[i][j] + u.U2[i][j+1]*u.U3[i][j+1]/u.U1[i][j+1])  -((u.U6[i][j]*u.U7[i][j])+u.U6[i][j+1]*u.U7[i][j+1])*0.5 - (u.U3[i][j+1]-u.U3[i][j])*c.lx[i][j]*0.5; 

		
       		   f.F4[i][j] = 0.5*(u.U2[i][j]*u.U4[i][j]/u.U1[i][j]+ u.U2[i][j+1]*u.U4[i][j+1]/u.U1[i][j+1])  -((u.U8[i][j]*u.U7[i][j])+(u.U8[i][j+1]*u.U7[i][j+1]))*0.5- (u.U4[i][j+1]-u.U4[i][j])*c.lx[i][j]*0.5; 


   f.F5[i][j] = 0.5*(u.U2[i][j]*q.h[i][j]+ u.U2[i][j+1]*q.h[i][j+1])- (u.U7[i][j]*(u.U7[i][j]*u.U2[i][j]/u.U1[i][j]+u.U6[i][j]*u.U3[i][j]/u.U1[i][j]+ u.U8[i][j]*u.U4[i][j]/u.U1[i][j])+ (u.U7[i][j+1]*(u.U7[i][j+1]*u.U2[i][j+1]/u.U1[i][j+1]+u.U6[i][j+1]*u.U3[i][j+1]/u.U1[i][j+1]+ u.U8[i][j+1]*u.U4[i][j+1]/u.U1[i][j+1])))*0.5 - (u.U5[i][j+1]-u.U5[i][j])*c.lx[i][j]*0.5;
	



                   f.F6[i][j]=((u.U6[i][j]*u.U2[i][j]/u.U1[i][j])+(u.U6[i][j+1]*u.U2[i][j+1]/u.U1[i][j+1]))*.5 - ((u.U7[i][j]*u.U3[i][j]/u.U1[i][j])+(u.U7[i][j+1]*u.U3[i][j+1]/u.U1[i][j+1]))*.5-(u.U6[i][j+1]-u.U6[i][j])*c.lx[i][j]*0.5; 

                  f.F7[i][j]= -(u.U7[i][j+1]-u.U7[i][j])*c.lx[i][j]*0.5; //bx
		   
                  f.F8[i][j]= ((u.U8[i][j]*u.U2[i][j]/u.U1[i][j])+(u.U8[i][j+1]*u.U2[i][j+1]/u.U1[i][j+1]))*.5 - ((u.U7[i][j]*u.U4[i][j]/u.U1[i][j])+(u.U7[i][j+1]*u.U4[i][j+1]/u.U1[i][j+1]))*.5-(u.U8[i][j+1]-u.U8[i][j])*c.lx[i][j]*0.5;


	// checked.......done...	/////////.............................................................../ for y direction ...////////

                   f.G1[i][j] = ((u.U1[i][j]*u.U3[i][j]/u.U1[i][j] + u.U1[i+1][j]*u.U3[i+1][j]/u.U1[i+1][j]) - (u.U1[i+1][j]-u.U1[i][j])*c.ly[i][j])*0.5;       
		  
    
 f.G2[i][j] = 0.5*(u.U3[i][j]*u.U2[i][j]/u.U1[i][j]+ u.U3[i+1][j]*u.U2[i+1][j]/u.U1[i+1][j]) -((u.U6[i][j]*u.U7[i][j])+u.U6[i+1]
[j]*u.U7[i+1][j])*0.5 - (u.U2[i+1][j]-u.U2[i][j])*c.ly[i][j]*0.5; 

    


	   f.G3[i][j] = 0.5*(u.U3[i][j]*u.U3[i][j]/u.U1[i][j]+ u.U3[i+1][j]*u.U3[i+1][j]/u.U1[i+1][j]+q.p[i][j]+q.p[i+1][j]) - (u.U3[i+1][j]-u.U3[i][j])*c.ly[i][j]*0.5 +0.25*(((u.U6[i][j]*u.U6[i][j])+(u.U7[i][j]*u.U7[i][j])+(u.U8[i][j]*u.U8[i][j]))+(u.U6[i+1][j]*u.U6[i+1][j])+(u.U7[i+1][j]*u.U7[i+1][j])+(u.U8[i+1][j]*u.U8[i+1][j])) 
                  - ((u.U6[i+1][j]*u.U6[i+1][j])+(u.U6[i][j]*u.U6[i][j]))/2;
 
		   
	   f.G4[i][j] = 0.5*(u.U3[i][j]*u.U4[i][j]/u.U1[i][j]+ u.U3[i+1][j]*u.U4[i+1][j]/u.U1[i+1][j]) -((u.U6[i][j]*u.U8[i][j])+u.U6[i+1][j]*u.U8[i+1][j])*0.5 - (u.U4[i+1][j]-u.U4[i][j])*c.ly[i][j]*0.5; 

 

	   f.G5[i][j] = 0.5*(u.U3[i][j]*q.h[i][j]+ u.U3[i+1][j]*q.h[i+1][j]) - (u.U5[i+1][j]-u.U5[i][j])*c.ly[i][j]*0.5 - (.5*(u.U6[i][j]*(u.U7[i][j]*u.U2[i][j]/u.U1[i][j]+(u.U6[i][j]*u.U3[i][j]/u.U1[i][j])+(u.U8[i][j]*u.U4[i][j]/u.U1[i][j]))+(u.U6[i+1][j]*(u.U7[i+1][j]*u.U2[i+1][j]/u.U1[i+1][j]+(u.U6[i+1][j]*u.U3[i+1][j]/u.U1[i+1][j])+(u.U8[i+1][j]*u.U4[i+1][j]/u.U1[i+1][j]))))) ;
	   

      	   f.G6[i][j]=-(u.U6[i+1][j]-u.U6[i][j])*c.ly[i][j]*0.5; //by
 
           f.G7[i][j]= ((u.U7[i][j]*u.U3[i][j]/u.U1[i][j])+(u.U7[i+1][j]*u.U3[i+1][j]/u.U1[i+1][j]))*.5 - ((u.U6[i][j]*u.U2[i][j]/u.U1[i][j])+(u.U6[i+1][j]*u.U2[i+1][j]/u.U1[i+1][j]))*.5- (u.U7[i+1][j]-u.U7[i][j])*c.ly[i][j]*0.5; 

           f.G8[i][j]= ((u.U8[i][i]*u.U3[i][j]/u.U1[i][j])+(u.U8[i+1][j]*u.U3[i+1][j]/u.U1[i+1][j]))*.5 - ((u.U6[i][j]*u.U4[i][j]/u.U1[i][j])+(u.U6[i+1][j]*u.U4[i+1][j]/u.U1[i+1][j]))*.5-(u.U8[i+1][j]-u.U8[i][j])*c.ly[i][j]*0.5;
    
//---------------------------------------------------------------------------------------// end ..//////////////////
           }

         }

        //......................................................................................................//  
   
     
double NetdivBx, NetdivBy;

//start= offset1(prcs);
//stop = offset2(prcs);

for(i=st; i<=stp; i++)
 {
 	
    for(j=1; j<=xsize-1; j++)
     {
 

  

    

                NetdivBx=   (u.U7[i][j+1]-u.U7[i][j]); 
   
                NetdivBy =  (u.U6[i+1][j]-u.U6[i][j]);


                S.s1[i][j]=0;
        
		S.s2[i][j] = (NetdivBx+ NetdivBy)*u.U7[i][j]; 
		
		S.s3[i][j] = (NetdivBx+ NetdivBy)*u.U6[i][j]; 
		 
		
		S.s4[i][j] = (NetdivBx+ NetdivBy)*u.U8[i][j]; 
		 
		
		S.s5[i][j] = ((u.U7[i][j]*q.vx[i][j])+(u.U6[i][j]*q.vy[i][j])+(u.U8[i][j]*q.vz[i][j]))*(NetdivBy +NetdivBx); 
		
		S.s8[i][j] = q.vz[i][j]*(NetdivBx+ NetdivBy);  
		
		S.s7[i][j] = q.vx[i][j]*(NetdivBx+ NetdivBy);
        
		S.s6[i][j] = q.vy[i][j]*(NetdivBx+ NetdivBy);
  
  
          }

       }
    
   }
    
void RK3()
{


   iteration();

//start= offset1(prcs);
//stop = offset2(prcs);

if( prcs != np-1 ){     
         
                         MPI_Send( f.G1[stp],  xpoints, MPI_DOUBLE, prcs+1, q7,   MPI_COMM_WORLD);  
                         MPI_Send( f.G2[stp],  xpoints, MPI_DOUBLE, prcs+1, q8,   MPI_COMM_WORLD);  
                         MPI_Send( f.G3[stp],  xpoints, MPI_DOUBLE, prcs+1, q9,   MPI_COMM_WORLD);  
                         MPI_Send( f.G4[stp],  xpoints, MPI_DOUBLE, prcs+1, q10,   MPI_COMM_WORLD);  
                         MPI_Send( f.G5[stp],  xpoints, MPI_DOUBLE, prcs+1, q11,   MPI_COMM_WORLD);  
                         MPI_Send( f.G6[stp],  xpoints, MPI_DOUBLE, prcs+1, q12,   MPI_COMM_WORLD);  
                         MPI_Send( f.G7[stp],  xpoints, MPI_DOUBLE, prcs+1, q13,   MPI_COMM_WORLD);  
                         MPI_Send( f.G8[stp],  xpoints, MPI_DOUBLE, prcs+1, q14,   MPI_COMM_WORLD);  
              
                         
                 }
      

  
   
    if( prcs != 0 ){       
                      MPI_Recv( f.G1[st-1], xpoints, MPI_DOUBLE, prcs-1,q7, MPI_COMM_WORLD, &status);
                      MPI_Recv( f.G2[st-1], xpoints, MPI_DOUBLE, prcs-1,q8, MPI_COMM_WORLD, &status);
                      MPI_Recv( f.G3[st-1], xpoints, MPI_DOUBLE, prcs-1,q9, MPI_COMM_WORLD, &status);    
                      MPI_Recv( f.G4[st-1], xpoints, MPI_DOUBLE, prcs-1,q10, MPI_COMM_WORLD, &status);   
                      MPI_Recv( f.G5[st-1], xpoints, MPI_DOUBLE, prcs-1,q11, MPI_COMM_WORLD, &status);   
                      MPI_Recv( f.G6[st-1], xpoints, MPI_DOUBLE, prcs-1,q12, MPI_COMM_WORLD, &status);   
                      MPI_Recv( f.G7[st-1], xpoints, MPI_DOUBLE, prcs-1,q13, MPI_COMM_WORLD, &status);   
                      MPI_Recv( f.G8[st-1], xpoints, MPI_DOUBLE, prcs-1,q14, MPI_COMM_WORLD, &status);   
 
                    


                    }

 
             


  for(i=st; i<=stp; i++)
        { 
 	
          for(j=1; j<=xsize-1; j++)
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
 
   for(i=st; i<=stp; i++)
        { 
 	
          for(j=1; j<=xsize-1; j++)
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
	
	 q.h[i][j] =  0.5*(q.vx[i][j]*q.vx[i][j]) + 0.5*(q.vy[i][j]*q.vy[i][j])+ (q.a[i][j]*q.a[i][j])/(GAMMA-1)+
	  
               ((q.bz[i][j]*q.bz[i][j] + q.by[i][j]*q.by[i][j]+q.bx[i][j]*q.bx[i][j])/q.rho[i][j]);
 
   

           }

  
      }
   

boundaryConditions();
 iteration();
 

  for(i=st; i<=stp; i++)
        { 
 	
          for(j=1; j<=xsize-1; j++)
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
	
	 q.h[i][j] =  0.5*(q.vx[i][j]*q.vx[i][j]) + 0.5*(q.vy[i][j]*q.vy[i][j])+ (q.a[i][j]*q.a[i][j])/(GAMMA-1)+
	  
               ((q.bz[i][j]*q.bz[i][j] + q.by[i][j]*q.by[i][j]+q.bx[i][j]*q.bx[i][j])/q.rho[i][j]);
 
         q.b[i][j] = sqrt(q.bx[i][j]*q.bx[i][j]+ q.by[i][j]*q.by[i][j] + q.bz[i][j]*q.bz[i][j]);
   
     }


  }

boundaryConditions();


  

}



void tracer()
   {

for(i=st; i<=stp; i++)
        { 
 	
          for(j=0; j<=xsize; j++)
            {


                if(i==stp-4 && j==xcells-5) {
                 cout<<"==============================================================================="<<endl;
                 cout<<"EVOLVING FOR ITERATION :: "<<t<<endl;
                 cout<<"DENSITY = "<<q.rho[i][j]<<"\n";
                 cout<<"PRESSURE = "<<q.p[i][j]<<"\n";
                 cout<<"VELOCITY = "<<q.vx[i][j]<<"\n";
                 cout<<"MAGNETIC-Y = "<<q.by[i][j]<<"\n";
                 cout<<"MAGNETIC-x = "<<q.bx[i][j]<<"\n";
                 cout<<"LAMDA_MAX ="<<lamdaxMAX<<"\n";
                 cout<<"DELTA_T ="<<dt<<"\n";
         
                  }    
   
             }

        }

 }

void time_step()
{

double temp=lamdaMAX;

for(i=1;i<=np-1;i++){
  MPI_Recv( &lamdaMAX, 1, MPI_DOUBLE, i,EIGEN, MPI_COMM_WORLD, &status);
     if(lamdaMAX>temp)
      temp=lamdaMAX; 
  
 
    }
	dt= 0.5*CFL*(dx/temp);
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


     ofstream vx,p,vy,rho,by,bx,Bt;
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


 cout<<"WRITING FOR PROCESS-"<<"\t"<<"1"<<endl;
 for(i=st; i<=stp; i++)
        { 
 	
          for(j=1; j<=xsize-1; j++)
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
             by<<endl;
             Bt<<endl; 
  
              
       }



 // PartitionSize +=1;
 // deAllocate(PartitionSize);
 
  stp=stp-1;

   for(int u=1; u<=np-1; u++) {	


   cout<<"WRITING FOR PROCESS-"<<"\t"<<u+1<<endl;

   for(i=st; i<=stp; i++)
  {
     
     MPI_Recv( q.rho[i], xpoints , MPI_DOUBLE, u, R  ,  MPI_COMM_WORLD, &status);
     MPI_Recv( q.p[i], xpoints , MPI_DOUBLE, u, P  ,  MPI_COMM_WORLD, &status);
     MPI_Recv( q.vx[i], xpoints , MPI_DOUBLE, u, VX  ,  MPI_COMM_WORLD, &status);
     MPI_Recv( q.vy[i], xpoints , MPI_DOUBLE, u, VY  ,  MPI_COMM_WORLD, &status);
     MPI_Recv( q.by[i], xpoints , MPI_DOUBLE, u, BY  ,  MPI_COMM_WORLD, &status);
     MPI_Recv( q.b[i], xpoints , MPI_DOUBLE, u,   tb ,  MPI_COMM_WORLD, &status);


            }

     MPI_Recv( y, PartitionSize , MPI_DOUBLE, u,   Y23 ,  MPI_COMM_WORLD, &status); 


  for(i=st; i<=stp; i++)
        { 
 	
          for(j=1; j<=xsize-1; j++)
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
             by<<endl;
             Bt<<endl; 
  
         
        }     
              
     }    

  
  
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
                                       // have to check ..// this has to be taken care of ../
   if(bc1==bc2 && bc1==1)
   {

     if( prcs == np-1 ){     
         
       
           // cout<<stp-1<<"\t"<<stp+1<<endl;
           	     
	   MPI_Sendrecv( u.U1[stp-1],  xpoints, MPI_DOUBLE, 0, BC93,u.U1[stp+1], xpoints, MPI_DOUBLE, 0, PC1, MPI_COMM_WORLD, &status );

	   MPI_Sendrecv( u.U2[stp-1],  xpoints, MPI_DOUBLE, 0, BC2, u.U2[stp+1], xpoints, MPI_DOUBLE, 0, PC2, MPI_COMM_WORLD, &status );

	   MPI_Sendrecv( u.U3[stp-1],  xpoints, MPI_DOUBLE, 0, BC3, u.U3[stp+1], xpoints, MPI_DOUBLE, 0, PC3, MPI_COMM_WORLD, &status );
	   MPI_Sendrecv( u.U4[stp-1],  xpoints, MPI_DOUBLE, 0, BC4, u.U4[stp+1], xpoints, MPI_DOUBLE, 0, PC4, MPI_COMM_WORLD, &status );
	   MPI_Sendrecv( u.U5[stp-1],  xpoints, MPI_DOUBLE, 0, BC5, u.U5[stp+1], xpoints, MPI_DOUBLE, 0, PC5, MPI_COMM_WORLD, &status );
 	   MPI_Sendrecv( u.U6[stp-1],  xpoints, MPI_DOUBLE, 0, BC6, u.U6[stp+1], xpoints, MPI_DOUBLE, 0, PC6, MPI_COMM_WORLD, &status );
	   MPI_Sendrecv( u.U7[stp-1],  xpoints, MPI_DOUBLE, 0, BC7, u.U7[stp+1], xpoints, MPI_DOUBLE, 0, PC7, MPI_COMM_WORLD, &status );
	   MPI_Sendrecv( u.U8[stp-1],  xpoints, MPI_DOUBLE, 0, BC8, u.U8[stp+1], xpoints, MPI_DOUBLE, 0, PC8, MPI_COMM_WORLD, &status );
	   MPI_Sendrecv( q.a[stp-1],  xpoints, MPI_DOUBLE, 0, BC9,  q.a[stp+1], xpoints, MPI_DOUBLE, 0, PC9, MPI_COMM_WORLD, &status );
	   MPI_Sendrecv( q.h[stp-1],  xpoints, MPI_DOUBLE, 0, BC91, q.h[stp+1], xpoints, MPI_DOUBLE, 0, PC91, MPI_COMM_WORLD, &status );
	   MPI_Sendrecv( q.p[stp-1],  xpoints, MPI_DOUBLE, 0, BC92, q.p[stp+1], xpoints, MPI_DOUBLE, 0, PC92, MPI_COMM_WORLD, &status );
 
            
              }
      

   ///        cout<<stp-1<<"\t"<<stp+1<<endl; 
    
    if( prcs == 0 ){       
        

	    MPI_Sendrecv( u.U1[2], xpoints, MPI_DOUBLE, np-1, PC1, u.U1[0], xpoints, MPI_DOUBLE, np-1, BC93, MPI_COMM_WORLD, &status );
	    MPI_Sendrecv( u.U2[2], xpoints, MPI_DOUBLE, np-1, PC2, u.U2[0], xpoints, MPI_DOUBLE, np-1, BC2, MPI_COMM_WORLD, &status );
	    MPI_Sendrecv( u.U3[2], xpoints, MPI_DOUBLE, np-1, PC3, u.U3[0], xpoints, MPI_DOUBLE, np-1, BC3, MPI_COMM_WORLD, &status );
	    MPI_Sendrecv( u.U4[2], xpoints, MPI_DOUBLE, np-1, PC4, u.U4[0], xpoints, MPI_DOUBLE, np-1, BC4, MPI_COMM_WORLD, &status );
	    MPI_Sendrecv( u.U5[2], xpoints, MPI_DOUBLE, np-1, PC5, u.U5[0], xpoints, MPI_DOUBLE, np-1, BC5, MPI_COMM_WORLD, &status );
	    MPI_Sendrecv( u.U6[2], xpoints, MPI_DOUBLE, np-1, PC6, u.U6[0], xpoints, MPI_DOUBLE, np-1, BC6, MPI_COMM_WORLD, &status );
            MPI_Sendrecv( u.U7[2], xpoints, MPI_DOUBLE, np-1, PC7, u.U7[0], xpoints, MPI_DOUBLE, np-1, BC7, MPI_COMM_WORLD, &status );
            MPI_Sendrecv( u.U8[2], xpoints, MPI_DOUBLE, np-1, PC8, u.U8[0], xpoints, MPI_DOUBLE, np-1, BC8, MPI_COMM_WORLD, &status );
            MPI_Sendrecv( q.a[2], xpoints, MPI_DOUBLE, np-1, PC9,  q.a[0], xpoints, MPI_DOUBLE, np-1, BC9, MPI_COMM_WORLD, &status );
            MPI_Sendrecv( q.h[2], xpoints, MPI_DOUBLE, np-1, PC91, q.h[0], xpoints, MPI_DOUBLE, np-1, BC91, MPI_COMM_WORLD, &status );
            MPI_Sendrecv( q.p[2], xpoints, MPI_DOUBLE, np-1, PC92, q.p[0], xpoints, MPI_DOUBLE, np-1, BC92, MPI_COMM_WORLD, &status );
  



                    }
           // cout<<stp-1<<"\t"<<stp+1<<endl;
   
         }

       
             //  cout<<stp-1<<"\t"<<stp+1<<endl;
 	
                if(bc4==1)   // left and right ..// not to be changed 
                {

                   for(i=st; i<=stp; i++)
                   { 

                              

                     q.p[i][xsize] = q.p[i][2];           
           
          
                     q.h[i][xsize] = q.h[i][2];            
          
                  
                     q.a[i][xsize] = q.a[i][2];            
           
                     u.U1[i][xsize] = u.U1[i][2];         
                     u.U2[i][xsize] = u.U2[i][2];         
	             u.U3[i][xsize] = u.U3[i][2];         
	             u.U4[i][xsize] = u.U4[i][2];         
 	             u.U5[i][xsize] = u.U5[i][2];         
 	             u.U6[i][xsize] = u.U6[i][2];         
 	             u.U7[i][xsize] = u.U7[i][2];         
 	             u.U8[i][xsize] = u.U8[i][2];         
           
                     }

                  }
             
                if(bc3==1)
                {
  
                  for(i=st; i<=stp; i++)
                   { 


                    q.p[i][0] = q.p[i][xsize-2];           
          
                    q.h[i][0]= q.h[i][xsize-2];           
          
                    q.a[i][0] = q.a[i][xsize-2] ;          
          
                 
  
                    u.U1[i][0] = u.U1[i][xsize-2] ;    
     
                    u.U2[i][0] = u.U2[i][xsize-2];    
     
                    u.U3[i][0] = u.U3[i][xsize-2];    
     
                    u.U4[i][0] = u.U4[i][xsize-2];    
     
                    u.U5[i][0] = u.U5[i][xsize-2];    
	            u.U6[i][0] = u.U6[i][xsize-2];    
                    u.U7[i][0] = u.U7[i][xsize-2];    
                    u.U8[i][0] = u.U8[i][xsize-2];  
                  

                      }
  
                  }
              
           
            

      }
  
  

void grid()
{

   int l=0,m=0, start, stop;
    

   y = new double[PartitionSize];


 

	if( prcs==0)
	start=1;

	else{

	 start= width*prcs+1;  
         m= start;

         }  
	   



     for(j=1; j<=xsize-1; j++)
          {
 	
          
	     
  x[j]=  -0.0+ dx*(l);
 
           l++;
 
   
          } 
  

//    m = start-1;

       for(i=st; i<=stp; i++)
         {
  
 
  
   y[i]= -.0 + dy*(m);
           m++;
            
         }
   
    }

int offset1(int prcs)
{

if(prcs==0)
   return 1;

  

 else 
   return st;


    
 }

int offset2(int prcs)
{

if(prcs==0)
 {

  return stp;

  }

 else if(prcs==np-1 )
 {

   return stp-1;


    }

 else 
 return stp;

   }




  



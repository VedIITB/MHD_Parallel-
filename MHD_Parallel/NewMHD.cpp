#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fstream>
#include <iostream>
#include "mesh.h"
#include <mpi.h>
#include <time.h>
#include <sstream>


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

double x[xpoints],y[ypoints], t1, t2;

int prcs,np, width , st, stp, ymax,start, stop;

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
   width=ymax/np;
   st=0;
   stp=width; 
 
   for(i=1; i<=np-1; i++){
  
   MPI_Send( &width,  1, MPI_INT, prcs+i, DOWN,   MPI_COMM_WORLD);
   

          }

    }

  grid();
  MPI_Bcast(&tf, 1, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);
  MPI_Bcast(&CFL, 1, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);
  MPI_Bcast(&bc1, 1, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);
  MPI_Bcast(&bc2, 1, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);
  MPI_Bcast(&bc3, 1, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);
  MPI_Bcast(&bc4, 1, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);
 
  
if(prcs!=0){


   MPI_Recv( &width, 1 , MPI_INT, MPI_ANY_SOURCE, DOWN  ,  MPI_COMM_WORLD, &status);

   if(prcs!=np-1) 
    { 
                                    // will be in a loop separately ./
   st=width*prcs + 1;

   stp=width*(prcs+1);

       }

  else if (prcs==np-1){

   st=width*prcs + 1;

   stp=width*(prcs+1)+ (ymax%np-1);

     }

 }


 
    initilization();
    boundaryConditions();

	   do{

             RK3();
	    
          if(prcs==0)
           time_step();
     
           MPI_Bcast(&dt, 1, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);
           MPI_Bcast(&t, 1, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);  
   
       
	  if(prcs==0)
           tracer();
 	   
          }	   while(t<=tf);
 

   t2 = clock();   	 
 
  MPI_Barrier( MPI_COMM_WORLD);
   cout<<"THE END OF THE PROGRAMME IN =\t"<<(t2-t1)/CLOCKS_PER_SEC<<"\tsecs"<<endl;


   for(i=st; i<=stp; i++)
    {

	     

     MPI_Send( &q.rho[i][0], xpoints  , MPI_DOUBLE, 0, R,   MPI_COMM_WORLD);
     MPI_Send( &q.p[i][0], xpoints  , MPI_DOUBLE, 0, P,   MPI_COMM_WORLD);
     MPI_Send( &q.vx[i][0], xpoints  , MPI_DOUBLE, 0, VX,   MPI_COMM_WORLD);
     MPI_Send( &q.vy[i][0], xpoints  , MPI_DOUBLE, 0, VY,   MPI_COMM_WORLD);
     MPI_Send( &q.by[i][0], xpoints  , MPI_DOUBLE, 0, BY,   MPI_COMM_WORLD);
     MPI_Send( &q.b[i][0], xpoints  , MPI_DOUBLE, 0, tb,   MPI_COMM_WORLD);
    
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
  
	   /*
  
            if(x[j]<=0.5)
      {
               q.p[i][j] = 1;
	       q.rho[i][j] = 1; 	         
                q.by[i][j] = 1;//5.0/(sqrt(4*3.142)); 
                q.vx[i][j] = 0;
	        q.bz[i][j]  = 0; // 5/(sqrt(4*3.142));
      q.bx[i][j]  =  .75; // 5/(sqrt(4*3.142));
      q.vy[i][j]  =  0; // 5/(sqrt(4*3.142));
      q.vz[i][j]  =  0; // 5/(sqrt(4*3.142));
      


		 }
	  else
	{
	
	
		
              q.p[i][j] = .1;
	       q.rho[i][j] = 0.125; 	         
                q.by[i][j] = -1;//5.0/(sqrt(4*3.142)); 
               q.vx[i][j]  =  0; // 5/(sqrt(4*3.142));
          q.bz[i][j]  = 0; // 5/(sqrt(4*3.142));
      q.bx[i][j]  =  .75; // 5/(sqrt(4*3.142));
      q.vy[i][j]  =  0; // 5/(sqrt(4*3.142));
      q.vz[i][j]  =  0; // 5/(sqrt(4*3.142));
      
 
     }
    }
	
   */
  
 /*
    r = sqrt((x[j]-0.5) * (x[j]-0.5) + (y[i]-0.5) * (y[i]-0.5));
	     
		    if(r<0.1)
	        {
	    		f=(0.115-r)/(0.015);
			
	     q.p[i][j] = 1;
	
	    q.rho[i][j] = 10; 	         
    
	    q.by[i][j] = 0; 
  	
	    q.vx[i][j] = -20*(y[i]-0.5);
      
	    q.vy[i][j] = 20*(x[j]-0.5);
      
	  	q.vz[i][j] = 0;
	  
	  	q.bx[i][j]         = 5/sqrt(4*3.142) ;
	  
	    q.bz[i][j]        = 0 ;
	  
}
 

	  
	 else if(r> r_1)
	 {
	    	  q.p[i][j] = 1;
	
	    q.rho[i][j] = 1; 	         
    
	    q.by[i][j] = 0.0; 
  	
	    q.vx[i][j] = 0;
      
	    q.vy[i][j] = 0;
      
	  	q.vz[i][j] = 0;
	  
	  	q.bz[i][j]         = 0.0 ;
	  
	    q.bx[i][j]        =  5/sqrt(4*3.142) ;
	  
	
	 }
		else 
	{
			
			f=(0.115-r)/(0.015);
	        
	     q.p[i][j] = 1.0;
	
	    q.rho[i][j] = 1+9.0*f; 	         
    
	    q.by[i][j] = 0.0; 
  	
	    q.vx[i][j] = -f*2*(y[i]-0.5)/r;
      
	    q.vy[i][j] = f*2*(x[j]-0.5)/r;
      
	  	q.vz[i][j] = 0;
	  
	  	q.bz[i][j]         = 0;
	  
	    q.bx[i][j]        =  5/sqrt(4*3.142);
	 
	 
	 
			
	   
	 } 
	 }	     
*/
    
 
          

	 

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
boundaryConditions();

  
  }

void wave_system()
{

double temp1,alfven;
for(i=st; i<=stp; i++)
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
     
 }



void eigen_value()
{
      
   
 for(i=st; i<=stp; i++)
 {
 	
    for(j=0; j<=xsize; j++)
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

     

      for(i=st; i<=stp; i++)
        { 
 	
          for(j=0; j<=xsize; j++)
            {


               if(fabs(c.ly2[i][j])>fabs(c.ly1[i][j]))
                 c.lMaxy[i][j] = fabs(c.ly2[i][j]);

               else 
                 c.lMaxy [i][j]= fabs(c.ly1[i][j]);
  
            
 
             }   
  

         }


stop=offset2(prcs);  



    if( prcs != 0 )    
         MPI_Send( &c.lMaxy[st][0],  xpoints, MPI_DOUBLE, prcs-1, SENDY,   MPI_COMM_WORLD);  




   
  if( prcs != np-1 )       
        MPI_Recv( &c.lMaxy[stop+1][0], xpoints, MPI_DOUBLE, MPI_ANY_SOURCE,SENDY, MPI_COMM_WORLD, &status);
     

      for(i=st; i<=stop; i++)
        { 
 	
          for(j=0; j<=xsize-1; j++)
            {


               if(c.lMaxx[i][j+1] > c.lMaxx[i][j] )
                 c.lx[i][j] = c.lMaxx[i][j+1];

               else
                 c.lx[i][j] = c.lMaxx[i][j];
     
 
             }   
  

         }

      for(i=st; i<=stop; i++)
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
                         MPI_Send( &u.U1[st][0],  xpoints, MPI_DOUBLE, prcs-1, S1,   MPI_COMM_WORLD);  
                         MPI_Send( &u.U2[st][0],  xpoints, MPI_DOUBLE, prcs-1, S2,   MPI_COMM_WORLD);   
                         MPI_Send( &u.U3[st][0],  xpoints, MPI_DOUBLE, prcs-1, S3,   MPI_COMM_WORLD);   
                         MPI_Send( &u.U4[st][0],  xpoints, MPI_DOUBLE, prcs-1, S4,   MPI_COMM_WORLD);   
                        
                         MPI_Send( &u.U5[st][0],  xpoints, MPI_DOUBLE, prcs-1, S5,   MPI_COMM_WORLD);  
                         MPI_Send( &u.U6[st][0],  xpoints, MPI_DOUBLE, prcs-1, S6,   MPI_COMM_WORLD);  
                         MPI_Send( &u.U7[st][0],  xpoints, MPI_DOUBLE, prcs-1, S7,   MPI_COMM_WORLD);  
                         MPI_Send( &u.U8[st][0],  xpoints, MPI_DOUBLE, prcs-1, S8,   MPI_COMM_WORLD);  
      
        
 
 
              }
      

  
   
    if( prcs != np-1 ){       
                      MPI_Recv( &u.U1[stop+1][0], xpoints, MPI_DOUBLE, MPI_ANY_SOURCE,S1, MPI_COMM_WORLD, &status);
                      MPI_Recv( &u.U2[stop+1][0], xpoints, MPI_DOUBLE, MPI_ANY_SOURCE,S2, MPI_COMM_WORLD, &status);
                      MPI_Recv( &u.U3[stop+1][0], xpoints, MPI_DOUBLE, MPI_ANY_SOURCE,S3, MPI_COMM_WORLD, &status);
                      MPI_Recv( &u.U4[stop+1][0], xpoints, MPI_DOUBLE, MPI_ANY_SOURCE,S4, MPI_COMM_WORLD, &status);
                     
                      MPI_Recv( &u.U5[stop+1][0], xpoints, MPI_DOUBLE, MPI_ANY_SOURCE,S5, MPI_COMM_WORLD, &status);   
                      MPI_Recv( &u.U6[stop+1][0], xpoints, MPI_DOUBLE, MPI_ANY_SOURCE,S6, MPI_COMM_WORLD, &status);   
                      MPI_Recv( &u.U7[stop+1][0], xpoints, MPI_DOUBLE, MPI_ANY_SOURCE,S7, MPI_COMM_WORLD, &status);   
                      MPI_Recv( &u.U8[stop+1][0], xpoints, MPI_DOUBLE, MPI_ANY_SOURCE,S8, MPI_COMM_WORLD, &status);   
               
 

 
                     }


    

 if( prcs != 0 ){     
         
                         MPI_Send( &q.h[st][0],  xpoints, MPI_DOUBLE, prcs-1, q4,   MPI_COMM_WORLD);  
                         MPI_Send( &q.p[st][0],  xpoints, MPI_DOUBLE, prcs-1, q5,   MPI_COMM_WORLD);                                      
                        
                     

 
              }

  
   
    if( prcs != np-1 ){       
          
                      MPI_Recv( &q.h[stop+1][0], xpoints, MPI_DOUBLE, MPI_ANY_SOURCE,q4, MPI_COMM_WORLD, &status);    
                      MPI_Recv( &q.p[stop+1][0], xpoints, MPI_DOUBLE, MPI_ANY_SOURCE,q5, MPI_COMM_WORLD, &status);   
                   
                     }
          

             
  for(i=st; i<=stop; i++)
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

for(i=st; i<=stp; i++)
 {
 	
    for(j=0; j<=xsize; j++)
     {
 

  

    

                NetdivBx=  -(u.U7[i][j+1]-u.U7[i][j]); 
   
                NetdivBy =  -(u.U6[i+1][j]-u.U6[i][j]);


                S.s1[i][j]=0;
        
		S.s1[i][j] = (NetdivBx+ NetdivBy)*u.U7[i][j]; 
		
		S.s3[i][j] = (NetdivBx+ NetdivBy)*u.U6[i][j]; 
		 
		
		S.s4[i][j] = (NetdivBx+ NetdivBy)*u.U8[i][j]; 
		 
		
		S.s5[i][j] = ((u.U7[i][j]*q.vx[i][j])+(u.U6[i][j]*q.vy[i][j])+(u.U8[i][j]*q.vz[i][j]))*(NetdivBx +NetdivBx); 
		
		S.s8[i][j] = q.vz[i][j]*(NetdivBx+ NetdivBy);  
		
		S.s7[i][j] = q.vx[i][j]*(NetdivBx+ NetdivBy);
        
		S.s6[i][j] =  q.vy[i][j]*(NetdivBx+ NetdivBy);
  
  
          }

       }
    
   }
    
void RK3()
{


   iteration();

start= offset1(prcs);
stop = offset2(prcs);

if( prcs != np-1 ){     
         
                         MPI_Send( &f.G1[stp][0],  xpoints, MPI_DOUBLE, prcs+1, q7,   MPI_COMM_WORLD);  
                         MPI_Send( &f.G2[stp][0],  xpoints, MPI_DOUBLE, prcs+1, q8,   MPI_COMM_WORLD);  
                         MPI_Send( &f.G3[stp][0],  xpoints, MPI_DOUBLE, prcs+1, q9,   MPI_COMM_WORLD);  
                         MPI_Send( &f.G4[stp][0],  xpoints, MPI_DOUBLE, prcs+1, q10,   MPI_COMM_WORLD);  
                         MPI_Send( &f.G5[stp][0],  xpoints, MPI_DOUBLE, prcs+1, q11,   MPI_COMM_WORLD);  
                         MPI_Send( &f.G6[stp][0],  xpoints, MPI_DOUBLE, prcs+1, q12,   MPI_COMM_WORLD);  
                         MPI_Send( &f.G7[stp][0],  xpoints, MPI_DOUBLE, prcs+1, q13,   MPI_COMM_WORLD);  
                         MPI_Send( &f.G8[stp][0],  xpoints, MPI_DOUBLE, prcs+1, q14,   MPI_COMM_WORLD);  
              
                         
                 }
      

  
   
    if( prcs != 0 ){       
                      MPI_Recv( &f.G1[st-1][0], xpoints, MPI_DOUBLE, MPI_ANY_SOURCE,q7, MPI_COMM_WORLD, &status);
                      MPI_Recv( &f.G2[st-1][0], xpoints, MPI_DOUBLE, MPI_ANY_SOURCE,q8, MPI_COMM_WORLD, &status);
                      MPI_Recv( &f.G3[st-1][0], xpoints, MPI_DOUBLE, MPI_ANY_SOURCE,q9, MPI_COMM_WORLD, &status);    
                      MPI_Recv( &f.G4[st-1][0], xpoints, MPI_DOUBLE, MPI_ANY_SOURCE,q10, MPI_COMM_WORLD, &status);   
                      MPI_Recv( &f.G5[st-1][0], xpoints, MPI_DOUBLE, MPI_ANY_SOURCE,q11, MPI_COMM_WORLD, &status);   
                      MPI_Recv( &f.G6[st-1][0], xpoints, MPI_DOUBLE, MPI_ANY_SOURCE,q12, MPI_COMM_WORLD, &status);   
                      MPI_Recv( &f.G7[st-1][0], xpoints, MPI_DOUBLE, MPI_ANY_SOURCE,q13, MPI_COMM_WORLD, &status);   
                      MPI_Recv( &f.G8[st-1][0], xpoints, MPI_DOUBLE, MPI_ANY_SOURCE,q14, MPI_COMM_WORLD, &status);   
 
                    


                    }

 
             


  for(i=start; i<=stop; i++)
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
 
   for(i=start; i<=stop; i++)
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
 

  for(i=start; i<=stop; i++)
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


                if(i==stp) {
                 cout<<"==============================================================================="<<endl;
                 cout<<"EVOLVING FOR TIME :: "<<t<<endl;
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
  MPI_Recv( &lamdaMAX, 1, MPI_DOUBLE, MPI_ANY_SOURCE,EIGEN, MPI_COMM_WORLD, &status);
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
 for(i=1; i<=stp; i++)
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


for(int u=1; u<=np-1; u++){ 

  
    st=width*u + 1;

   stp=width*(u+1);

    
  if(u==np-1)
   stp=width*(u+1)+(ymax%np-1);


   cout<<"WRITING FOR PROCESS-"<<"\t"<<u+1<<endl;

   for(i=st; i<=stp; i++)
  {
     
     MPI_Recv( &q.rho[i][0], xpoints , MPI_DOUBLE, u, R  ,  MPI_COMM_WORLD, &status);
     MPI_Recv( &q.p[i][0], xpoints , MPI_DOUBLE, u, P  ,  MPI_COMM_WORLD, &status);
     MPI_Recv( &q.vx[i][0], xpoints , MPI_DOUBLE, u, VX  ,  MPI_COMM_WORLD, &status);
     MPI_Recv( &q.vy[i][0], xpoints , MPI_DOUBLE, u, VY  ,  MPI_COMM_WORLD, &status);
     MPI_Recv( &q.by[i][0], xpoints , MPI_DOUBLE, u, BY  ,  MPI_COMM_WORLD, &status);
     MPI_Recv( &q.b[i][0], xpoints , MPI_DOUBLE, u,   tb ,  MPI_COMM_WORLD, &status);
          
            }

   

if(u==np-1)
stp=stp-1;

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

if(bc1==bc2 && bc1==1){

if( prcs == np-1 ){     
         
                     
   MPI_Sendrecv( &u.U1[xsize-2][0],  xpoints, MPI_DOUBLE, 0, BC93, &u.U1[xsize][0], xpoints, MPI_DOUBLE, 0, PC1, MPI_COMM_WORLD, &status );

   MPI_Sendrecv( &u.U2[xsize-2][0],  xpoints, MPI_DOUBLE, 0, BC2, &u.U2[xsize][0], xpoints, MPI_DOUBLE, 0, PC2, MPI_COMM_WORLD, &status );

   MPI_Sendrecv( &u.U3[xsize-2][0],  xpoints, MPI_DOUBLE, 0, BC3, &u.U3[xsize][0], xpoints, MPI_DOUBLE, 0, PC3, MPI_COMM_WORLD, &status );
   MPI_Sendrecv( &u.U4[xsize-2][0],  xpoints, MPI_DOUBLE, 0, BC4, &u.U4[xsize][0], xpoints, MPI_DOUBLE, 0, PC4, MPI_COMM_WORLD, &status );
   MPI_Sendrecv( &u.U5[xsize-2][0],  xpoints, MPI_DOUBLE, 0, BC5, &u.U5[xsize][0], xpoints, MPI_DOUBLE, 0, PC5, MPI_COMM_WORLD, &status );
   MPI_Sendrecv( &u.U6[xsize-2][0],  xpoints, MPI_DOUBLE, 0, BC6, &u.U6[xsize][0], xpoints, MPI_DOUBLE, 0, PC6, MPI_COMM_WORLD, &status );
   MPI_Sendrecv( &u.U7[xsize-2][0],  xpoints, MPI_DOUBLE, 0, BC7, &u.U7[xsize][0], xpoints, MPI_DOUBLE, 0, PC7, MPI_COMM_WORLD, &status );
   MPI_Sendrecv( &u.U8[xsize-2][0],  xpoints, MPI_DOUBLE, 0, BC8, &u.U8[xsize][0], xpoints, MPI_DOUBLE, 0, PC8, MPI_COMM_WORLD, &status );
   MPI_Sendrecv( &q.a[xsize-2][0],  xpoints, MPI_DOUBLE, 0, BC9, &q.a[xsize][0], xpoints, MPI_DOUBLE, 0, PC9, MPI_COMM_WORLD, &status );
   MPI_Sendrecv( &q.h[xsize-2][0],  xpoints, MPI_DOUBLE, 0, BC91, &q.h[xsize][0], xpoints, MPI_DOUBLE, 0, PC91, MPI_COMM_WORLD, &status );
   MPI_Sendrecv( &q.p[xsize-2][0],  xpoints, MPI_DOUBLE, 0, BC92, &q.p[xsize][0], xpoints, MPI_DOUBLE, 0, PC92, MPI_COMM_WORLD, &status );
 
  }
      

  
   
   else if( prcs == 0 ){       
        

    MPI_Sendrecv( &u.U1[2][0], xpoints, MPI_DOUBLE, np-1, PC1, &u.U1[0][0], xpoints, MPI_DOUBLE, np-1, BC93, MPI_COMM_WORLD, &status );
    MPI_Sendrecv( &u.U2[2][0], xpoints, MPI_DOUBLE, np-1, PC2, &u.U2[0][0], xpoints, MPI_DOUBLE, np-1, BC2, MPI_COMM_WORLD, &status );
  MPI_Sendrecv( &u.U3[2][0], xpoints, MPI_DOUBLE, np-1, PC3, &u.U3[0][0], xpoints, MPI_DOUBLE, np-1, BC3, MPI_COMM_WORLD, &status );
  MPI_Sendrecv( &u.U4[2][0], xpoints, MPI_DOUBLE, np-1, PC4, &u.U4[0][0], xpoints, MPI_DOUBLE, np-1, BC4, MPI_COMM_WORLD, &status );
  MPI_Sendrecv( &u.U5[2][0], xpoints, MPI_DOUBLE, np-1, PC5, &u.U5[0][0], xpoints, MPI_DOUBLE, np-1, BC5, MPI_COMM_WORLD, &status );
  MPI_Sendrecv( &u.U6[2][0], xpoints, MPI_DOUBLE, np-1, PC6, &u.U6[0][0], xpoints, MPI_DOUBLE, np-1, BC6, MPI_COMM_WORLD, &status );
  MPI_Sendrecv( &u.U7[2][0], xpoints, MPI_DOUBLE, np-1, PC7, &u.U7[0][0], xpoints, MPI_DOUBLE, np-1, BC7, MPI_COMM_WORLD, &status );
  MPI_Sendrecv( &u.U8[2][0], xpoints, MPI_DOUBLE, np-1, PC8, &u.U8[0][0], xpoints, MPI_DOUBLE, np-1, BC8, MPI_COMM_WORLD, &status );
  MPI_Sendrecv( &q.a[2][0], xpoints, MPI_DOUBLE, np-1, PC9, &q.a[0][0], xpoints, MPI_DOUBLE, np-1, BC9, MPI_COMM_WORLD, &status );
  MPI_Sendrecv( &q.h[2][0], xpoints, MPI_DOUBLE, np-1, PC91, &q.h[0][0], xpoints, MPI_DOUBLE, np-1, BC91, MPI_COMM_WORLD, &status );
  MPI_Sendrecv( &q.p[2][0], xpoints, MPI_DOUBLE, np-1, PC92, &q.p[0][0], xpoints, MPI_DOUBLE, np-1, BC92, MPI_COMM_WORLD, &status );
  



                    }

   
         }

 for(i=st; i<=stp; i++)
        { 
 	
          for(j=0; j<=xsize; j++)
            {

              
 
          if(bc1==0)
                {
  
                    q.rho[xsize][j] = q.rho[xsize-1][j] ;           
           
                    q.p[xsize][j] = q.p[xsize-1][j] ;          
          
                    q.vx[xsize][j] = q.vx[xsize-1][j];            
          
                    q.vy[xsize][j] = q.vy[xsize-1][j] ;           
          
                    q.h[xsize][j] = q.h[xsize-1][j]   ;         
          
                  //  q.M[xsize][j] = q.M[xsize-1][j]  ;          
          
                    q.a[xsize][j] = q.a[xsize-1][j]  ;  
                  
                    q.vz[xsize][j] = q.vz[xsize-1][j]  ;  
           
                  q.bx[xsize][j] = q.bx[xsize-1][j]  ;  
          
                   q.by[xsize][j] = q.by[xsize-1][j]  ;  
                   
                    q.bz[xsize][j] = q.bz[xsize-1][j]  ;  
        
          
                 
                    u.U1[xsize][j] = u.U1[xsize-1][j]  ;          

                    u.U2[xsize][j] = u.U2[xsize-1][j]  ;          
              
                    u.U3[xsize][j] = u.U3[xsize-1][j]  ;          

                    u.U4[xsize][j] = u.U4[xsize-1][j]  ;          
                   
                     u.U5[xsize][j] = u.U5[xsize-1][j]  ;          
                   
                     u.U6[xsize][j] = u.U6[xsize-1][j]  ;          
               
                    u.U7[xsize][j] = u.U7[xsize-1][j]  ;          
                   
                   u.U8[xsize][j] = u.U8[xsize-1][j]  ;          
                   


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
  
                    q.rho[i][xsize] = q.rho[i][xsize-1] ;           
           
                    q.p[i][xsize] = q.p[i][xsize-1]  ;         
          
                    q.vx[i][xsize] = q.vx[i][xsize-1]  ;          
          
                    q.vy[i][xsize] = q.vy[i][xsize-1]  ;          
          
                    q.h[i][xsize] = q.h[i][xsize-1]   ;         
          
                  //  q.M[i][xsize] = q.M[i][xsize-1]  ;          
          
                    q.a[i][xsize] = q.a[i][xsize-1]  ;          
          
                     q.vz[i][xsize] = q.vz[i][xsize-1]  ;          
             q.by[i][xsize] = q.by[i][xsize-1]  ;          
             q.bx[i][xsize] = q.bx[i][xsize-1]  ;          
             q.bz[i][xsize] = q.bz[i][xsize-1]  ;          
             
                 
                    u.U1[i][xsize] = u.U1[i][xsize-1]  ;          
           
                    u.U2[i][xsize] = u.U2[i][xsize-1]  ;  
  
                    u.U3[i][xsize] = u.U3[i][xsize-1]  ;  

                    u.U4[i][xsize] = u.U4[i][xsize-1]  ;  
 
                   u.U5[i][xsize] = u.U5[i][xsize-1]  ;  
  u.U6[i][xsize] = u.U6[i][xsize-1]  ;  
  u.U7[i][xsize] = u.U7[i][xsize-1]  ;  
  u.U8[i][xsize] = u.U8[i][xsize-1]  ;  
 

                  }
               
             else if(bc4==1)
                {

            
                // not here  ..//
                    q.p[i][xsize] = q.p[i][2]     ;      
          
          
                    q.h[i][xsize] = q.h[i][2]     ;       
          
                   // q.M[xsize][j] = q.M[1][j]   ;         
          
                    q.a[i][xsize] = q.a[i][2]   ;         
           
                     u.U1[i][xsize] = u.U1[i][2]   ;         
                     u.U2[i][xsize] = u.U2[i][2]   ;         
             u.U3[i][xsize] = u.U3[i][2]   ;         
             u.U4[i][xsize] = u.U4[i][2]   ;         
             u.U5[i][xsize] = u.U5[i][2]   ;         
             u.U6[i][xsize] = u.U6[i][2]   ;         
             u.U7[i][xsize] = u.U7[i][2]   ;         
             u.U8[i][xsize] = u.U8[i][2]   ;         
           
            
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
  
                    q.rho[0][j] = q.rho[1][j]   ;         
           
                    q.p[0][j] = q.p[1][j]   ;        
          
                    q.vx[0][j] = q.vx[1][j]   ;         
          
                    q.vy[0][j] = q.vy[1][j]   ;         
          
                    q.h[0][j] = q.h[1][j]  ;          
          
                //    q.M[0][j] = q.M[1][j]   ;         
          
                    q.a[0][j] = q.a[1][j]   ;    

                     q.vz[0][j] = q.vz[1][j]   ;    
                    q.by[0][j] = q.by[1][j]   ;    
                   q.bx[0][j] = q.bx[1][j]   ;    
                   q.bz[0][j] = q.bz[1][j]   ;    

                    u.U1[0][j] = u.U1[1][j]   ;    
     
                    u.U2[0][j] = u.U2[1][j]   ;    
     
                    u.U3[0][j] = u.U3[1][j]   ;    
     
                    u.U4[0][j] = u.U4[1][j]   ;    
     
                     u.U5[0][j] = u.U5[1][j]   ;    
                     u.U6[0][j] = u.U6[1][j]   ;    
                     u.U7[0][j] = u.U7[1][j]   ;    
                     u.U8[0][j] = u.U8[1][j]   ;    
     


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
  
                 
                    q.p[i][0] = q.p[i][1];           
          
                 
                    q.h[i][0]= q.h[i][1] ;           
          
                 //   q.M[i][0] = q.M[i][1] ;           
          
                    q.a[i][0] = q.a[i][1]  ;          
          
                 
  
                  u.U1[i][0] = u.U1[i][1]   ;    
     
                  u.U2[i][0] = u.U2[i][1]   ;    
     
                 u.U3[i][0] = u.U3[i][1]  ;    
     
                 u.U4[i][0] = u.U4[i][1]   ;    
     
                 u.U5[i][0] = u.U5[i][1]   ;    
       u.U6[i][0] = u.U6[i][1]   ;    
       u.U7[i][0] = u.U7[i][1]   ;    
       u.U8[i][0] = u.U8[i][1]   ;    
     

               }
               
             else if(bc3==1)
                {
  
                  q.p[i][0] = q.p[i][xsize-2];           
          
                 
                    q.h[i][0]= q.h[i][xsize-2] ;           
          
                 //   q.M[i][0] = q.M[i][1] ;           
          
                    q.a[i][0] = q.a[i][xsize-2]  ;          
          
                 
  
                  u.U1[i][0] = u.U1[i][xsize-2]   ;    
     
                  u.U2[i][0] = u.U2[i][xsize-2]   ;    
     
                 u.U3[i][0] = u.U3[i][xsize-2]  ;    
     
                 u.U4[i][0] = u.U4[i][xsize-2]   ;    
     
                 u.U5[i][0] = u.U5[i][xsize-2]   ;    
       u.U6[i][0] = u.U6[i][xsize-2]   ;    
       u.U7[i][0] = u.U7[i][xsize-2]   ;    
       u.U8[i][0] = u.U8[i][xsize-2]   ;  
                  
  
                  }
              
           
            }

      }
  
  }

void grid()
{

   int l=0,m=0;
   
     for(j=1; j<=xsize-1; j++)
          {
 	
          
	     
  x[j]=  -.0+ dx*(l);
 
           l++;
 
   
          } 
  
  
       for(i=1; i<=ysize-1; i++)
         {
  
 
  
   y[i]= -.0 + dy*(m);
           m++;
            
         }
   
    }

int offset1(int prcs)
{

if(prcs==0)
 {

  return 1;

  }

 else if(prcs!=0 )
 {

   return st;


    }



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

 else {

 return stp;

   }

  }




  



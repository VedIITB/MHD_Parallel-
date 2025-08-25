#include "mesh.h"
#include "math.h"


void Primitive(double p[][xpoints], double rho[][xpoints],  double vx[][xpoints],double vy[][xpoints],double vz[][xpoints],double bx[][xpoints],double by[][xpoints], double bz[][xpoints] , double x[], double y[]);




void Primitive(double p[][xpoints], double rho[][xpoints],double vx[][xpoints],double vy[][xpoints],double vz[][xpoints], double bx[][xpoints],double by[][xpoints],double bz[][xpoints],double x[] , double y[] ){

 double r,f,r_0=0.1, r_1=0.115; 
for(i=3; i<=ysize-3; i++)
    {
 	
    for(j=3; j<=xsize-3; j++)
          
        {
	  
            
          p[i][j] = 1/GAMMA;
	
	    rho[i][j] = 1; 	         
    
	    by[i][j] = sin(4*3.142*x[j])/GAMMA; 
  	
	    vx[i][j] = -sin(2*3.142*y[i]);
      
	    vy[i][j] = sin(2*3.142*x[j]);
      
	  	vz[i][j] = 0;
	  
	  	bz[i][j]         = 0 ;
	  
	    bx[i][j]        = -sin(2*3.142*y[i])/GAMMA; 
        }
   
    }

 
   }




#include "mesh.h"
#include "math.h"


void Primitive(double p[][xpoints], double rho[][xpoints],  double vx[][xpoints],double vy[][xpoints],double vz[][xpoints],double bx[][xpoints],double by[][xpoints], double bz[][xpoints] , double x[], double y[]);




void Primitive(double p[][xpoints], double rho[][xpoints],double vx[][xpoints],double vy[][xpoints],double vz[][xpoints], double bx[][xpoints],double by[][xpoints],double bz[][xpoints],double x[] , double y[] ){

double r,f,r_0=0.1, r_1=0.115;
for(i=3; i<=ysize-3; i++)
    {
 	
    for(j=3; j<=xsize-3; j++)
          
        {


        r = sqrt((x[j]-0.5) * (x[j]-0.5) + (y[i]-0.5) * (y[i]-0.5));
	     
		    if(r<0.1)
	        {
	    		f=(0.115-r)/(0.015);
			
	     p[i][j] = 1;
	
	    rho[i][j] = 10; 	         
    
	    by[i][j] = 0; 
  	
	    vx[i][j] = -20*(y[i]-0.5);
      
	    vy[i][j] = 20*(x[j]-0.5);
      
	  	vz[i][j] = 0;
	  
	  	bx[i][j]         = 5/sqrt(4*3.142) ;
	  
	    vz[i][j]        = 0 ;
	  
}
 

	  
	 else if(r> r_1)
	 {
	    	  p[i][j] = 1;
	
	    rho[i][j] = 1; 	         
    
	    by[i][j] = 0.0; 
  	
	    vx[i][j] = 0;
      
	    vy[i][j] = 0;
      
	  	vz[i][j] = 0;
	  
	  	bz[i][j]         = 0.0 ;
	  
	    bx[i][j]        =  5/sqrt(4*3.142) ;
	  
	
	 }
		else 
	{
			
			f=(0.115-r)/(0.015);
	        
	     p[i][j] = 1.0;
	
	    rho[i][j] = 1+9.0*f; 	         
    
	    by[i][j] = 0.0; 
  	
	    vx[i][j] = -f*2*(y[i]-0.5)/r;
      
	    vy[i][j] = f*2*(x[j]-0.5)/r;
      
	  	vz[i][j] = 0;
	  
	  	bz[i][j]         = 0;
	  
	    bx[i][j]        =  5/sqrt(4*3.142);
	 
	 
	 
			
	   
	 } 

 } 

}







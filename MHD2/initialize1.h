#include "mesh.h"
#include "math.h"


void Primitive(double p[][xpoints], double rho[][xpoints],  double vx[][xpoints],double vy[][xpoints],double vz[][xpoints],double bx[][xpoints],double by[][xpoints], double bz[][xpoints] , double x[], double y[]);




void Primitive(double p[][xpoints], double rho[][xpoints],double vx[][xpoints],double vy[][xpoints],double vz[][xpoints], double bx[][xpoints],double by[][xpoints],double bz[][xpoints],double x[] , double y[] ){

double r,f,r_0=0.1, r_1=0.115;
for(i=3; i<=ysize-3; i++)
    {
 	
    for(j=3; j<=xsize-3; j++)
          
        {
/*

	   if(x[j]>=0 && y[i]>=0)
	  {
	  
	          rho[i][j] = 1; 	         
         vx[i][j] = .75;
         vy[i][j] = -0.5;
        bz[i][j]         = 1/sqrt(8*3.142);
	  	    bx[i][j]        =  2/sqrt(8*3.142);
	        by[i][j] = 0.0; 
  	          p[i][j] = 1.0;
	     vz[i][j] = 0;
	  
	 
	}
	
	 else  if(x[j] <0 && y[i] >= 0)
	  {
	  
	          rho[i][j] = 2; 	         
         vx[i][j] = .75;
         vy[i][j] = 0.5;
     bz[i][j]         = 1/sqrt(8*3.142);
	  	    bx[i][j]        =  2/sqrt(8*3.142);
	        by[i][j] = 0.0; 
  	          p[i][j] = 1.0;
	     vz[i][j] = 0;
	  
	    
	}
	 else  if(x[j]< 0 && y[i] <0)
	  {
	  
	          rho[i][j] = 1; 	         
         vx[i][j] = -.75;
         vy[i][j] = 0.5;
     bz[i][j]         = 1/sqrt(8*3.142);
	  	    bx[i][j]        =  2/sqrt(8*3.142);
	        by[i][j] = 0.0; 
  	          p[i][j] = 1.0;
	     vz[i][j] = 0;
	  
	    
	}
	 
	else if(x[j]>= 0 && y[i] <0)
	  {
	  
	          rho[i][j] = 3; 	         
         vx[i][j] = -.75;
         vy[i][j] = -0.5;
     bz[i][j]         = 1/sqrt(8*3.142);
	  	    bx[i][j]        =  2/sqrt(8*3.142);
	        by[i][j] = 0.0; 
  	          p[i][j] = 1.0;
	     vz[i][j] = 0;
	  
	    
	}
	
	
	
	 } 
	 }	         

  }

/*
 r = sqrt((x[j]-0.5) * (x[j]-0.5) + (y[i]-0.5) * (y[i]-0.5));
	     
		    if(r<=0.1)
	        {
	    	
			
	     p[i][j] = 10;
	
	    rho[i][j] = 1; 	         
    
	    by[i][j] = 1/sqrt(2); 
  	
	    vx[i][j] = 0;
      
	    vy[i][j] = 0;
      
	  	vz[i][j] = 0;
	  
	  	bx[i][j]         = 1/sqrt(2) ;
	  
	    bz[i][j]        = 0 ;
	  
}
 		else 
	{
			
		
	        
	     p[i][j] = 0.1;
	
	    rho[i][j] = 1; 	         
    
	    by[i][j] = 1/sqrt(2); 
  	
	    vx[i][j] = 0;
      
	    vy[i][j] = 0;
      
	  	vz[i][j] = 0;
	  
	  	bz[i][j]         = 0;
	  
	    bx[i][j]        =  1/sqrt(2);
	 
	 
	 
			
	   
	 } 


/*
 r = sqrt((x[j]-0.8) * (x[j]-0.8) + (y[i]-0.5) * (y[i]-0.5));
	     
		    if(x[j]<0.6)
	        {
			
	     p[i][j] = 167.35;
	
	    rho[i][j] = 3.87; 	         
    
	    by[i][j] = 2.18; 
  	
	    vx[i][j] = 0;
      
	    vy[i][j] = 0;
      
	  	vz[i][j] = 0;
	  
	  	bz[i][j]         = -2.18 ;
	  
	    bx[i][j]        = 0 ;
	  
}
 

	  
	 else if(r<0.15)
	 {
	     p[i][j] = 1.0;
	
	    rho[i][j] = 10.0; 	         
    
	    by[i][j] = 0.56; 
  	
	    vx[i][j] = -11.25;
      
	    vy[i][j] = 0;
      
	  	vz[i][j] = 0;
	  
	  	bz[i][j]         =0.56 ;
	  
	    bx[i][j]        = 0 ;
	 
	 
	 
	 }
		else 
	{
	
			
	     p[i][j] = 1;
	
	    rho[i][j] = 1; 	         
    
	    by[i][j] = 0.56; 
  	
	    vx[i][j] = -11.25;
      
	    vy[i][j] = 0;
      
	  	vz[i][j] = 0;
	  
	  	bz[i][j]         = 0.56 ;
	  
	    bx[i][j]        = 0 ;
	  
	
 }
*/

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







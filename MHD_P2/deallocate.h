#ifndef DEALLOCATE_H_
#define DEALLOCATE_H_


#include "mesh.h"
 


void  deAllocate(int PartitionGrid);

 void  deAllocate(int PartitionGrid){




    for (int i = 0; i <= PartitionGrid; i++){
      

                       delete[] q.rho[i];
    		       delete[] q.vx[i];
                       delete[] q.vy[i];
                       delete[] q.vz[i];
                       delete[] q.bx[i];
                       delete[] q.by[i];
                       delete[] q.bz[i];
                       delete[] q.b[i];
                       delete[] q.a[i];
                       delete[] q.p[i];
                       delete[] q.h[i];
  

            }

             delete[] q.rho; delete[] q.vx; delete[] q.vy; delete[] q.vz;
    	     delete[] q.bx;  delete[] q.by; delete[] q.bz;
    	     delete[] q.b;   delete[] q.a;  delete[] q.p;  delete[] q.h;        



       
              PartitionGrid++;
		      
              cout<<"\tPG\t"<<PartitionGrid<<endl;

	      q.rho = new double*[PartitionGrid];
 	      q.vx  = new double*[PartitionGrid];
 	      q.vy  = new double*[PartitionGrid];
 	      q.vz  = new double*[PartitionGrid];
 	      q.bx  = new double*[PartitionGrid];
	      q.by  = new double*[PartitionGrid];
 	      q.bz  = new double*[PartitionGrid];
 	      q.b   = new double*[PartitionGrid];
 	      q.a   = new double*[PartitionGrid];
 	      q.p   = new double*[PartitionGrid];
	      q.h   = new double*[PartitionGrid];


            

             for (int i = 0; i <= PartitionGrid; i++){
      

  		          q.rho[i] = new double[xpoints];
 	                  q.vx[i]  = new double[xpoints];
	         	  q.vy[i]  = new double[xpoints];
	         	  q.vz[i]  = new double[xpoints];
   	   		  q.bx[i]  = new double[xpoints];
	        	  q.by[i]  = new double[xpoints];
	        	  q.bz[i]  = new double[xpoints];
       	    	          q.b[i]   = new double[xpoints];
	                  q.a[i]   = new double[xpoints];
	                  q.p[i]   = new double[xpoints];
	                  q.h[i]   = new double[xpoints];



                   }


           }

#endif

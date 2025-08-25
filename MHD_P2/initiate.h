#ifndef INITIATE_H_
#define INITIATE_H_


#include "mesh.h"
 


void  distribution(int PartitionGrid , int st , int stp);

 void  distribution(int PartitionGrid , int st , int stp){


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

      w.fwx   = new double*[PartitionGrid];
      w.fwy   = new double*[PartitionGrid];
      w.afx   = new double*[PartitionGrid];
      w.afy   = new double*[PartitionGrid];
 
      S.s1   = new double*[PartitionGrid];
      S.s2   = new double*[PartitionGrid];
      S.s3   = new double*[PartitionGrid];
      S.s4   = new double*[PartitionGrid];
      S.s5   = new double*[PartitionGrid];
      S.s6   = new double*[PartitionGrid];
      S.s7   = new double*[PartitionGrid];
      S.s8   = new double*[PartitionGrid];
 
      k.u1   = new double*[PartitionGrid];
      k.u2   = new double*[PartitionGrid];
      k.u3   = new double*[PartitionGrid];
      k.u4   = new double*[PartitionGrid];
      k.u5   = new double*[PartitionGrid];
      k.u6   = new double*[PartitionGrid];
      k.u7   = new double*[PartitionGrid];
      k.u8   = new double*[PartitionGrid];

      k2.u1   = new double*[PartitionGrid];
      k2.u2   = new double*[PartitionGrid];
      k2.u3   = new double*[PartitionGrid];
      k2.u4   = new double*[PartitionGrid];
      k2.u5   = new double*[PartitionGrid];
      k2.u6   = new double*[PartitionGrid];
      k2.u7   = new double*[PartitionGrid];
      k2.u8   = new double*[PartitionGrid];

      u.U1   = new double*[PartitionGrid];
      u.U2   = new double*[PartitionGrid];
      u.U3   = new double*[PartitionGrid];
      u.U4   = new double*[PartitionGrid];
      u.U5   = new double*[PartitionGrid];
      u.U6   = new double*[PartitionGrid];
      u.U7   = new double*[PartitionGrid];
      u.U8   = new double*[PartitionGrid];

      f.F1   = new double*[PartitionGrid];
      f.F2   = new double*[PartitionGrid];
      f.F3   = new double*[PartitionGrid];
      f.F4   = new double*[PartitionGrid];
      f.F5   = new double*[PartitionGrid];
      f.F6   = new double*[PartitionGrid];
      f.F7   = new double*[PartitionGrid];
      f.F8   = new double*[PartitionGrid];

      f.G1   = new double*[PartitionGrid];
      f.G2   = new double*[PartitionGrid];
      f.G3   = new double*[PartitionGrid];
      f.G4   = new double*[PartitionGrid];
      f.G5   = new double*[PartitionGrid];
      f.G6   = new double*[PartitionGrid];
      f.G7   = new double*[PartitionGrid];
      f.G8   = new double*[PartitionGrid];

      c.lx1   = new double*[PartitionGrid];
      c.lx2   = new double*[PartitionGrid];
      c.ly1   = new double*[PartitionGrid];
      c.ly2   = new double*[PartitionGrid];
      c.lMaxx = new double*[PartitionGrid];
      c.lMaxy = new double*[PartitionGrid];
      c.lx    = new double*[PartitionGrid];
      c.ly    = new double*[PartitionGrid];


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

    // waves w
	    w.fwx[i] = new double[xpoints];
	    w.fwy[i] = new double[xpoints];
	    w.afx[i] = new double[xpoints];
	    w.afy[i] = new double[xpoints];

    // Powell S
	    S.s1[i] = new double[xpoints];
	    S.s2[i] = new double[xpoints];
	    S.s3[i] = new double[xpoints];
	    S.s4[i] = new double[xpoints];
	    S.s5[i] = new double[xpoints];
	    S.s6[i] = new double[xpoints];
	    S.s7[i] = new double[xpoints];
	    S.s8[i] = new double[xpoints];

    // RungeKutta k
	    k.u1[i] = new double[xpoints];
	    k.u2[i] = new double[xpoints];
	    k.u3[i] = new double[xpoints];
	    k.u4[i] = new double[xpoints];
	    k.u5[i] = new double[xpoints];
	    k.u6[i] = new double[xpoints];
	    k.u7[i] = new double[xpoints];
	    k.u8[i] = new double[xpoints];

    // RungeKutta k2
	    k2.u1[i] = new double[xpoints];
	    k2.u2[i] = new double[xpoints];
	    k2.u3[i] = new double[xpoints];
	    k2.u4[i] = new double[xpoints];
	    k2.u5[i] = new double[xpoints];
	    k2.u6[i] = new double[xpoints];
	    k2.u7[i] = new double[xpoints];
	    k2.u8[i] = new double[xpoints];

    // conservative u
	    u.U1[i] = new double[xpoints];
	    u.U2[i] = new double[xpoints];
	    u.U3[i] = new double[xpoints];
	    u.U4[i] = new double[xpoints];
	    u.U5[i] = new double[xpoints];
	    u.U6[i] = new double[xpoints];
	    u.U7[i] = new double[xpoints];
            u.U8[i] = new double[xpoints];

    // flux f
	    f.F1[i] = new double[xpoints];
	    f.F2[i] = new double[xpoints];
	    f.F3[i] = new double[xpoints];
	    f.F4[i] = new double[xpoints];
	    f.F5[i] = new double[xpoints];
 	    f.F6[i] = new double[xpoints];
            f.F7[i] = new double[xpoints];
            f.F8[i] = new double[xpoints];

            f.G1[i] = new double[xpoints];
	    f.G2[i] = new double[xpoints];
	    f.G3[i] = new double[xpoints];
	    f.G4[i] = new double[xpoints];
	    f.G5[i] = new double[xpoints];
	    f.G6[i] = new double[xpoints];
	    f.G7[i] = new double[xpoints];
	    f.G8[i] = new double[xpoints];

    // cellinfo c
	    c.lx1[i]   = new double[xpoints];
	    c.lx2[i]   = new double[xpoints];	
	    c.ly1[i]   = new double[xpoints];
	    c.ly2[i]   = new double[xpoints];
	    c.lMaxx[i] = new double[xpoints];
	    c.lMaxy[i] = new double[xpoints];
	    c.lx[i]    = new double[xpoints];
	    c.ly[i]    = new double[xpoints];



            }

        
       }

#endif

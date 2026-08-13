#ifndef DC3DQ_H
#define DC3DQ_H

#include <stdio.h>
#include <math.h>
#include <quadmath.h>

typedef struct
{
  __float128 alp1, alp2, alp3, alp4, alp5;
  __float128 sd;
  __float128 cd; 
  __float128 sdsd, cdcd, sdcd, s2d, c2d;
} c0_str_q;

typedef struct
{
  __float128 p,q,s,t,xy,x2,y2,d2,r,r2,r3,r5,qr,qrx,a3,a5,b3,c3;  
  __float128 uy,vy,wy,uz,vz,wz;
} c1_str_q;

typedef struct 
{
  __float128 xi2,et2,q2,r,r2,r3,r5,y,d,tt,alx,ale,x11,y11,x32,y32,  
         ey,ez,fy,fz,gy,gz,hy,hz;
} c2_str_q;                                


int dc3dq(__float128 alpha, __float128 x, __float128 y, __float128 z,
         __float128 depth, __float128 dip, __float128 al1, __float128 al2,
         __float128 aw1, __float128 aw2, __float128 disl1, __float128 disl2,
         __float128 disl3, __float128 *ux, __float128 *uy, __float128 *uz,
         __float128 *uxx, __float128 *uyx, __float128 *uzx, __float128 *uxy,
         __float128 *uyy, __float128 *uzy, __float128 *uxz, __float128 *uyz,
         __float128 *uzz);
void uaq(__float128 xi, __float128 et, __float128 q, __float128 disl1, __float128 disl2, __float128 disl3, 
        __float128 *u, c0_str_q c0, c1_str_q c1, c2_str_q c2);

void ubq(__float128 xi, __float128 et, __float128 q, __float128 disl1, __float128 disl2, __float128 disl3, 
        __float128 *u, c0_str_q c0, c1_str_q c1, c2_str_q c2);
void ucq(__float128 xi, __float128 et, __float128 q, __float128 z, __float128 disl1, __float128 disl2, __float128 disl3, 
        __float128 *u,  c0_str_q c0, c1_str_q c1, c2_str_q c2);
void dccon0q(__float128 alpha, __float128 dip,  c0_str_q *c0, c1_str_q *c1, c2_str_q *c2);  
void dccon1q(__float128 x, __float128 y, __float128 d,  c0_str_q *c0, c1_str_q *c1, c2_str_q *c2);
void dccon2q(__float128 xi, __float128 et, __float128 q, __float128 sd, __float128 cd, __float128 kxi, 
           __float128 ket,  c0_str_q *c0, c1_str_q *c1, c2_str_q *c2);
void dccon3q(__float128 x, __float128 p, __float128 q, __float128 al1, __float128 al2, __float128 aw1, 
            __float128 aw2, __float128 *kxi, __float128 *ket);
  
int dc3d0q(__float128 alpha, __float128 x, __float128 y, __float128 z, __float128 depth, __float128 dip,
          __float128 pot1, __float128 pot2, __float128 pot3, __float128 pot4,      
          __float128 *ux, __float128 *uy, __float128 *uz, 
          __float128 *uxx, __float128 *uyx, __float128 *uzx,
          __float128 *uxy, __float128 *uyy, __float128 *uzy,
          __float128 *uxz, __float128 *uyz, __float128 *uzz) ;          
void ua0q(__float128 x, __float128 y, __float128 d, __float128 pot1, __float128 pot2, 
         __float128 pot3, __float128 pot4, __float128 *u, c0_str_q *c0, c1_str_q *c1);
void ub0q(__float128 x, __float128 y, __float128 d, __float128 z,
	     __float128 pot1, __float128 pot2, __float128 pot3, __float128 pot4,
	     __float128 *u, c0_str_q *c0, c1_str_q *c1);           
void uc0q(__float128 x, __float128 y, __float128 d, __float128 z, 
         __float128 pot1, __float128 pot2, __float128 pot3, __float128 pot4, 
         __float128 u[12], c0_str_q *c0, c1_str_q *c1);

void dc3d_passByRefq(__float128 *alpha_p, __float128 *x_p, __float128 *y_p, __float128 *z_p,
          __float128 *depth_p, __float128 *dip_p, __float128 *al1_p, __float128 *al2_p,
          __float128 *aw1_p, __float128 *aw2_p, __float128 *disl1_p, __float128 *disl2_p,
          __float128 *disl3_p, __float128 *ux, __float128 *uy, __float128 *uz,
          __float128 *uxx, __float128 *uyx, __float128 *uzx, __float128 *uxy,
          __float128 *uyy, __float128 *uzy, __float128 *uxz, __float128 *uyz,
          __float128 *uzz);

void dc3d_passByRef_vq(
           __float128 *alpha, __float128 *x, __float128 *y, __float128 *z,
           __float128 *depth, __float128 *dip, __float128 *al1, __float128 *al2,
           __float128 *aw1, __float128 *aw2, __float128 *disl1, __float128 *disl2,
           __float128 *disl3, __float128 *ux, __float128 *uy, __float128 *uz,
           __float128 *uxx, __float128 *uyx, __float128 *uzx, __float128 *uxy,
           __float128 *uyy, __float128 *uzy, __float128 *uxz, __float128 *uyz,
           __float128 *uzz, int *npts);

#endif /* ndef DC3DQ_H */

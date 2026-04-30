#ifndef DC3D_H
#define DC3D_H

#include <stdio.h>
#include <math.h>

typedef struct
{
  double alp1, alp2, alp3, alp4, alp5;
  double sd;
  double cd; 
  double sdsd, cdcd, sdcd, s2d, c2d;
} c0_str;

typedef struct
{
  double p,q,s,t,xy,x2,y2,d2,r,r2,r3,r5,qr,qrx,a3,a5,b3,c3;  
  double uy,vy,wy,uz,vz,wz;
} c1_str;

typedef struct 
{
  double xi2,et2,q2,r,r2,r3,r5,y,d,tt,alx,ale,x11,y11,x32,y32,  
         ey,ez,fy,fz,gy,gz,hy,hz;
} c2_str;                                


int dc3d(double alpha, double x, double y, double z,
         double depth, double dip, double al1, double al2,
         double aw1, double aw2, double disl1, double disl2,
         double disl3, double *ux, double *uy, double *uz,
         double *uxx, double *uyx, double *uzx, double *uxy,
         double *uyy, double *uzy, double *uxz, double *uyz,
         double *uzz);
void ua(double xi, double et, double q, double disl1, double disl2, double disl3, 
        double *u, c0_str c0, c1_str c1, c2_str c2);

void ub(double xi, double et, double q, double disl1, double disl2, double disl3, 
        double *u, c0_str c0, c1_str c1, c2_str c2);
void uc(double xi, double et, double q, double z, double disl1, double disl2, double disl3, 
        double *u,  c0_str c0, c1_str c1, c2_str c2);
void dccon0(double alpha, double dip,  c0_str *c0, c1_str *c1, c2_str *c2);  
void dccon1(double x, double y, double d,  c0_str *c0, c1_str *c1, c2_str *c2);
void dccon2(double xi, double et, double q, double sd, double cd, double kxi, 
           double ket,  c0_str *c0, c1_str *c1, c2_str *c2);
void dccon3(double x, double p, double q, double al1, double al2, double aw1, 
            double aw2, double *kxi, double *ket);
  
int dc3d0(double alpha, double x, double y, double z, double depth, double dip,
          double pot1, double pot2, double pot3, double pot4,      
          double *ux, double *uy, double *uz, 
          double *uxx, double *uyx, double *uzx,
          double *uxy, double *uyy, double *uzy,
          double *uxz, double *uyz, double *uzz) ;          
void ua0(double x, double y, double d, double pot1, double pot2, 
         double pot3, double pot4, double *u, c0_str *c0, c1_str *c1);
void ub0(double x, double y, double d, double z,
	     double pot1, double pot2, double pot3, double pot4,
	     double *u, c0_str *c0, c1_str *c1);           
void uc0(double x, double y, double d, double z, 
         double pot1, double pot2, double pot3, double pot4, 
         double u[12], c0_str *c0, c1_str *c1);

void dc3d_passByRef(double *alpha_p, double *x_p, double *y_p, double *z_p,
          double *depth_p, double *dip_p, double *al1_p, double *al2_p,
          double *aw1_p, double *aw2_p, double *disl1_p, double *disl2_p,
          double *disl3_p, double *ux, double *uy, double *uz,
          double *uxx, double *uyx, double *uzx, double *uxy,
          double *uyy, double *uzy, double *uxz, double *uyz,
          double *uzz);

void dc3d_passByRef_v(
           double *alpha, double *x, double *y, double *z,
           double *depth, double *dip, double *al1, double *al2,
           double *aw1, double *aw2, double *disl1, double *disl2,
           double *disl3, double *ux, double *uy, double *uz,
           double *uxx, double *uyx, double *uzx, double *uxy,
           double *uyy, double *uzy, double *uxz, double *uyz,
           double *uzz, int *npts);

#endif /* ndef DC3D_H */

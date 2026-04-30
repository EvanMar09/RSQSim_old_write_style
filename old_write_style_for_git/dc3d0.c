/* nearly line-by-line translation of Okada [1992]'s Fortran into C
   Keith Ricards-Dinger, 2013
*/

#include "dc3d.h"

/*
C                                                                       00060000
C********************************************************************   00070000
C*****                                                          *****   00080000
C*****    DISPLACEMENT AND STRAIN AT DEPTH                      *****   00090000
C*****    DUE TO BURIED POINT SOURCE IN A SEMIINFINITE MEDIUM   *****   00100000
C*****                         CODED BY  Y.OKADA ... SEP.1991   *****   00110002
C*****                         REVISED   Y.OKADA ... NOV.1991   *****   00120002
C*****                                                          *****   00130000
C********************************************************************   00140000
C                                                                       00150000
C***** INPUT                                                            00160000
C*****   ALPHA : MEDIUM CONSTANT  (LAMBDA+MYU)/(LAMBDA+2*MYU)           00170000
C*****   X,Y,Z : COORDINATE OF OBSERVING POINT                          00180000
C*****   DEPTH : SOURCE DEPTH                                           00190000
C*****   DIP   : DIP-ANGLE (DEGREE)                                     00200000
C*****   POT1-POT4 : STRIKE-, DIP-, TENSILE- AND INFLATE-POTENCY        00210000
C*****       POTENCY=(  MOMENT OF DOUBLE-COUPLE  )/MYU     FOR POT1,2   00220000
C*****       POTENCY=(INTENSITY OF ISOTROPIC PART)/LAMBDA  FOR POT3     00230000
C*****       POTENCY=(INTENSITY OF LINEAR DIPOLE )/MYU     FOR POT4     00240000
C                                                                       00250000
C***** OUTPUT                                                           00260000
C*****   UX, UY, UZ  : DISPLACEMENT ( UNIT=(UNIT OF POTENCY) /          00270000
C*****               :                     (UNIT OF X,Y,Z,DEPTH)**2  )  00280000
C*****   UXX,UYX,UZX : X-DERIVATIVE ( UNIT= UNIT OF POTENCY) /          00290000
C*****   UXY,UYY,UZY : Y-DERIVATIVE        (UNIT OF X,Y,Z,DEPTH)**3  )  00300000
C*****   UXZ,UYZ,UZZ : Z-DERIVATIVE                                     00310000
C*****   IRET        : RETURN CODE  ( =0....NORMAL,   =1....SINGULAR )  00320002
C                                                                       00330000
*/
int dc3d0(double alpha, double x, double y, double z, double depth, double dip,
          double pot1, double pot2, double pot3, double pot4,      
          double *ux, double *uy, double *uz, 
          double *uxx, double *uyx, double *uzx,
          double *uxy, double *uyy, double *uzy,
          double *uxz, double *uyz, double *uzz) 
{
  c0_str c0;
  c1_str c1;
  c2_str c2;
  
  int i;

  double  u[12], dua[12], dub[12], duc[12]; 
  double dd, du;
  
  /* value returned if singular */
  *ux = *uy = *uz = *uxx = *uxy = *uxz = *uyx = *uyy = *uyz = 
        *uzx = *uzy = *uzz = 0.0;  

  if (z > 0.0) fprintf(stderr, "POSITIVE Z WAS GIVEN IN DC3D0\n");

  for (i=0; i<12; i++)
  {
	u[i] = 0.0;                                                         
	dua[i] = 0.0;                                                       
	dub[i] = 0.0;                                                       
	duc[i] = 0.0;                                                       
  }

  dccon0(alpha, dip, &c0, &c1, &c2);                              

/*
C======================================                                 
C=====  REAL-SOURCE CONTRIBUTION  =====                                 
C====================================== 
*/                                
  
  dd = depth + z;
  dccon1(x, y, dd, &c0, &c1, &c2); 
											  
  if (c1.r == 0.0) return(1);  
		
  ua0(x, y, dd, pot1, pot2, pot3, pot4, dua, &c0, &c1);                         

  for (i=0; i<12; i++)
  {                                                    
	if (i < 9) u[i] = u[i] - dua[i];                                    
	if (i >= 9) u[i] = u[i] + dua[i]; 
  }                                   

/*
C=======================================                                
C=====  IMAGE-SOURCE CONTRIBUTION  =====                                
C=======================================                                
*/

  dd = depth - z;                                                        
  dccon1(x, y, dd, &c0, &c1, &c2);                                             
  ua0(x, y, dd, pot1, pot2, pot3, pot4, dua, &c0, &c1);                            
  ub0(x, y, dd, z, pot1, pot2, pot3, pot4, dub, &c0, &c1);                         
  uc0(x, y, dd, z, pot1, pot2, pot3, pot4, duc, &c0, &c1);                         

  for (i=0; i<12; i++)
  {
        du = dua[i] + dub[i] + z*duc[i];                                      
        if (i >= 9) du = du + duc[i-9];                                     
        u[i] = u[i] + du;                                                    
  }
  
  
  *ux = u[0];                                                           
  *uy = u[1];                                                           
  *uz = u[2];                                                           
  *uxx = u[3];                                                          
  *uyx = u[4];                                                          
  *uzx = u[5];                                                          
  *uxy = u[6];                                                          
  *uyy = u[7];                                                          
  *uzy = u[8];                                                          
  *uxz = u[9];                                                         
  *uyz = u[10];                                                         
  *uzz = u[11]; 
                                                                                                                    
  return(0);                                                          
}

/*
C                                                                       01170000
C********************************************************************   01180000
C*****    DISPLACEMENT AND STRAIN AT DEPTH (PART-A)             *****   01190000
C*****    DUE TO BURIED POINT SOURCE IN A SEMIINFINITE MEDIUM   *****   01200000
C********************************************************************   01210000
C                                                                       01220000
C***** INPUT                                                            01230000
C*****   X,Y,D : STATION COORDINATES IN FAULT SYSTEM                    01240000
C*****   POT1-POT4 : STRIKE-, DIP-, TENSILE- AND INFLATE-POTENCY        01250000
C***** OUTPUT                                                           01260000
C*****   U(12) : DISPLACEMENT AND THEIR DERIVATIVES                     01270000
C                                                                       01280000
*/

void ua0(double x, double y, double d, double pot1, double pot2, 
         double pot3, double pot4, double *u, c0_str *c0, c1_str *c1)
{                      
  double du[12];                                            

  double pi2 = 6.283185307179586; 
  
  int i;   
  
  for (i=0; i<12; i++) u[i] = 0.0;                                                    

/*
C======================================                                 
C=====  STRIKE-SLIP CONTRIBUTION  =====                                 
C====================================== 
*/                                

  if (pot1 != 0.0) 
  {                                             
	du[ 0] = c0->alp1*c1->q/c1->r3    + c0->alp2*c1->x2*c1->qr;                                
	du[ 1] = c0->alp1*x/c1->r3*c0->sd + c0->alp2*c1->xy*c1->qr;                              
	du[ 2] =-c0->alp1*x/c1->r3*c0->cd + c0->alp2*x*d*c1->qr;                               
	du[ 3] = x*c1->qr*(-c0->alp1 + c0->alp2*(1.0 + c1->a5) );                             
	du[ 4] = c0->alp1*c1->a3/c1->r3*c0->sd + c0->alp2*y*c1->qr*c1->a5;                             
	du[ 5] =-c0->alp1*c1->a3/c1->r3*c0->cd + c0->alp2*d*c1->qr*c1->a5;                             
	du[ 6] = c0->alp1*(c0->sd/c1->r3 - y*c1->qr) + c0->alp2*3.0*c1->x2/c1->r5*c1->uy;                     
	du[ 7] = 3.0*x/c1->r5*(-c0->alp1*y*c0->sd + c0->alp2*(y*c1->uy + c1->q) );                    
	du[ 8] = 3.0*x/c1->r5*( c0->alp1*y*c0->cd + c0->alp2*d*c1->uy );                        
	du[ 9] = c0->alp1*(c0->cd/c1->r3 + d*c1->qr) + c0->alp2*3.0*c1->x2/c1->r5*c1->uz;                     
	du[10] = 3.0*x/c1->r5*( c0->alp1*d*c0->sd + c0->alp2*y*c1->uz );                        
	du[11] = 3.0*x/c1->r5*(-c0->alp1*d*c0->cd + c0->alp2*(d*c1->uz - c1->q) );   
	                 
	for (i=0; i<12; i++) u[i] += pot1/pi2*du[i];                                        
  }                                                           
/*
C===================================                                    
C=====  DIP-SLIP CONTRIBUTION  =====                                    
C===================================                                    
*/

  if (pot2 != 0.0) 
  {                                              
	du[ 0] =                         c0->alp2*x*c1->p*c1->qr;                                  
	du[ 1] = c0->alp1*c1->s/c1->r3 + c0->alp2*y*c1->p*c1->qr;                                  
	du[ 2] =-c0->alp1*c1->t/c1->r3 + c0->alp2*d*c1->p*c1->qr;                                  
	du[ 3] =                         c0->alp2*c1->p*c1->qr*c1->a5;                            
	du[ 4] =-c0->alp1*3.0*x*c1->s/c1->r5 - c0->alp2*y*c1->p*c1->qrx;                            
	du[ 5] = c0->alp1*3.0*x*c1->t/c1->r5 - c0->alp2*d*c1->p*c1->qrx;                            
	du[ 6] =                          c0->alp2*3.0*x/c1->r5*c1->vy;                
	du[ 7] = c0->alp1*(c0->s2d/c1->r3 - 3.0*y*c1->s/c1->r5) + 
			 c0->alp2*(3.0*y/c1->r5*c1->vy + c1->p*c1->qr);         
	du[ 8] =-c0->alp1*(c0->c2d/c1->r3-3.0*y*c1->t/c1->r5) + c0->alp2*3.0*d/c1->r5*c1->vy;                
	du[ 9] =                          c0->alp2*3.0*x/c1->r5*c1->vz;                
	du[10] = c0->alp1*(c0->c2d/c1->r3 + 3.0*d*c1->s/c1->r5) + 
	         c0->alp2*3.0*y/c1->r5*c1->vz;                
	du[11] = c0->alp1*(c0->s2d/c1->r3 - 3.0*d*c1->t/c1->r5) +
	         c0->alp2*(3.0*d/c1->r5*c1->vz - c1->p*c1->qr);         

	for (i=0; i<12; i++) u[i] = u[i] + pot2/pi2*du[i];                                        

  }                                                             

/*
C========================================                               
C=====  TENSILE-FAULT CONTRIBUTION  =====                               
C========================================                               
*/

  if (pot3 != 0.0) 
  {                                              
	du[ 0] = c0->alp1*x/c1->r3 - c0->alp2*x*c1->q*c1->qr;                                 
	du[ 1] = c0->alp1*c1->t/c1->r3 - c0->alp2*y*c1->q*c1->qr;                                  
	du[ 2] = c0->alp1*c1->s/c1->r3 - c0->alp2*d*c1->q*c1->qr;                                  
	du[ 3] = c0->alp1*c1->a3/c1->r3 - c0->alp2*c1->q*c1->qr*c1->a5;                            
	du[ 4] =-c0->alp1*3.0*x*c1->t/c1->r5 + c0->alp2*y*c1->q*c1->qrx;                           
	du[ 5] =-c0->alp1*3.0*x*c1->s/c1->r5 + c0->alp2*d*c1->q*c1->qrx;                            
	du[ 6] =-c0->alp1*3.0*c1->xy/c1->r5  - c0->alp2*x*c1->qr*c1->wy;                   
	du[ 7] = c0->alp1*(c0->c2d/c1->r3-3.0*y*c1->t/c1->r5) - 
	         c0->alp2*(y*c1->wy + c1->q)*c1->qr;             
	du[ 8] = c0->alp1*(c0->s2d/c1->r3 - 3.0*y*c1->s/c1->r5) - c0->alp2*d*c1->qr*c1->wy;                   
	du[ 9] = c0->alp1*3.0*x*d/c1->r5 - c0->alp2*x*c1->qr*c1->wz;                   
	du[10] =-c0->alp1*(c0->s2d/c1->r3 - 3.0*d*c1->t/c1->r5) - c0->alp2*y*c1->qr*c1->wz;                   
	du[11] = c0->alp1*(c0->c2d/c1->r3+3.0*d*c1->s/c1->r5) - 
	         c0->alp2*(d*c1->wz - c1->q)*c1->qr;              

	for (i=0; i<12; i++)  u[i] = u[i] + pot3/pi2*du[i];                                        
  }                                                             

/*
C=========================================                              
C=====  INFLATE SOURCE CONTRIBUTION  =====                              
C========================================= 
*/                             

  if (pot4 != 0.0)
  {                                               
	du[ 0] =-c0->alp1*x/c1->r3;                                              
	du[ 1] =-c0->alp1*y/c1->r3;                                              
	du[ 2] =-c0->alp1*d/c1->r3;                                               
	du[ 3] =-c0->alp1*c1->a3/c1->r3;                                              
	du[ 4] = c0->alp1*3.0*c1->xy/c1->r5;                                           
	du[ 5] = c0->alp1*3.0*x*d/c1->r5;                                          
	du[ 6] = du[4];                                        
	du[ 7] =-c0->alp1*c1->b3/c1->r3;                                              
	du[ 8] = c0->alp1*3.0*y*d/c1->r5;                                          
	du[ 9] =-du[5];                                                   
	du[10] =-du[8];                                                  
	du[11] = c0->alp1*c1->c3/c1->r3;                                              

	for (i=0; i<12; i++) u[i] = u[i] + pot4/pi2*du[i];                                        

  }                                                             

  return; 
}                                                           


/*
C                                                                       02180000
C********************************************************************   02190000
C*****    DISPLACEMENT AND STRAIN AT DEPTH (PART-B)             *****   02200000
C*****    DUE TO BURIED POINT SOURCE IN A SEMIINFINITE MEDIUM   *****   02210000
C********************************************************************   02220000
C                                                                       02230000
C***** INPUT                                                            02240000
C*****   X,Y,D,Z : STATION COORDINATES IN FAULT SYSTEM                  02250000
C*****   POT1-POT4 : STRIKE-, DIP-, TENSILE- AND INFLATE-POTENCY        02260000
C***** OUTPUT                                                           02270000
C*****   U(12) : DISPLACEMENT AND THEIR DERIVATIVES                     02280000
C                                                                       02290000
*/

      void ub0(double x, double y, double d, double z,
               double pot1, double pot2, double pot3, double pot4,
               double *u, c0_str *c0, c1_str *c1)
{                   

      double du[12];
      double c, rd, d12, d32, d33, d53, d54;   
      double fi1, fi2, fi3, fi4, fi5, fj1, fj2, fj3, fj4, fk1, fk2, fk3;                                        

      double pi2 = 6.283185307179586;  
      
      int i;                                  

      c = d + z;                                                           
      rd = c1->r + d;                                                        
      d12 = 1.0/(c1->r*rd*rd);                                                
      d32 = d12*(2.0*c1->r + d)/c1->r2;                                            
      d33 = d12*(3.0*c1->r + d)/(c1->r2*rd);                                         
      d53 = d12*(8.0*c1->r2 + 9.0*c1->r*d + 3.0*c1->d2)/(c1->r2*c1->r2*rd);                          
      d54 = d12*(5.0*c1->r2 + 4.0*c1->r*d + c1->d2)/c1->r3*d12;                                 

      fi1 = y*(d12 - c1->x2*d33);                                              
      fi2 = x*(d12 - c1->y2*d33);                                              
      fi3 = x/c1->r3 - fi2;                                                    
      fi4 = -c1->xy*d32;                                                      
      fi5 = 1.0/(c1->r*rd) - c1->x2*d32;                                            
      fj1 = -3.0*c1->xy*(d33 - c1->x2*d54);                                          
      fj2 = 1.0/c1->r3 - 3.0*d12 + 3.0*c1->x2*c1->y2*d54;                                   
      fj3 = c1->a3/c1->r3 - fj2;                                                   
      fj4 = -3.0*c1->xy/c1->r5 - fj1;                                                
      fk1 = -y*(d32 - c1->x2*d53);                                              
      fk2 = -x*(d32 - c1->y2*d53);                                              
      fk3 = -3.0*x*d/c1->r5 - fk2;                                               

      for (i=0; i<12; i++) u[i] = 0.0;                                                          

/*
C======================================                                
C=====  STRIKE-SLIP CONTRIBUTION  =====                                
C======================================                                
*/
  if(pot1 != 0.0)
  {                                              
	du[ 0] = -c1->x2*c1->qr - c0->alp3*fi1*c0->sd;                                    
	du[ 1] = -c1->xy*c1->qr - c0->alp3*fi2*c0->sd;                                    
	du[ 2] = -c*x*c1->qr - c0->alp3*fi4*c0->sd;                                    
	du[ 3] = -x*c1->qr*(1.0 + c1->a5) - c0->alp3*fj1*c0->sd;                              
	du[ 4] = -y*c1->qr*c1->a5 - c0->alp3*fj2*c0->sd;                              
	du[ 5] = -c*c1->qr*c1->a5 - c0->alp3*fk1*c0->sd;                              
	du[ 6] = -3.0*c1->x2/c1->r5*c1->uy - c0->alp3*fj2*c0->sd;                          
	du[ 7] = -3.0*c1->xy/c1->r5*c1->uy - x*c1->qr - c0->alp3*fj4*c0->sd;                          
	du[ 8] = -3.0*c*x/c1->r5*c1->uy - c0->alp3*fk2*c0->sd;                          
	du[ 9] = -3.0*c1->x2/c1->r5*c1->uz + c0->alp3*fk1*c0->sd;                              
	du[10] = -3.0*c1->xy/c1->r5*c1->uz + c0->alp3*fk2*c0->sd;                              
	du[11] = 3.0*x/c1->r5*(-c*c1->uz + c0->alp3*y*c0->sd);                             

	for (i=0; i<12; i++) u[i] = u[i] + pot1/pi2*du[i];   
  }                                                         

/*
C===================================                                   
C=====  DIP-SLIP CONTRIBUTION  =====                                   
C===================================                                   
*/

  if (pot2 != 0.0) 
  {                                           
	du[ 0] = -x*c1->p*c1->qr + c0->alp3*fi3*c0->sdcd;                                  
	du[ 1] = -y*c1->p*c1->qr + c0->alp3*fi1*c0->sdcd;                                  
	du[ 2] = -c*c1->p*c1->qr + c0->alp3*fi5*c0->sdcd;                                  
	du[ 3] = -c1->p*c1->qr*c1->a5 + c0->alp3*fj3*c0->sdcd;                                 
	du[ 4] = y*c1->p*c1->qrx + c0->alp3*fj1*c0->sdcd;                                 
	du[ 5] = c*c1->p*c1->qrx + c0->alp3*fk3*c0->sdcd;                                 
	du[ 6] = -3.0*x/c1->r5*c1->vy + c0->alp3*fj1*c0->sdcd;                         
	du[ 7] = -3.0*y/c1->r5*c1->vy - c1->p*c1->qr + c0->alp3*fj2*c0->sdcd;                         
	du[ 8] = -3.0*c/c1->r5*c1->vy + c0->alp3*fk1*c0->sdcd;                         
	du[ 9] = -3.0*x/c1->r5*c1->vz - c0->alp3*fk3*c0->sdcd;                              
	du[10] = -3.0*y/c1->r5*c1->vz - c0->alp3*fk1*c0->sdcd;                              
	du[11] = -3.0*c/c1->r5*c1->vz + c0->alp3*c1->a3/c1->r3*c0->sdcd;                            

	for (i=0; i<12; i++) u[i] = u[i] + pot2/pi2*du[i];                                       
  }                                                          

/*
C========================================                              
C=====  TENSILE-FAULT CONTRIBUTION  =====                              
C========================================                              
*/

  if (pot3 != 0.0)
  {                                              
	du[ 0] = x*c1->q*c1->qr - c0->alp3*fi3*c0->sdsd;                                  
	du[ 1] = y*c1->q*c1->qr - c0->alp3*fi1*c0->sdsd;                                  
	du[ 2] = c*c1->q*c1->qr - c0->alp3*fi5*c0->sdsd;                                  
	du[ 3] = c1->q*c1->qr*c1->a5 - c0->alp3*fj3*c0->sdsd;                                 
	du[ 4] =-y*c1->q*c1->qrx - c0->alp3*fj1*c0->sdsd;                                 
	du[ 5] =-c*c1->q*c1->qrx - c0->alp3*fk3*c0->sdsd;                                 
	du[ 6] = x*c1->qr*c1->wy - c0->alp3*fj1*c0->sdsd;                             
	du[ 7] = c1->qr*(y*c1->wy + c1->q) - c0->alp3*fj2*c0->sdsd;                             
	du[ 8] = c*c1->qr*c1->wy - c0->alp3*fk1*c0->sdsd;                             
	du[ 9] = x*c1->qr*c1->wz + c0->alp3*fk3*c0->sdsd;                                 
	du[10] = y*c1->qr*c1->wz + c0->alp3*fk1*c0->sdsd;                                 
	du[11] = c*c1->qr*c1->wz - c0->alp3*c1->a3/c1->r3*c0->sdsd;                               

	for (i=0; i<12; i++) u[i] = u[i] + pot3/pi2*du[i];
  }                                                         

/*
C=========================================                             
C=====  INFLATE SOURCE CONTRIBUTION  =====                             
C=========================================                             
*/

  if (pot4 != 0.0)
  {                                              
	du[ 0] = c0->alp3*x/c1->r3;                                              
	du[ 1] = c0->alp3*y/c1->r3;                                              
	du[ 2] = c0->alp3*d/c1->r3;                                              
	du[ 3] = c0->alp3*c1->a3/c1->r3;                                             
	du[ 4] = -c0->alp3*3.0*c1->xy/c1->r5;                                          
	du[ 5] = -c0->alp3*3.0*x*d/c1->r5;                                         
	du[ 6] = du[4];                                                  
	du[ 7] = c0->alp3*c1->b3/c1->r3;                                             
	du[ 8] = -c0->alp3*3.0*y*d/c1->r5;                                         
	du[ 9] = -du[5];                                                  
	du[10] = -du[8];                                                  
	du[11] = -c0->alp3*c1->c3/c1->r3;                                             

	for (i=0; i<12; i++) u[i] = u[i] + pot4/pi2*du[i];  
  }                                                         

  return; 
}                                                          

/*
C********************************************************************   
C*****    DISPLACEMENT AND STRAIN AT DEPTH (PART-B)             *****   
C*****    DUE TO BURIED POINT SOURCE IN A SEMIINFINITE MEDIUM   *****   
C********************************************************************   
C                                                                       
C***** INPUT                                                            
C*****   X,Y,D,Z : STATION COORDINATES IN FAULT SYSTEM                  
C*****   POT1-POT4 : STRIKE-, DIP-, TENSILE- AND INFLATE-POTENCY        
C***** OUTPUT                                                           
C*****   U(12) : DISPLACEMENT AND THEIR DERIVATIVES                                                                                            
*/

  void uc0(double x, double y, double d, double z, 
           double pot1, double pot2, double pot3, double pot4, 
           double u[12], c0_str *c0, c1_str *c1)                    
{
 
  double du[12];                                            
  
  double pi2 = 6.283185307179586;
  
  double c, q2, r7, a7, b5, b7, c5, c7, d7, qr5, qr7, dr5;   
  
  int i;                                  

  c = d + z;                                                           
  q2 = c1->q*c1->q;                                                          
  r7 = c1->r5*c1->r2;                                                          
  a7 = 1.0 - 7.0*c1->x2/c1->r2;                                                    
  b5 = 1.0 - 5.0*c1->y2/c1->r2;                                                    
  b7 = 1.0 - 7.0*c1->y2/c1->r2;                                                    
  c5 = 1.0 - 5.0*c1->d2/c1->r2;                                                    
  c7 = 1.0 - 7.0*c1->d2/c1->r2;                                                    
  d7 = 2.0 - 7.0*q2/c1->r2;                                                    
  qr5 = 5.0*c1->q/c1->r2;                                                       
  qr7 = 7.0*c1->q/c1->r2;                                                       
  dr5 = 5.0*d/c1->r2;                                                       


  for (i=0; i<12; i++) u[i] = 0.0;                                                   

/*
C======================================                                 
C=====  STRIKE-SLIP CONTRIBUTION  =====                                 
C======================================                                 
*/

  if (pot1 != 0.0)
  {                                               
	du[ 0] = -c0->alp4*c1->a3/c1->r3*c0->cd  + c0->alp5*c*c1->qr*c1->a5;                            
	du[ 1] = 3.0*x/c1->r5*(c0->alp4*y*c0->cd + c0->alp5*c*(c0->sd - y*qr5) );                
	du[ 2] = 3.0*x/c1->r5*(-c0->alp4*y*c0->sd + c0->alp5*c*(c0->cd + d*qr5) );                
	du[ 3] = c0->alp4*3.0*x/c1->r5*(2.0 + c1->a5)*c0->cd - c0->alp5*c*c1->qrx*(2.0 + a7);           
	du[ 4] = 3.0/c1->r5*(c0->alp4*y*c1->a5*c0->cd + 
	         c0->alp5*c*(c1->a5*c0->sd - y*qr5*a7) );          
	du[ 5] = 3.0/c1->r5*(-c0->alp4*y*c1->a5*c0->sd + 
	         c0->alp5*c*(c1->a5*c0->cd + d*qr5*a7) );        
	du[ 6] = du[4];                                                   
	du[ 7] = 3.0*x/c1->r5*(c0->alp4*b5*c0->cd - 
	         c0->alp5*5.0*c/c1->r2*(2.0*y*c0->sd+c1->q*b7) ); 
	du[ 8] = 3.0*x/c1->r5*(-c0->alp4*b5*c0->sd + 
	         c0->alp5*5.0*c/c1->r2*(d*b7*c0->sd - y*c7*c0->cd) );
	du[ 9] = 3.0/c1->r5*(-c0->alp4*d*c1->a5*c0->cd + 
	         c0->alp5*c*(c1->a5*c0->cd + d*qr5*a7) );      
	du[10] = 15.0*x/r7*(c0->alp4*y*d*c0->cd  + c0->alp5*c*(d*b7*c0->sd - y*c7*c0->cd) );    
	du[11] = 15.0*x/r7*(-c0->alp4*y*d*c0->sd  + c0->alp5*c*(2.0*d*c0->cd - c1->q*c7) );        

	for (i=0; i<12; i++) u[i] = u[i] + pot1/pi2*du[i];                                        
  }                                                             

/*
C===================================                                    
C=====  DIP-SLIP CONTRIBUTION  =====                                    
C===================================                                    
*/

  if (pot2 != 0.0) 
  {                                              
	du[ 0] = c0->alp4*3.0*x*c1->t/c1->r5 - c0->alp5*c*c1->p*c1->qrx;                   
	du[ 1] = -c0->alp4/c1->r3*(c0->c2d-3.0*y*c1->t/c1->r2) + 
	         c0->alp5*3.0*c/c1->r5*(c1->s - y*c1->p*qr5);       
	du[ 2] = -c0->alp4*c1->a3/c1->r3*c0->sdcd + c0->alp5*3.0*c/c1->r5*(c1->t + d*c1->p*qr5);       
	du[ 3] = c0->alp4*3.0*c1->t/c1->r5*c1->a5 - c0->alp5*5.0*c*c1->p*c1->qr/c1->r2*a7;      
	du[ 4] = 3.0*x/c1->r5*(c0->alp4*(c0->c2d - 5.0*y*c1->t/c1->r2) -
	                       c0->alp5*5.0*c/c1->r2*(c1->s - y*c1->p*qr7)); 
	du[ 5] = 3.0*x/c1->r5*(c0->alp4*(2.0 + c1->a5)*c0->sdcd -
	                       c0->alp5*5.0*c/c1->r2*(c1->t+d*c1->p*qr7)); 
	du[ 6] = du[4];               
	du[ 7] = 3.0/c1->r5*(c0->alp4*(2.0*y*c0->c2d + c1->t*b5) +                             
                         c0->alp5*c*(c0->s2d - 10.0*y*c1->s/c1->r2 - c1->p*qr5*b7)); 
	du[ 8] = 3.0/c1->r5*(c0->alp4*y*c1->a5*c0->sdcd - 
	                     c0->alp5*c*((3.0 + c1->a5)*c0->c2d + y*c1->p*dr5*qr7)); 
	du[ 9] = 3.0*x/c1->r5*(-c0->alp4*(c0->s2d - c1->t*dr5) - 
	                       c0->alp5*5.0*c/c1->r2*(c1->t+d*c1->p*qr7));   
	du[10] = 3.0/c1->r5*(-c0->alp4*(d*b5*c0->c2d + y*c5*c0->s2d) -                        
                        c0->alp5*c*((3.0 + c1->a5)*c0->c2d + y*c1->p*dr5*qr7));
	du[11] = 3.0/c1->r5*(-c0->alp4*d*c1->a5*c0->sdcd - 
	                     c0->alp5*c*(c0->s2d - 10.0*d*c1->t/c1->r2 + c1->p*qr5*c7));

	for (i=0; i<12; i++) u[i] =u[i] + pot2/pi2*du[i];                                        
  }                                                             

/*
C========================================                               
C=====  TENSILE-FAULT CONTRIBUTION  =====                               
C========================================                               
*/

  if(pot3 != 0.0)
  {                                               
	du[ 0] = 3.0*x/c1->r5*(-c0->alp4*c1->s + c0->alp5*(c*c1->q*qr5 - z));                     
	du[ 1] = c0->alp4/c1->r3*(c0->s2d - 3.0*y*c1->s/c1->r2) + 
	         c0->alp5*3.0/c1->r5*(c*(c1->t - y + y*c1->q*qr5) - y*z);
	du[ 2] =-c0->alp4/c1->r3*(1.0 - c1->a3*c0->sdsd) -
	                          c0->alp5*3.0/c1->r5*(c*(c1->s - d + d*c1->q*qr5) - d*z);
	du[ 3] =-c0->alp4*3.0*c1->s/c1->r5*c1->a5 +
	        c0->alp5*(c*c1->qr*qr5*a7 - 3.0*z/c1->r5*c1->a5);          
	du[ 4] = 3.0*x/c1->r5*(-c0->alp4*(c0->s2d - 5.0*y*c1->s/c1->r2) -                          
                           c0->alp5*5.0/c1->r2*(c*(c1->t - y + y*c1->q*qr7) - y*z)); 
	du[ 5] = 3.0*x/c1->r5*(c0->alp4*(1.0 - (2.0 + c1->a5)*c0->sdsd) +                        
                           c0->alp5*5.0/c1->r2*(c*(c1->s - d + d*c1->q*qr7) - d*z)); 
	du[ 6] = du[4];                                                   
	du[ 7] = 3.0/c1->r5*(-c0->alp4*(2.0*y*c0->s2d + c1->s*b5) -                            
                         c0->alp5*(c*(2.0*c0->sdsd + 10.0*y*(c1->t - y)/c1->r2 - 
                                      c1->q*qr5*b7) + z*b5)); 
	du[ 8] = 3.0/c1->r5*(c0->alp4*y*(1.0 - c1->a5*c0->sdsd) +                            
                         c0->alp5*(c*(3.0 + c1->a5)*c0->s2d - y*dr5*(c*d7 + z)));             
	du[ 9] = 3.0*x/c1->r5*(-c0->alp4*(c0->c2d + c1->s*dr5) +                              
                           c0->alp5*(5.0*c/c1->r2*(c1->s - d + d*c1->q*qr7) - 1.0 - z*dr5));            
	du[10] = 3.0/c1->r5*(c0->alp4*(d*b5*c0->s2d - y*c5*c0->c2d) +                        
                         c0->alp5*(c*((3.0 + c1->a5)*c0->s2d - y*dr5*d7) -
                                   y*(1.0 + z*dr5)));     
	du[11] = 3.0/c1->r5*(-c0->alp4*d*(1.0 - c1->a5*c0->sdsd) -                             
                         c0->alp5*(c*(c0->c2d + 10.0*d*(c1->s-d)/c1->r2 - c1->q*qr5*c7) +
                                   z*(1.0 + c5))); 

	for (i=0; i<12; i++) u[i] = u[i] + pot3/pi2*du[i];                                        
  }                                                             

/*
C=========================================                              
C=====  INFLATE SOURCE CONTRIBUTION  =====                              
C=========================================                              
*/

  if (pot4 != 0.0)
  {                                               
	du[ 0] = c0->alp4*3.0*x*d/c1->r5;                                          
	du[ 1] = c0->alp4*3.0*y*d/c1->r5;                                          
	du[ 2] = c0->alp4*c1->c3/c1->r3;                                              
	du[ 3] = c0->alp4*3.0*d/c1->r5*c1->a5;                                         
	du[ 4] = -c0->alp4*15.0*c1->xy*d/r7;                                        
	du[ 5] = -c0->alp4*3.0*x/c1->r5*c5;                                         
	du[ 6] = du[4];                                                   
	du[ 7] = c0->alp4*3.0*d/c1->r5*b5;                                         
	du[ 8] = -c0->alp4*3.0*y/c1->r5*c5;                                         
	du[ 9] = du[5];                                                   
	du[10] = du[8];                                                  
	du[11] = c0->alp4*3.0*d/c1->r5*(2.0 + c5);                                    

	for (i=0; i<12; i++) u[i] = u[i] + pot4/pi2*du[i];                                        
  }                                                             

  return;  
}                                                          

/* this version is so R can call it */
void dc3d0_passByRef(double *alpha_p, double *x_p, double *y_p, double *z_p,
          double *depth_p, double *dip_p, double *pot1_p, double *pot2_p,
          double *pot3_p, double *pot4_p, double *ux, double *uy, double *uz,
          double *uxx, double *uyx, double *uzx, double *uxy,
          double *uyy, double *uzy, double *uxz, double *uyz,
          double *uzz)
{
  double alpha, x, y, z, depth, dip, pot1, pot2, pot3, pot4;
  int err;
  
  alpha = *alpha_p;
  x = *x_p;
  y = *y_p;
  z = *z_p;
  depth = *depth_p;
  dip = *dip_p;
  pot1 = *pot1_p;
  pot2 = *pot2_p;
  pot3 = *pot3_p;
  pot4 = *pot4_p;
    
  err = dc3d0(alpha, x, y, z, depth, dip, pot1, pot2, pot3, pot4,
              ux, uy, uz, uxx, uyx, uzx, uxy, uyy, uzy, uxz, uyz, uzz);
  
    
  return;
}
  
/* dc3d0_passByRef_v is for use from R, so that it can do the loop
   over observation points here in C instead of in R.
   Currently it just loops over the the obs. points calling
   the single-point dc3d)() for each one.  Eventually, should
   save all the set up work instead of having dc3d0() do it
   each time 


   Args are same as for dc3d0_passByRef() except: x, y, and z are
   vectors of length npts (which is an additional arg),
   and the output args are now also vectors of length
   npts   
*/

void dc3d0_passByRef_v(
           double *alpha, double *x, double *y, double *z,
           double *depth, double *dip, double *pot1, double *pot2,
           double *pot3, double *pot4, double *ux, double *uy, double *uz,
           double *uxx, double *uyx, double *uzx, double *uxy,
           double *uyy, double *uzy, double *uxz, double *uyz,
           double *uzz, int *npts)  
{
  int ipt;
  
  for (ipt=0; ipt < *npts; ipt++)
    dc3d0(*alpha, x[ipt], y[ipt], z[ipt], *depth, *dip,  
          *pot1, *pot2, *pot3, *pot4, &(ux[ipt]), &(uy[ipt]), (&uz[ipt]),
          &(uxx[ipt]), &(uyx[ipt]), &(uzx[ipt]), &(uxy[ipt]), &(uyy[ipt]), &(uzy[ipt]),
          &(uxz[ipt]), &(uyz[ipt]), &(uzz[ipt]));
  
  return;
}
                                                     

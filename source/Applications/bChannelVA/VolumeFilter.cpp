#include "Postprocessing.h"
/*
* To change this license header, choose License Headers in Project Properties.
* To change this template file, choose Tools | Templates
* and open the template in the editor.
*/


/**
*
* @author kraft
*/


// Filtering of velocity field
// velocities are stored in ux[],uy[],uz[]
// coordinates are stored in coordx[],coordy[],coordz[]

//////////////////////////////////////////////////////////////////////////
double G(double x, double x_min, double x_max)
{
   double x_mean=(x_min+x_max)/2.0;
   if(x<0.0)
   {
      if(fabs(x)>x_min)
         return 0.0;
      else
         return  x*x_mean/x_min+x_mean;
   }
   else
   {
      if(x>x_max)
         return 0.0;
      else
         return -x*x_mean/x_max+x_mean;
   }
}
//////////////////////////////////////////////////////////////////////////
double m(double x, double y, double z, double x_min, double x_max, double y_min, double y_max, double z_min, double z_max, double norm)
{
   return G(x,x_min,x_max)*G(y,y_min,y_max)*G(z,z_min,z_max)/norm;
}
//////////////////////////////////////////////////////////////////////////
void VolumenFilter() 
{
   double norm, lx_mean,ly_mean,lz_mean;
   double steps=100.0;
   double  L=5,
      x_min, x_max,y_min,y_max,z_min,z_max,
      dx,dy,dz,
      dx3,
      X,
      Y,
      Z,
      integral=0.0;

   x_min=x_max=L/2.0;
   z_min=z_max=L/2.0;
   y_min=L/4.0;
   y_max=L/2.0;
   lx_mean=(x_min+x_max)/2.0;
   ly_mean=(y_min+y_max)/2.0;
   lz_mean=(z_min+z_max)/2.0;

   dx=(x_min+x_max)/steps;
   dy=(y_min+y_max)/steps;
   dz=(z_min+z_max)/steps;
   dx3=dx*dy*dz;
   norm=	 (lx_mean*lx_mean)
      *(ly_mean*ly_mean)
      *(lz_mean*lz_mean);
   X=30.0;
   Y=30.0;
   Z=30.0;


   for(double x=X-x_min;x<X+x_max;x+=dx){
      for(double y=Y-y_min;y<Y+y_max;y+=dy){
         for(double z=Z-z_min;z<Z+z_max;z+=dz){
            integral+=m(x+dx/2.0-X,y+dy/2.0-Y,z+dz/2.0-Z,x_min,x_max,y_min,y_max,z_min,z_max, norm)*dx3;
         }
      }
   }

   std::cout<<"Integralwert= " << integral<<std::endl;
}


#ifndef FEPOINT3D_H
#define FEPOINT3D_H

#include <sstream>
#include <numerics/geometry3d/GbPoint3D.h>


/*=========================================================================*/
/* GbPoint3D                                                               */
/*                                                                         */
/**
 * This Class provides basic 3D point objects.
 * <BR><BR><HR>
 * @author <A HREF="mailto:geller@bauinf.tu-cottbus.de">S. Geller</A>
 * @version 1.0 - 10.02.07
 * 
	*/
class FePoint3D : public GbPoint3D 
{
private:

	double Fx;
	double Fy;
	double Fz;
    //double oldFx;
    //double oldFy;
    double sumFx;
    double sumFy;
    double sumFz;
	double velocityX;
	double velocityY;
	double velocityZ;
    //double accelerationX;
    //double accelerationY;
   UbTupleDouble6 stresses; 

   double area;
public:
	FePoint3D():GbPoint3D()
   {
      this->init();
	}
	FePoint3D(double x, double y, double z):GbPoint3D(x,y,z)
	{
      this->init();
	}
   FePoint3D(GbPoint3D *point):GbPoint3D(point)
   {
      this->init();
   }

	FePoint3D(FePoint3D *point):GbPoint3D(point)
	{
      this->Fx = point->Fx;
      this->Fy = point->Fy;
      this->Fz = point->Fz;
      throw UbException(__FILE__,__LINE__,"spaeter fertig machen...");
   }

	virtual ~FePoint3D()
	{
	}

   void init()
   {
      this->Fx = 0.0;
      this->Fy = 0.0;
      this->Fz = 0.0;
      //this->oldFx = 0.0;
      //this->oldFy = 0.0;
      this->sumFx = 0.0;
      this->sumFy = 0.0;
      this->sumFz = 0.0;
      val<1>(stresses) = 0.0;
      val<2>(stresses) = 0.0;
      val<3>(stresses) = 0.0;
      val<4>(stresses) = 0.0;
      val<5>(stresses) = 0.0;
      val<6>(stresses) = 0.0;
      this->area = 0.0;
      this->velocityX = 0.0;
      this->velocityY = 0.0;
      this->velocityZ = 0.0;
   }
	/*======================================================================*/
   FePoint3D* clone()   
   { 
      return(new FePoint3D(this)); 
   }

	double getVelocityX()   { return(this->velocityX); }
	double getVelocityY()   { return(this->velocityY); }
	double getVelocityZ()   { return(this->velocityZ); }
	void setVelocityX(double x)   { this->velocityX = x; }
	void setVelocityY(double y)   { this->velocityY = y; }
	void setVelocityZ(double z)   { this->velocityZ = z; }

   double getFX()   { return(this->Fx); }
	double getFY()   { return(this->Fy); }
	double getFZ()   { return(this->Fz); }
	void setFX(double FX)   { this->Fx = FX; }
	void setFY(double FY)   { this->Fy = FY; }
	void setFZ(double FZ)   { this->Fz = FZ; }
	void addFX(double FX)   { this->Fx += FX; }
	void addFY(double FY)   { this->Fy += FY; }
	void addFZ(double FZ)   { this->Fz += FZ; }

   double getSumFX()   { return(this->sumFx); }
   double getSumFY()   { return(this->sumFy); }
   double getSumFZ()   { return(this->sumFz); }
   void setSumFX(double FX)   { this->sumFx = FX; }
   void setSumFY(double FY)   { this->sumFy = FY; }
   void setSumFZ(double FZ)   { this->sumFz = FZ; }
   void addSumFX(double FX)   { this->sumFx += FX; }
   void addSumFY(double FY)   { this->sumFy += FY; }
   void addSumFZ(double FZ)   { this->sumFz += FZ; }

   UbTupleDouble6& getStresses() { return this->stresses; }
//   void setS11(double S11) { this->S11 = S11; }
//   void setS12(double S12) { this->S12 = S12; }
//   void setS22(double S22) { this->S22 = S22; }
   double getArea() { return this->area; }
   void setArea(double area) { this->area = area; }
   void addArea(double area) { this->area += area; }

   /**
    * Returns a string representation of this 3D fe-point.
    * @return a string representation of this 3D fe-point
    */
   std::string toString()
   {
      std::stringstream ss; 
      ss<<"FePoint3D[";
		ss<<"x1="<<this->x1;
		ss<<", x2="<<this->x2;		
		ss<<", x3="<<this->x3;		
		ss<<", Fx="<<this->Fx;
		ss<<", Fy="<<this->Fy;
		ss<<", Fz="<<this->Fz<<"]";
		return((ss.str()).c_str());
   }
   /*======================================================================*/
};
/*=========================================================================*/
#endif




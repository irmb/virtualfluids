//  _    ___      __              __________      _     __
// | |  / (_)____/ /___  ______ _/ / ____/ /_  __(_)___/ /____
// | | / / / ___/ __/ / / / __ `/ / /_  / / / / / / __  / ___/
// | |/ / / /  / /_/ /_/ / /_/ / / __/ / / /_/ / / /_/ (__  )
// |___/_/_/   \__/\__,_/\__,_/_/_/   /_/\__,_/_/\__,_/____/
//
#ifndef COORDINATETRANSFORMATION3D_H
#define COORDINATETRANSFORMATION3D_H

#ifdef RCF_USE_BOOST_SERIALIZATION
   #include <boost/archive/text_oarchive.hpp>
   #include <boost/archive/text_iarchive.hpp>	
#endif //RCF_USE_BOOST_SERIALIZATION

#include <cmath>
#include <string>
#include <sstream>

#include <basics/utilities/UbException.h>
#include <basics/utilities/UbFileInput.h>
#include <basics/utilities/UbFileOutput.h>

#include <basics/memory/MbSharedPointerDefines.h>
class CoordinateTransformation3D;
typedef std::shared_ptr<CoordinateTransformation3D> CoordinateTransformation3DPtr;


//description:     x1/x2/x3 = alt, x1*/x2*/x3* = neu
//   x2      
//   ^             x*
//   |            /
//   |           2*
//   4          /
//   |         /
//   3        1*                     => neues coordsys ist um originX1=originX2=originX3=2 verschoben
//   |       /                          neues dx1=dx2=dx2=2 -> skalierung um 2 in x1-,x2- und x3-richtung
//   2      /                           ERST verdrehung um alpha um "x1" achse
//   |       \                          DANN verdrehung um beta  um "x2" achse
//   1         \                        DANN verdrehung um gamma um "x3" achse
//   |           x1*
//   |--1--2--3--4--5------------- > x1
//
// Bemerkung: kann sein, dass die Verdrehung um x1 und x3 vertauschst sind 
//            - muss mal einer prüfen ...



class CoordinateTransformation3D
{
public:
   CoordinateTransformation3D();
   CoordinateTransformation3D(const double& originX1, const double& originX2, const double& originX3, const double& dx1, const double& dx2, const double& dx3, const double& alpha, const double& beta, const double& gamma);
   CoordinateTransformation3D(const double& originX1, const double& originX2, const double& originX3, const double& dx1, const double& dx2, const double& dx3);
   CoordinateTransformation3D(CoordinateTransformation3D* transformation);
   
   void setTransformationValues(const double& originX1, const double& originX2, const double& originX3, const double& dx1, const double& dx2, const double& dx3, const double& alpha, const double& beta, const double& gamma);
   double getX1CoordinateOffset()  const { return this->Tx1;   }   //Translation
   double getX2CoordinateOffset()  const { return this->Tx2;   }
   double getX3CoordinateOffset()  const { return this->Tx3;   }
   double getX1CoordinateScaling() const { return this->Sx1;   }	 //Scaling
   double getX2CoordinateScaling() const { return this->Sx2;   }
   double getX3CoordinateScaling() const { return this->Sx3;   }
   double getRotationX1Angle()     const { return this->alpha; }
   double getRotationX2Angle()     const { return this->beta;  }
   double getRotationX3Angle()     const { return this->gamma; }	 //Rotation

   //Achtung die Winkel passen nicht überein -siehe setTransformationValues 
   void setRotationX1Angle(double alpha) { this->setTransformationValues(this->Tx1, this->Tx2, this->Tx3, this->Sx1, this->Sx2, this->Sx3, alpha, this->beta, this->gamma); }
   void setRotationX2Angle(double beta ) { this->setTransformationValues(this->Tx1, this->Tx2, this->Tx3, this->Sx1, this->Sx2, this->Sx3, this->alpha, beta, this->gamma); }
   void setRotationX3Angle(double gamma) { this->setTransformationValues(this->Tx1, this->Tx2, this->Tx3, this->Sx1, this->Sx2, this->Sx3, this->alpha, this->beta, gamma); }

   void setActive(const bool& active);
   bool isActive()          const { return this->active; }
   bool isTransformation()  const { return this->transformation; }

   double transformForwardToX1Coordinate(const double& x1, const double& x2, const double& x3) const;
   double transformForwardToX2Coordinate(const double& x1, const double& x2, const double& x3) const;
   double transformForwardToX3Coordinate(const double& x1, const double& x2, const double& x3) const;
   double transformForwardToX1CoordinateIgnoringRotation(const double& x1) const;
   double transformForwardToX2CoordinateIgnoringRotation(const double& x2) const;
   double transformForwardToX3CoordinateIgnoringRotation(const double& x3) const;
   double transformBackwardToX1Coordinate(const double& x1, const double& x2, const double& x3) const;
   double transformBackwardToX2Coordinate(const double& x1, const double& x2, const double& x3) const;
   double transformBackwardToX3Coordinate(const double& x1, const double& x2, const double& x3) const;
   double transformBackwardToX1CoordinateIgnoringRotation(const double& x1) const;
   double transformBackwardToX2CoordinateIgnoringRotation(const double& x2) const;
   double transformBackwardToX3CoordinateIgnoringRotation(const double& x3) const;
   std::string toString() const;

   //------------- implements CAB serialization ----- start
   void write(UbFileOutput* out) const
   {
      out->writeString("Coordtransfomartion3D");
      out->writeDouble(this->Tx1);
      out->writeDouble(this->Tx2);
      out->writeDouble(this->Tx3);
      out->writeDouble(this->Sx1);
      out->writeDouble(this->Sx2);
      out->writeDouble(this->Sx3);
      out->writeDouble(this->alpha);
      out->writeDouble(this->beta );
      out->writeDouble(this->gamma);
   }
   void read(UbFileInput* in)
   {
      in->readString();
      this->Tx1   = in->readDouble();
      this->Tx2   = in->readDouble();
      this->Tx3   = in->readDouble();
      this->Sx1   = in->readDouble();
      this->Sx2   = in->readDouble();
      this->Sx3   = in->readDouble();
      this->alpha = in->readDouble();
      this->beta  = in->readDouble();
      this->gamma = in->readDouble();

      this->setTransformationValues(Tx1,Tx2,Tx3,Sx1,Sx2,Sx3,alpha,beta,gamma);
   }
   //------------- implements CAB serialization ----- end

#ifdef CAB_RCF
   template<class Archive>
   void SF_SERIALIZE(Archive & ar)
   {
      ar & Tx1;   ar & Tx2; ar & Tx3; 
      ar & Sx1;   ar & Sx2; ar & Sx3; 
      ar & alpha; ar & beta; ar & gamma;
      
      ar & toX1factorX1; ar & toX1factorX2; ar & toX1factorX3; ar & toX1delta;
      ar & toX2factorX1; ar & toX2factorX2; ar & toX2factorX3; ar & toX2delta;
      ar & toX3factorX1; ar & toX3factorX2; ar & toX3factorX3; ar & toX3delta;

      ar & fromX1factorX1; ar & fromX1factorX2; ar & fromX1factorX3; ar & fromX1delta;
      ar & fromX2factorX1; ar & fromX2factorX2; ar & fromX2factorX3; ar & fromX2delta;
      ar & fromX3factorX1; ar & fromX3factorX2; ar & fromX3factorX3; ar & fromX3delta;

      ar & active;
      ar & transformation;
   }
#endif //CAB_RCF

private:
   double Tx1, Tx2, Tx3, Sx1, Sx2, Sx3, alpha, beta, gamma;

   double toX1factorX1, toX1factorX2, toX1factorX3, toX1delta;
   double toX2factorX1, toX2factorX2, toX2factorX3, toX2delta;
   double toX3factorX1, toX3factorX2, toX3factorX3, toX3delta;

   double fromX1factorX1, fromX1factorX2, fromX1factorX3, fromX1delta;
   double fromX2factorX1, fromX2factorX2, fromX2factorX3, fromX2delta;
   double fromX3factorX1, fromX3factorX2, fromX3factorX3, fromX3delta;

   bool   active;
   bool   transformation;

   friend class MPIIORestartCoProcessor;
   friend class MPIIOMigrationCoProcessor;
};

#endif //COORDINATETRANSFORMATION3D_H

//  _    ___      __              __________      _     __
// | |  / (_)____/ /___  ______ _/ / ____/ /_  __(_)___/ /____
// | | / / / ___/ __/ / / / __ `/ / /_  / / / / / / __  / ___/
// | |/ / / /  / /_/ /_/ / /_/ / / __/ / / /_/ / / /_/ (__  )
// |___/_/_/   \__/\__,_/\__,_/_/_/   /_/\__,_/_/\__,_/____/
//
#ifndef GBTRIANGLE3D_H
#define GBTRIANGLE3D_H

#include <sstream>

#include <numerics/geometry3d/GbObject3D.h>
#include <numerics/geometry3d/GbVector3D.h>
#include <numerics/geometry3d/GbPoint3D.h>

#ifdef CAB_RCF
   #include <3rdParty/rcf/RcfSerializationIncludes.h>
#endif //CAB_RCF

#include <basics/memory/MbSharedPointerDefines.h>
class GbTriangle3D;
typedef VFSharedPtr<GbTriangle3D> GbTriangle3DPtr;


class GbCuboid3D;
class GbPolygon3D;
class GbObject3DCreator;

/*=========================================================================*/
/* GbTriangle3D                                                            */
/*                                                                         */
/*                                                               
* This Class provides basic 3D triangle objects.
*/
//class GbLine2D;

class GbTriangle3D : public GbObject3D , public UbObserver
{
public:
   /*======================================================================*/
   /*  Konstruktoren                                                       */
   /*                                                                      */
   GbTriangle3D();
   GbTriangle3D(GbPoint3D* point1, GbPoint3D* point2, GbPoint3D* point3);
   GbTriangle3D(GbTriangle3D* triangle);
   ~GbTriangle3D();
   /*======================================================================*/
   /*  Methoden                                                            */
   /*                                                                      */
   GbTriangle3D* clone();
   void finalize()
   {
      this->deletePoints();
   }

   GbPoint3D* getPoint1()   { return this->points[0]; }
   GbPoint3D* getPoint2()   { return this->points[1]; }
   GbPoint3D* getPoint3()   { return this->points[2]; }

   GbVector3D getNormal();
   void       calculateNormal();

   void deletePoints();

   int contains(GbPoint3D* point);
   int containsEqual(GbPoint3D* point);
   GbPoint3D* getPoint(const int& index);
   std::vector<GbPoint3D> getPoints();
   double getArea();
   double getX1Centroid();
   double getX1Minimum();
   double getX1Maximum();           
   double getX2Centroid();
   double getX2Minimum();
   double getX2Maximum();
   double getX3Centroid();
   double getX3Minimum();
   double getX3Maximum();

   void setInconsistent() { this->consistent = false;}

   void setPoint(GbPoint3D *point, int index);

   //bool equals(GbObject3D *object)
   std::vector<GbTriangle3D*> getSurfaceTriangleSet();
   bool isPointInGbObject3D(const double& x1, const double& x2, const double& x3)   
   {
      //der einfachheit halber ... 
      return false;
      //throw UbException(__FILE__, __LINE__, "GbTriangle3D::isPointInObject3D- not implemented");
   }
   bool isPointInGbObject3D(const double& x1, const double& x2, const double& x3, bool& pointIsOnBoundary)   
   {
      //der einfachheit halber ... 
      pointIsOnBoundary = false;
      return false;
      //throw UbException(__FILE__, __LINE__, "GbTriangle3D::isPointInObject3D- not implemented");
   }
   bool isCellInsideGbObject3D(const double& x11,const double& x21,const double& x31,const double& x12,const double& x22,const double& x23) { return false; }


   // get distance from a point to the triangle
   //todo CHANGE...
   double getDistanceFromPoint(GbVector3D punct);

   std::string toString();
   ObObjectCreator* getCreator();
   void write(UbFileOutput* out);
   void read(UbFileInput* in);

   /*======================================================================*/
   /*  Calculation                                                         */
   /*                                                                      */
//   std::vector<GbPoint3D> calculateIntersectionPoints3D(GbLine3D *line);
   bool hasRaytracing() { return true; }
   /*|r| must be 1! einheitsvector!!*/
   double getIntersectionRaytraceFactor(const double& x1, const double& x2, const double& x3, const double& rx1, const double& rx2, const double& rx3);
//   bool isPointOnEdge(GbVector3D& q);
   
   GbPoint3D* calculateIntersectionPoints3D(GbLine3D* line);
   GbPoint3D* calculateIntersectionPoints3D(GbPoint3D* linePoint1, GbPoint3D* linePoint2);
   double calculateDistanceToPoint3D(GbPoint3D *point);
   double calculateDistanceToPoint3D(const double& x1, const double& x2, const double& x3);      
   double calculateNormalizedDistanceToPoint3D(const double& x1, const double& y1, const double& z1, const double& x2, const double& y2, const double& z2);

   bool enclosesPoint2D(double x1, double x2);
   GbPolygon3D* createClippedPolygon3D(GbCuboid3D* cube);   
   GbLine3D* createClippedLine3D (GbPoint3D& point1, GbPoint3D& point2);
   //public GbPolygon2D createClippedPolygon2D(GbPoint2D p1, GbPoint2D p2);
   GbPolygon3D* createClippedPolygon3D(const double& p1x1, const double& p1x2, const double& p1x3, const double& p2x1, const double& p2x2, const double& p2x3);
   //bool enclosesRectangle2D(GbRectangle2D *rectangle);
   //public boolean enclosesRectangle2D(GbPoint2D p1, GbPoint2D p2);
   //public boolean enclosesRectangle2D(double p1x1, double p1x2, double p2x1, double p2x2);
   //public boolean crossesRectangle2D(GbRectangle2D rectangle);
   //public boolean crossesRectangle2D(GbPoint2D p1, GbPoint2D p2);
   //public boolean crossesRectangle2D(double p1x1, double p1x2, double p2x1, double p2x2);
   //public boolean enclosesOrCrossesRectangle2D(GbRectangle2D rectangle);
   //public boolean enclosesOrCrossesRectangle2D(GbPoint2D p1, GbPoint2D p2);
   //public boolean enclosesOrCrossesRectangle2D(double p1x1, double p1x2, double p2x1, double p2x2);
   /*======================================================================*/
   /*======================================================================*/
   /*  Private Methoden                                                    */
   /*                                                                      */
   virtual void calculateValues();

   /*======================================================================*/
   //class PointObserver : public UbObserver
   //{
   //    GbTriangle3D *triangle;

   //    PointObserver(GbTriangle3D *triangle)
   //    {
   //      this->triangle = triangle;
   //    }

   //public:
   //   void objectChanged(GbObject3D *object)
   //    {
   //      if(object == this->triangle->points[0] || object == this->triangle->points[1]  || object == this->triangle->points[2])
   //      {
   //         this->triangle->consistent = false;
   //         this->triangle->notifyObservers();
   //      }
   //    }
   //};
   /*======================================================================*/

   //virtuelle Methoden von UbObserver
   //!! quick and dirty von sirann !!
   void objectChanged(UbObservable* changedObject)
   {
      GbPoint3D* point = dynamic_cast<GbPoint3D*>(changedObject);
      if(!point || (  this->points[0]!=point && this->points[1]!=point && this->points[2]!=point) ) 
         return;
      
      this->consistent = false;
   }
   void objectWillBeDeleted(UbObservable* objectForDeletion)
   {
      if(this->points[0])
      {
         UbObservable* observedObj = dynamic_cast<UbObservable*>(this->points[0]);
         if(objectForDeletion == observedObj) { this->points[0] = NULL; }
      }
      if(this->points[1])
      {
         UbObservable* observedObj = dynamic_cast<UbObservable*>(this->points[1]);
         if(objectForDeletion == observedObj) { this->points[1] = NULL; }
      }
      if(this->points[2])
      {
         UbObservable* observedObj = dynamic_cast<UbObservable*>(this->points[2]);
         if(objectForDeletion == observedObj) { this->points[2] = NULL; }
      }
      //ACHTUNG: eigentlich muessten in allen methoden von GbLine if abfragen fuer NULL pointer hin... toDo
   }
   using GbObject3D::isPointInGbObject3D; //Grund: dadurch muss man hier  isPointInGbObject3D(GbPoint3D*) nicht ausprogrammieren, welche sonst hier "ueberdeckt" waere

#ifdef CAB_RCF
   template<class Archive>
   void SF_SERIALIZE(Archive & ar)
   {
      SF_SERIALIZE_PARENT<GbObject3D>(ar, *this);
      ar & points;
      ar & normal;
      ar & x1s;
      ar & x2s;
      ar & x3s;
      ar & x1min;
      ar & x1max;
      ar & x2min;
      ar & x2max;
      ar & x3min;
      ar & x3max;
      ar & area;
      ar & consistent;
      if( ArchiveTools::isReading(ar) ) this->calculateNormal();
   }
#endif //CAB_RCF

protected:
   bool   consistent;
   double x1s;
   double x2s;
   double x3s;
   double x1min;
   double x1max;
   double x2min;
   double x2max;
   double x3min;
   double x3max;
   double area;
   
   GbVector3D normal;
   std::vector<GbPoint3D*> points;
   
private:
   void init();
};
/*=========================================================================*/

#if defined(RCF_USE_SF_SERIALIZATION) && !defined(SWIG)
   UB_AUTO_RUN_NAMED(   SF::registerType<GbTriangle3D  >("GbTriangle3D  ")        , SF_GbTriangle3D     );
   UB_AUTO_RUN_NAMED( ( SF::registerBaseAndDerived< GbObject3D, GbTriangle3D >() ), SF_GbTriangle3D_BD1 );
   UB_AUTO_RUN_NAMED( ( SF::registerBaseAndDerived< UbObserver, GbTriangle3D>()  ), SF_GbTriangle3D_BD2 );
#endif //RCF_USE_SF_SERIALIZATION

#endif

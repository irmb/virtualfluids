//  _    ___      __              __________      _     __
// | |  / (_)____/ /___  ______ _/ / ____/ /_  __(_)___/ /____
// | | / / / ___/ __/ / / / __ `/ / /_  / / / / / / __  / ___/
// | |/ / / /  / /_/ /_/ / /_/ / / __/ / / /_/ / / /_/ (__  )
// |___/_/_/   \__/\__,_/\__,_/_/_/   /_/\__,_/_/\__,_/____/
//
#ifndef GBPOLYGON3D_H
#define GBPOLYGON3D_H

#include <sstream>
#include <iostream>


#include <numerics/geometry3d/GbObject3D.h>
#include <numerics/geometry3d/GbLine3D.h>
#include <numerics/geometry3d/GbTriangle3D.h>
#include <numerics/geometry3d/GbSystem3D.h>

#include <basics/memory/MbSharedPointerDefines.h>
class GbPolygon3D;
typedef VFSharedPtr<GbPolygon3D> GbPolygon3DPtr;


class GbObject3DCreator;

/*=========================================================================*/
/* GbPolygon2D                                                             */
/*                                                                         */
/*
* This Class provides basic 3D polygon objects.
*/
class GbPolygon3D : public GbObject3D
{
public:
   using GbObject3D::isPointInGbObject3D; //Grund: dadurch muss man hier  isPointInGbObject3D(GbPoint3D*) nicht ausprogrammieren, welche sonst hier "ueberdeckt" waere
private:
   /*======================================================================*/
   double            x1s  ;
   double            x2s  ;
   double            x3s  ;
   double            x1min;
   double            x1max;
   double            x2min;
   double            x2max;
   double            x3min;
   double            x3max;

   std::vector<GbPoint3D> points;
   bool                   consistent;

   GbSystem3D::PointSet3 *ps;
   //private PointObserver     po         = null;

   void init();

   /*======================================================================*/


   /*======================================================================*/
   /*  Konstruktoren                                                       */
   /*                                                                      */
   /*
   * Creates an empty 2D polygon.
   */
public:
   static int counter;
   GbPolygon3D();
   /*
   * Creates an empty 2D polygon with the specified capacity.
   * @param capacity the initial capacity
   */
   GbPolygon3D(int capacity);
   /*
   * Creates a 2D polygon with the specified points.
   * @param points the initial points of the polygon
   */
   GbPolygon3D(std::vector<GbPoint3D> &points);
   /*
   * Creates a 2D polygon as clone of the specified 2D polygon.
   * @param polygon the 2D polygon to be cloned
   */
   GbPolygon3D(GbPolygon3D *polygon);

   ~GbPolygon3D();

   /*======================================================================*/


   /*======================================================================*/
   /*  Methoden                                                            */
   /*                                                                      */
   /*
   * Creates a 2D polygon as clone of this 2D polygon.
   */
   GbPolygon3D* clone() {   return(new GbPolygon3D(this)); }
   void finalize()
   {
      throw UbException(UB_EXARGS,"toDo");
   }

   /*
   * Returns the number of points.
   * @return the number of points
   */
   int size();
   /*
   * Returns the number of times this 2D polygon contains the specified point.
   * @param point the point
   * @return the number of times this 2D polygon contains the specified point
   */
   int contains(GbPoint3D *point);
   /*
   * Returns the number of times this 2D polygon contains a point equal to the specified point.
   * @param point the point
   * @return the number of times this 2D polygon contains a point equal to the specified point
   */
   int containsEqual(GbPoint3D* point);
   /*
   * Returns true, if this 2D polygon contains the specified line.
   * @param point1 the first point
   * @param point2 the second point
   * @return true, if this 2D polygon contains the specified line
   */
   bool containsLine(GbPoint3D* point1, GbPoint3D* point2);
   /*
   * Returns true, if this 2D polygon contains the specified line.
   * @param line the line
   * @return true, if this 2D polygon contains the specified line
   */
   bool containsLine(GbLine3D* line);
   /*
   * Returns the first point.
   * @return the first point
   */
   GbPoint3D* getFirstPoint();
   /*
   * Returns the last point.
   * @return the last point
   */
   GbPoint3D* getLastPoint();
   /*
   * Returns the specified point.
   * @param index the index
   * @return the specified point
   * @exception ArrayIndexOutOfBoundsException if the specified index is not valid
   */
   GbPoint3D* getPoint(const int& index);
   /*
   * Returns the points.
   * @return the points
   */
   std::vector<GbPoint3D> getPoints();
   /*
   * Returns the points within the specified rectangle.
   * @param rectangle the 2D rectangle
   * @return the points within the specified rectangle
   */
   //public GbPoint2D[] getPoints(GbRectangle2D rectangle)
   //{
   //   return(this.getPoints(rectangle.p1.x1, rectangle.p1.x2, rectangle.p2.x1, rectangle.p2.x2));
   //}
   /*
   * Returns the points within the specified rectangle.
   * @param p1 the 1st point of the rectangle
   * @param p2 the 2nd point of the rectangle
   * @return the points within the specified rectangle
   */
   std::vector<GbPoint3D> getPoints(GbPoint3D* p1, GbPoint3D* p2);
   /*
   * Returns the points within the specified rectangle.
   * @param p1x1 the 1st x1 coordinate of the rectangle
   * @param p1x2 the 1st x2 coordinate of the rectangle
   * @param p2x1 the 2nd x1 coordinate of the rectangle
   * @param p2x2 the 2nd x2 coordinate of the rectangle
   * @return the points within the specified rectangle
   */
   std::vector<GbPoint3D> getPoints(const double& p1x1, const double& p1x2, const double& p1x3, const double& p2x1, const double& p2x2, const double& p2x3);
   /*
   * Returns the area of this polygon.
   * The area is positive for positive ordered points, otherwise negative.
   * @return the area of this polygon
   */
   //double getArea()
   //{
   //   if(!this.consistent) this.calculateValues();
   //   return(this.area);
   //}
   double getX1Centroid();
   double getX1Minimum();
   double getX1Maximum();
   double getX2Centroid();
   double getX2Minimum();
   double getX2Maximum();
   double getX3Centroid();
   double getX3Minimum();
   double getX3Maximum();

   /*
   * Adds a point to the end of this polygon. Notifies the observers of this 2D polygon.
   * @param point the point
   */
   void addPoint(GbPoint3D* point);
   /*
   * Adds a number of points to the end of this polygon. Notifies the observers of this 2D polygon.
   * @param points the points
   */
   void addPoints(std::vector<GbPoint3D>& points);
   /*
   * Inserts a point at the specified position within this polygon. Notifies the observers of this 2D polygon.
   * @param point the point
   * @param index the index
   * @exception ArrayIndexOutOfBoundsException if the specified index is not valid
   */
   //public void insertPoint(GbPoint2D point, int index) throws ArrayIndexOutOfBoundsException
   //{
   //   if((this instanceof GbPolygon3D) && !(point instanceof GbPoint3D)) throw new IllegalArgumentException("GbPolygon2D.insertPoint(): points of 3D polygons have to be 3D points!");
   //   if(index < 0 || index > this.ps.size()) throw new ArrayIndexOutOfBoundsException("GbPolygon2D.insert(): invalid index specified: "+index);

   //   this.ps.insert(point, index);
   //   point.addObserver(this.po);
   //   this.consistent = false;
   //   super.notifyObservers();
   //}
   /*
   * Removes all points from this polygon identical to the specified one. Notifies the observers of this 2D polygon.
   * @param point the point
   */
   //public void deletePoint(GbPoint2D point)
   //{
   //   this.ps.delete(point);
   //   point.removeObserver(this.po);
   //   this.consistent = false;
   //   super.notifyObservers();
   //}
   /*
   * Removes all points from this polygon equal to the specified one. Notifies the observers of this 2D polygon.
   * @param point the point
   */
   //public void deleteEqualPoint(GbPoint2D point)
   //{
   //   this.ps.deleteEqual(point);
   //   point.removeObserver(this.po);
   //   this.consistent = false;
   //   super.notifyObservers();
   //}
   /*
   * Removes all points from this polygon. Notifies the observers of this 2D polygon.
   */
   void clear();

   /*
   * Returns true if this 2D polygon equals the specified object.
   * Two polygon are equal, if their points are equal.
   * <BR>Note that the order of points is recognized!
   * @return true if this 2D polygon equals the specified object
   * @see GbPoint2D#equals(java.lang.Object)
   * @see GbPoint3D#equals(java.lang.Object)
   */
   // bool equals(Object object)
   // {
   //    try
   //    {
   //	GbPolygon2D polygon = (GbPolygon2D) object;
   //int         n       = this.size();

   //if(n != polygon.size()) return(false);
   //for(int i=0; i<n; i++) if(!this.getPoint(i).equals(polygon.getPoint(i))) return(false);
   //return(true);
   //    }
   //    catch(Exception e){ return(false); }
   // }
   std::vector<GbTriangle3D*> getSurfaceTriangleSet()
   {
      std::cout<<"GbPolygon3D::getSurfaceTriangleSet() - not implemented\n";
      std::vector<GbTriangle3D*> tmp;
      return tmp;
   }
   bool isPointInGbObject3D(const double& x1, const double& x2, const double& x3)
   {
      throw UbException(__FILE__, __LINE__, "GbPolygon3D::isPointInObject3D- not implemented");
   }
   bool isPointInGbObject3D(const double& x1, const double& x2, const double& x3, bool& pointIsOnBoundary)
   {
      throw UbException(__FILE__, __LINE__, "GbPolygon3D::isPointInObject3D- not implemented");
   }
   bool isCellInsideGbObject3D(double x11,double x21,double x31,double x12,double x22,double x32) { return false; }

   GbLine3D* createClippedLine3D (GbPoint3D& point1, GbPoint3D &point2)
   {
      throw UbException(__FILE__, __LINE__, "GbPolygon3D::createClippedLine3D - not implemented");
   }                        
/*
   * Returns a string representation of this 2D polygon.
   * @return a string representation of this 2D polygon
   */
   std::string toString();
   ObObjectCreator* getCreator();
   void write(UbFileOutput* out);
   void read(UbFileInput* in);
   /*======================================================================*/


   /*======================================================================*/
   /*  Calculation                                                         */
   /*                                                                      */
   /*
   * Returns the intersection points of this 2D polygon and the specified 2D line.
   * @param line the 2D line to intersect
   * @return the intersection points of this 2D polygon and the specified 2D line
   */
   // public GbPoint2D[] calculateIntersectionPoints2D(GbLine2D line)
   // {
   //    GbSystem.PointSet pointSet = new GbSystem.PointSet(0);
   //    GbPoint2D         points[] = this.getPoints();
   //    GbPoint2D         pCrossed = null;
   //    int               n        = points.length;
   //    if(n < 2)         return(pointSet.getPoints());

   //    for(int i=1; i<n; i++)
   //    {
   //pCrossed = GbSystem.calculateIntersectionPoint2D(points[i-1], points[i], line.p1, line.p2);
   //if(pCrossed != null) pointSet.add(pCrossed);
   //    }
   //    pCrossed = GbSystem.calculateIntersectionPoint2D(points[n-1], points[0], line.p1, line.p2);
   //    if(pCrossed != null) pointSet.add(pCrossed);

   //    return(pointSet.getPoints());
   // }

   /*
   * Returns true if the specified 2D point lies within (or on the border of) this 2D polygon.
   * @param point the 2D point to check
   * @return true if the specified 2D point lies within (or on the border of) this 2D polygon
   */
   // public boolean enclosesPoint2D(GbPoint2D point)
   // {
   //    if(GbSystem.less(point.x1, this.x1min))    return(false);
   //    if(GbSystem.less(point.x2, this.x2min))    return(false);
   //    if(GbSystem.greater(point.x1, this.x1max)) return(false);
   //    if(GbSystem.greater(point.x2, this.x2max)) return(false);
   //    if(this.containsEqual(point) > 0)          return(true);

   //    QbList    ltest    = new QbList(GbPoint2D.class, QbList.NOEQUALOBJECTS);
   //    GbPoint2D points[] = this.getPoints();
   //    GbPoint2D ptest;
   //    int       n        = points.length;
   //    if(n < 2) return(false);

   //    if(GbSystem.equal(point.x2, this.x2min)) ptest = new GbPoint2D(point.x1, this.x2min-1.0);
   //    else                                     ptest = new GbPoint2D(point.x1, this.x2max+1.0);

   //    for(int i=1; i<n; i++)
   //    {
   //try { ltest.append(GbSystem.calculateIntersectionPoint2D(points[i-1], points[i], point, ptest)); }
   //catch(Exception e){}
   //    }
   //    try { ltest.append(GbSystem.calculateIntersectionPoint2D(points[n-1], points[0], point, ptest)); }
   //    catch(Exception e){}
   //    return((ltest.size()%2)==1);
   // }

   /*
   * Returns a new 2D polygon clipped by the specified 2D rectangle (result may be null!).
   * @param rectangle the 2D rectangle
   * @return a new 2D polygon clipped by the specified 2D rectangle
   */
   // GbPolygon3D *createClippedPolygon3D(GbCuboid3D *cube)
   // {
   //return(GbSystem::clipPolygon3D(this->getPoints(), cube->p1->x1, cube->p1->x2, cube->p1->x3, , cube->p2->x1, cube->p2->x2, cube->p2->x3));
   // };
   /*
   * Returns a new 2D polygon clipped by the specified 2D rectangle (result may be null!).
   * @param p1 the 1st point of the rectangle
   * @param p2 the 2nd point of the rectangle
   * @return a new 2D polygon clipped by the specified 2D rectangle
   */
   // GbPolygon3D *createClippedPolygon2D(GbPoint3D *p1, GbPoint3D *p2)
   // {
   //return(GbSystem::clipPolygon3D(this->getPoints(), p1->x1, p1->x2, p1->x3, p2->x1, p2->x2, p2->x3));
   // };
   /*
   * Returns a new 2D polygon clipped by the specified 2D rectangle (result may be null!).
   * @param p1x1 the 1st x1 coordinate of the rectangle
   * @param p1x2 the 1st x2 coordinate of the rectangle
   * @param p2x1 the 2nd x1 coordinate of the rectangle
   * @param p2x2 the 2nd x2 coordinate of the rectangle
   * @return a new 2D polygon clipped by the specified 2D rectangle
   */
   // GbPolygon3D *createClippedPolygon3D(double p1x1, double p1x2, double p1x3, double p2x1, double p2x2, double p2x3)
   // {
   //return(GbSystem::clipPolygon3D(this.getPoints(), p1x1, p1x2, p1x3, p2x1, p2x2. p2x3));
   // };

   /*
   * Returns true if the specified 2D rectangle lies completely within this 2D polygon.
   * @param rectangle the 2D rectangle to check
   * @return true if the specified 2D rectangle lies completely within this 2D polygon
   */
   //public boolean enclosesRectangle2D(GbRectangle2D rectangle)
   //{
   //   GbPolygon2D p = GbSystem.clipPolygon2D(this.getPoints(), rectangle.p1.x1, rectangle.p1.x2, rectangle.p2.x1, rectangle.p2.x2);
   //   return(p!=null && GbSystem.equal(Math.abs(p.getArea()), rectangle.getArea()));
   //}
   /*
   * Returns true if the specified 2D rectangle lies completely within this 2D polygon.
   * @param p1 the 1st point of the rectangle to check
   * @param p2 the 2nd point of the rectangle to check
   * @return true if the specified 2D rectangle lies completely within this 2D polygon
   */
   //public boolean enclosesRectangle2D(GbPoint2D p1, GbPoint2D p2)
   //{
   //   GbPolygon2D p = GbSystem.clipPolygon2D(this.getPoints(), p1.x1, p1.x2, p2.x1, p2.x2);
   //   return(p!=null && GbSystem.equal(Math.abs(p.getArea()), Math.abs((p1.x1-p2.x1)*(p1.x2-p2.x2))));
   //}
   /*
   * Returns true if the specified 2D rectangle lies completely within this 2D polygon.
   * @param p1x1 the 1st x1 coordinate of the rectangle to check
   * @param p1x2 the 1st x2 coordinate of the rectangle to check
   * @param p2x1 the 2nd x1 coordinate of the rectangle to check
   * @param p2x2 the 2nd x2 coordinate of the rectangle to check
   * @return true if the specified 2D rectangle lies completely within this 2D polygon
   */
   //public boolean enclosesRectangle2D(double p1x1, double p1x2, double p2x1, double p2x2)
   //{
   //   GbPolygon2D p = GbSystem.clipPolygon2D(this.getPoints(), p1x1, p1x2, p2x1, p2x2);
   //   return(p!=null && GbSystem.equal(Math.abs(p.getArea()), Math.abs((p1x1-p2x1)*(p1x2-p2x2))));
   //}

   /*
   * Returns true if the specified 2D rectangle is crossed by this 2D polygon.
   * @param rectangle the 2D rectangle to check
   * @return true if the specified 2D rectangle is crossed by this 2D polygon
   */
   //public boolean crossesRectangle2D(GbRectangle2D rectangle)
   //{
   //   GbPolygon2D p = GbSystem.clipPolygon2D(this.getPoints(), rectangle.p1.x1, rectangle.p1.x2, rectangle.p2.x1, rectangle.p2.x2);
   //   return(p!=null && GbSystem.inOpenInterval(Math.abs(p.getArea()), 0.0, rectangle.getArea()));
   //}
   /*
   * Returns true if the specified 2D rectangle is crossed by this 2D polygon.
   * @param p1 the 1st point of the rectangle to check
   * @param p2 the 2nd point of the rectangle to check
   * @return true if the specified 2D rectangle is crossed by this 2D polygon
   */
   //public boolean crossesRectangle2D(GbPoint2D p1, GbPoint2D p2)
   //{
   //   GbPolygon2D p = GbSystem.clipPolygon2D(this.getPoints(), p1.x1, p1.x2, p2.x1, p2.x2);
   //   return(p!=null && GbSystem.inOpenInterval(Math.abs(p.getArea()), 0.0, Math.abs((p1.x1-p2.x1)*(p1.x2-p2.x2))));
   //}
   /*
   * Returns true if the specified 2D rectangle is crossed by this 2D polygon.
   * @param p1x1 the 1st x1 coordinate of the rectangle to check
   * @param p1x2 the 1st x2 coordinate of the rectangle to check
   * @param p2x1 the 2nd x1 coordinate of the rectangle to check
   * @param p2x2 the 2nd x2 coordinate of the rectangle to check
   * @return true if the specified 2D rectangle is crossed by this 2D polygon
   */
   //public boolean crossesRectangle2D(double p1x1, double p1x2, double p2x1, double p2x2)
   //{
   //   GbPolygon2D p = GbSystem.clipPolygon2D(this.getPoints(), p1x1, p1x2, p2x1, p2x2);
   //   return(p!=null && GbSystem.inOpenInterval(Math.abs(p.getArea()), 0.0, Math.abs((p1x1-p2x1)*(p1x2-p2x2))));
   //}

   /*
   * Returns true if the specified 2D rectangle lies (at least partly) within this 2D polygon.
   * @param rectangle the 2D rectangle to check
   * @return true if the specified 2D rectangle lies (at least partly) within this 2D polygon
   */
   //public boolean enclosesOrCrossesRectangle2D(GbRectangle2D rectangle)
   //{
   //   GbPolygon2D p = GbSystem.clipPolygon2D(this.getPoints(), rectangle.p1.x1, rectangle.p1.x2, rectangle.p2.x1, rectangle.p2.x2);
   //   return(p!=null && GbSystem.greater(Math.abs(p.getArea()), 0.0));
   //}
   /*
   * Returns true if the specified 2D rectangle lies (at least partly) within this 2D polygon.
   * @param p1 the 1st point of the rectangle to check
   * @param p2 the 2nd point of the rectangle to check
   * @return true if the specified 2D rectangle lies (at least partly) within this 2D polygon
   */
   //public boolean enclosesOrCrossesRectangle2D(GbPoint2D p1, GbPoint2D p2)
   //{
   //   GbPolygon2D p = GbSystem.clipPolygon2D(this.getPoints(), p1.x1, p1.x2, p2.x1, p2.x2);
   //   return(p!=null && GbSystem.greater(Math.abs(p.getArea()), 0.0));
   //}
   /*
   * Returns true if the specified 2D rectangle lies (at least partly) within this 2D polygon.
   * @param p1x1 the 1st x1 coordinate of the rectangle to check
   * @param p1x2 the 1st x2 coordinate of the rectangle to check
   * @param p2x1 the 2nd x1 coordinate of the rectangle to check
   * @param p2x2 the 2nd x2 coordinate of the rectangle to check
   * @return true if the specified 2D rectangle lies (at least partly) within this 2D polygon
   */
   //public boolean enclosesOrCrossesRectangle2D(double p1x1, double p1x2, double p2x1, double p2x2)
   //{
   //   GbPolygon2D p = GbSystem.clipPolygon2D(this.getPoints(), p1x1, p1x2, p2x1, p2x2);
   //   return(p!=null && GbSystem.greater(Math.abs(p.getArea()), 0.0));
   //}
   /*======================================================================*/


   /*======================================================================*/
   /*  Private Methoden                                                    */
   /*                                                                      */
   void calculateValues();
   /*======================================================================*/


   /*======================================================================*/
   // private class PointObserver implements TiObserver
   // {
   //    GbPolygon2D polygon = null;

   //    PointObserver(GbPolygon2D polygon)
   //    {
   //this.polygon = polygon;
   //    }

   //    public void objectChanged(Object object)
   //    {
   //if((object instanceof GbPoint2D) && this.polygon.contains((GbPoint2D)object)>0)
   //{
   //   this.polygon.consistent = false;
   //   this.polygon.notifyObservers();
   //}
   //    }
   // }
   /*======================================================================*/
};
/*=========================================================================*/
#endif







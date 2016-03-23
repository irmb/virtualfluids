#include <numerics/geometry3d/GbPolygon3D.h>
#include <numerics/geometry3d/creator/GbPolygon3DCreator.h>

using namespace std;

ObObjectCreator* GbPolygon3D::getCreator()
{
   return GbPolygon3DCreator::getInstance();
}

int GbPolygon3D::counter = 0;

GbPolygon3D::GbPolygon3D()
{
   init();
   counter++;
   this->ps = new GbSystem3D::PointSet3(0);
}
void GbPolygon3D::init()
{
   x1s        = 0.0;
   x2s        = 0.0;
   x1min      = 0.0;
   x1max      = 0.0;
   x2min      = 0.0;
   x2max      = 0.0;
   //		points   = NULL;
   consistent = false;
   ps         = NULL;
}

/**
* Creates an empty 3D polygon with the specified capacity.
* @param capacity the initial capacity
*/
GbPolygon3D::GbPolygon3D(int capacity)
{
   init();
   counter++;
   this->ps = new GbSystem3D::PointSet3(capacity);
   //     this.po = new PointObserver(this);
}
/**
* Creates a 3D polygon with the specified points.
* @param points the initial points of the polygon
*/
GbPolygon3D::GbPolygon3D(vector<GbPoint3D>& points)
{
   init();
   counter++;
   this->ps = new GbSystem3D::PointSet3((int)points.size());
   this->addPoints(points);
}
/**
* Creates a 3D polygon as clone of the specified 3D polygon.
* @param polygon the 3D polygon to be cloned
*/
GbPolygon3D::GbPolygon3D(GbPolygon3D* polygon)
{
   this->init();
   counter++;
   this->ps = new GbSystem3D::PointSet3((int)polygon->size());
   vector<GbPoint3D> temp = polygon->getPoints();
   this->addPoints( temp  );
}

GbPolygon3D::~GbPolygon3D()
{
   counter--;
   //if(points)
   //for(unsigned u=0; u<points->size(); u++)
   //{
   //	delete (*points)[u];
   //}
   //		delete this->points;
   delete this->ps;
}

/*======================================================================*/
/**
* Returns the number of points.
* @return the number of points
*/
int GbPolygon3D::size()
{
   return(this->ps->size());
}
/**
* Returns the number of times this 3D polygon contains the specified point.
* @param point the point
* @return the number of times this 3D polygon contains the specified point
*/
int GbPolygon3D::contains(GbPoint3D* point)
{
   return(this->ps->contains(point));
}
/**
* Returns the number of times this 3D polygon contains a point equal to the specified point.
* @param point the point
* @return the number of times this 3D polygon contains a point equal to the specified point
*/
int GbPolygon3D::containsEqual(GbPoint3D* point)
{
   return(this->ps->containsEqual(point));
}
/**
* Returns true, if this 3D polygon contains the specified line.
* @param point1 the first point
* @param point2 the second point
* @return true, if this 3D polygon contains the specified line
*/
bool GbPolygon3D::containsLine(GbPoint3D* point1, GbPoint3D* point2)
{
   return(this->ps->containsLine(point1, point2));
}
/**
* Returns true, if this 3D polygon contains the specified line.
* @param line the line
* @return true, if this 3D polygon contains the specified line
*/
bool GbPolygon3D::containsLine(GbLine3D* line)
{
   return(this->ps->containsLine(line->getPoint1(), line->getPoint2()));
}
/**
* Returns the first point.
* @return the first point
*/
GbPoint3D* GbPolygon3D::getFirstPoint()
{
   return(this->ps->getFirstPoint());
}
/**
* Returns the last point.
* @return the last point
*/
GbPoint3D* GbPolygon3D::getLastPoint()
{
   return(this->ps->getLastPoint());
}
/**
* Returns the specified point.
* @param index the index
* @return the specified point
* @exception ArrayIndexOutOfBoundsException if the specified index is not valid
*/
GbPoint3D* GbPolygon3D::getPoint(const int& index) 
{
   if(index < 0 || index > this->ps->size()) throw UbException(UB_EXARGS,"ArrayIndexOutOfBoundsException-GbPolygon3D.getPoint()");
   return(this->ps->getPoint(index));
}
/**
* Returns the points.
* @return the points
*/
vector<GbPoint3D> GbPolygon3D::getPoints()
{
   if(!this->consistent) this->calculateValues();
   return(this->points);
}
/**
* Returns the points within the specified rectangle.
* @param p1 the 1st point of the rectangle
* @param p2 the 2nd point of the rectangle
* @return the points within the specified rectangle
*/
vector<GbPoint3D> GbPolygon3D::getPoints(GbPoint3D* p1, GbPoint3D* p2)
{
   return(this->getPoints(p1->x1, p1->x2, p1->x3, p2->x1, p2->x2, p2->x3));
}
/**
* Returns the points within the specified rectangle.
* @param p1x1 the 1st x1 coordinate of the rectangle
* @param p1x2 the 1st x2 coordinate of the rectangle
* @param p1x3 the 1st x3 coordinate of the rectangle
* @param p2x1 the 2nd x1 coordinate of the rectangle
* @param p2x2 the 2nd x2 coordinate of the rectangle
* @param p2x3 the 2nd x3 coordinate of the rectangle
* @return the points within the specified rectangle
*/
vector<GbPoint3D> GbPolygon3D::getPoints(const double& p1x1, const double& p1x2, const double& p1x3, const double& p2x1, const double& p2x2, const double& p2x3)
{
   double x1min, x1max, x2min, x2max, x3min, x3max;

   if(UbMath::less(p1x1, p2x1)) { x1min = p1x1; x1max = p2x1; }
   else                           { x1min = p2x1; x1max = p1x1; }
   if(UbMath::less(p1x2, p2x2)) { x2min = p1x2; x2max = p2x2; }
   else                           { x2min = p2x2; x2max = p1x2; }
   if(UbMath::less(p1x3, p2x3)) { x3min = p1x3; x3max = p2x3; }
   else                           { x3min = p2x3; x3max = p1x3; }

   GbSystem3D::PointSet3 *pts = new GbSystem3D::PointSet3(1);

   if(!this->consistent) this->calculateValues();
   for(int i=this->size()-1; i>=0; i--)
   {
      if(UbMath::lessEqual(x1min, (this->points)[i].x1) && UbMath::greaterEqual(x1max, (this->points)[i].x1) &&
         UbMath::lessEqual(x2min, (this->points)[i].x2) && UbMath::greaterEqual(x2max, (this->points)[i].x2) &&
         UbMath::lessEqual(x3min, (this->points)[i].x3) && UbMath::greaterEqual(x3max, (this->points)[i].x3))     pts->add((this->points)[i]);
   }
   return(pts->getPoints());
}
/**
* Returns the area of this polygon.
* The area is positive for positive ordered points, otherwise negative.
* @return the area of this polygon
*/
//double getArea()
//{
//   if(!this.consistent) this.calculateValues();
//   return(this.area);
//}
double GbPolygon3D::getX1Centroid()
{
   if(!this->consistent) this->calculateValues();
   return(this->x1s);
}
double GbPolygon3D::getX1Minimum()
{
   if(!this->consistent) this->calculateValues();
   return(this->x1min);
}
double GbPolygon3D::getX1Maximum()
{
   if(!this->consistent) this->calculateValues();
   return(this->x1max);
}
double GbPolygon3D::getX2Centroid()
{
   if(!this->consistent) this->calculateValues();
   return(this->x2s);
}
double GbPolygon3D::getX2Minimum()
{
   if(!this->consistent) this->calculateValues();
   return(this->x2min);
}
double GbPolygon3D::getX2Maximum()
{
   if(!this->consistent) this->calculateValues();
   return(this->x2max);
}
double GbPolygon3D::getX3Centroid()
{
   if(!this->consistent) this->calculateValues();
   return(this->x3s);
}
double GbPolygon3D::getX3Minimum()
{
   if(!this->consistent) this->calculateValues();
   return(this->x3min);
}
double GbPolygon3D::getX3Maximum()
{
   if(!this->consistent) this->calculateValues();
   return(this->x3max);
}

/**
* Adds a point to the end of this polygon. Notifies the observers of this 3D polygon.
* @param point the point
*/
void GbPolygon3D::addPoint(GbPoint3D* point)
{
   //if((this instanceof GbPolygon3D) && !(point instanceof GbPoint3D)) throw new IllegalArgumentException("GbPolygon3D.addPoint(): points of 3D polygons have to be 3D points!");

   this->ps->add(point);
   //point.addObserver(this.po);
   this->consistent = false;
   //super.notifyObservers();
}
/**
* Adds a number of points to the end of this polygon. Notifies the observers of this 3D polygon.
* @param points the points
*/
void GbPolygon3D::addPoints(vector<GbPoint3D>& points)
{
   //if((this instanceof GbPolygon3D) && (points.getClass().getComponentType() != GbPoint3D.class)) throw new IllegalArgumentException("GbPolygon3D.addPoints(): points of 3D polygons have to be 3D points!");

   this->ps->add(points);
   //for(int i=0; i<points.length; i++) points[i].addObserver(this.po);
   this->consistent = false;
   //super.notifyObservers();
}
/**
* Inserts a point at the specified position within this polygon. Notifies the observers of this 3D polygon.
* @param point the point
* @param index the index
* @exception ArrayIndexOutOfBoundsException if the specified index is not valid
*/
//public void insertPoint(GbPoint3D point, int index) throws ArrayIndexOutOfBoundsException
//{
//   if((this instanceof GbPolygon3D) && !(point instanceof GbPoint3D)) throw new IllegalArgumentException("GbPolygon3D.insertPoint(): points of 3D polygons have to be 3D points!");
//   if(index < 0 || index > this.ps.size()) throw new ArrayIndexOutOfBoundsException("GbPolygon3D.insert(): invalid index specified: "+index);

//   this.ps.insert(point, index);
//   point.addObserver(this.po);
//   this.consistent = false;
//   super.notifyObservers();
//}
/**
* Removes all points from this polygon identical to the specified one. Notifies the observers of this 3D polygon.
* @param point the point
*/
//public void deletePoint(GbPoint3D point)
//{
//   this.ps.delete(point);
//   point.removeObserver(this.po);
//   this.consistent = false;
//   super.notifyObservers();
//}
/**
* Removes all points from this polygon equal to the specified one. Notifies the observers of this 3D polygon.
* @param point the point
*/
//public void deleteEqualPoint(GbPoint3D point)
//{
//   this.ps.deleteEqual(point);
//   point.removeObserver(this.po);
//   this.consistent = false;
//   super.notifyObservers();
//}
/**
* Removes all points from this polygon. Notifies the observers of this 3D polygon.
*/
void GbPolygon3D::clear()
{
   //		delete this->points;
   this->ps->clearAndTrim();
   delete this->ps;

   //for(int i=points.length-1; i>=0; i--) points[i].removeObserver(this.po);
   this->consistent = false;
   //super.notifyObservers();
}

/**
* Returns true if this 3D polygon equals the specified object.
* Two polygon are equal, if their points are equal.
* <BR>Note that the order of points is recognized!
* @return true if this 3D polygon equals the specified object
* @see GbPoint3D#equals(java.lang.Object)
*/
// bool equals(Object object)
// {
//    try
//    {
//	GbPolygon2D polygon = (GbPolygon3D) object;
//int         n       = this.size();

//if(n != polygon.size()) return(false);
//for(int i=0; i<n; i++) if(!this.getPoint(i).equals(polygon.getPoint(i))) return(false);
//return(true);
//    }
//    catch(Exception e){ return(false); }
// }
/**
* Returns a string representation of this 3D polygon.
* @return a string representation of this 3D polygon
*/
string GbPolygon3D::toString()
{
   stringstream ss;
   ss<<"GbPolygon3D[";
   ss<<this->size()<<" points";

   //    ss<<", x1s="<<this->x1s;
   //    ss<<", x2s="<<this->x2s;
   //ss<<", x3s="<<this->x3s;
   //    ss<<", x1min="<<this->x1min;
   //    ss<<", x1max="<<this->x1max;
   //    ss<<", x2min="<<this->x2min;
   //    ss<<", x2max="<<this->x2max;
   //ss<<", x3min="<<this->x3min;
   //    ss<<", x3max="<<this->x3max;
   ss<<"]"<<endl;
   for(int u=0; u<this->size(); u++)
      ss<<this->ps->getPoint(u)->toString()<<endl;

   return(ss.str());
}
/*======================================================================*/


/*======================================================================*/
/*  Calculation                                                         */
/*                                                                      */
/*
* Returns the intersection points of this 3D polygon and the specified 3D line.
* @param line the 3D line to intersect
* @return the intersection points of this 3D polygon and the specified 3D line
*/
// public GbPoint3D[] calculateIntersectionPoints3D(GbLine3D line)
// {
//    GbSystem.PointSet pointSet = new GbSystem.PointSet(0);
//    GbPoint3D         points[] = this.getPoints();
//    GbPoint3D         pCrossed = null;
//    int               n        = points.length;
//    if(n < 2)         return(pointSet.getPoints());

//    for(int i=1; i<n; i++)
//    {
//pCrossed = GbSystem.calculateIntersectionPoint3D(points[i-1], points[i], line.p1, line.p2);
//if(pCrossed != null) pointSet.add(pCrossed);
//    }
//    pCrossed = GbSystem.calculateIntersectionPoint3D(points[n-1], points[0], line.p1, line.p2);
//    if(pCrossed != null) pointSet.add(pCrossed);

//    return(pointSet.getPoints());
// }

/*
* Returns true if the specified 3D point lies within (or on the border of) this 3D polygon.
* @param point the 3D point to check
* @return true if the specified 3D point lies within (or on the border of) this 3D polygon
*/
// public boolean enclosesPoint3D(GbPoint3D point)
// {
//    if(GbSystem.less(point.x1, this.x1min))    return(false);
//    if(GbSystem.less(point.x2, this.x2min))    return(false);
//    if(GbSystem.greater(point.x1, this.x1max)) return(false);
//    if(GbSystem.greater(point.x2, this.x2max)) return(false);
//    if(this.containsEqual(point) > 0)          return(true);

//    QbList    ltest    = new QbList(GbPoint2D.class, QbList.NOEQUALOBJECTS);
//    GbPoint3D points[] = this.getPoints();
//    GbPoint3D ptest;
//    int       n        = points.length;
//    if(n < 2) return(false);

//    if(GbSystem.equal(point.x2, this.x2min)) ptest = new GbPoint3D(point.x1, this.x2min-1.0);
//    else                                     ptest = new GbPoint3D(point.x1, this.x2max+1.0);

//    for(int i=1; i<n; i++)
//    {
//try { ltest.append(GbSystem.calculateIntersectionPoint2D(points[i-1], points[i], point, ptest)); }
//catch(Exception e){}
//    }
//    try { ltest.append(GbSystem.calculateIntersectionPoint3D(points[n-1], points[0], point, ptest)); }
//    catch(Exception e){}
//    return((ltest.size()%2)==1);
// }

/*
* Returns a new 3D polygon clipped by the specified 3D rectangle (result may be null!).
* @param rectangle the 3D rectangle
* @return a new 3D polygon clipped by the specified 3D rectangle
*/
// GbPolygon3D *createClippedPolygon3D(GbCuboid3D *cube)
// {
//return(GbSystem::clipPolygon3D(this->getPoints(), cube->p1->x1, cube->p1->x2, cube->p1->x3, , cube->p2->x1, cube->p2->x2, cube->p2->x3));
// }
/*                                          
* Returns a new 3D polygon clipped by the specified 3D rectangle (result may be null!).
* @param p1 the 1st point of the rectangle
* @param p2 the 2nd point of the rectangle
* @return a new 3D polygon clipped by the specified 3D rectangle
*/
// GbPolygon3D *createClippedPolygon3D(GbPoint3D *p1, GbPoint3D *p2)
// {
//return(GbSystem::clipPolygon3D(this->getPoints(), p1->x1, p1->x2, p1->x3, p2->x1, p2->x2, p2->x3));
// }
/*
* Returns a new 3D polygon clipped by the specified 3D rectangle (result may be null!).
* @param p1x1 the 1st x1 coordinate of the rectangle
* @param p1x2 the 1st x2 coordinate of the rectangle
* @param p2x1 the 2nd x1 coordinate of the rectangle
* @param p2x2 the 2nd x2 coordinate of the rectangle
* @return a new 3D polygon clipped by the specified 3D rectangle
*/
// GbPolygon3D *createClippedPolygon3D(double p1x1, double p1x2, double p1x3, double p2x1, double p2x2, double p2x3)
// {
//return(GbSystem::clipPolygon3D(this.getPoints(), p1x1, p1x2, p1x3, p2x1, p2x2. p2x3));
// }

/*
* Returns true if the specified 3D rectangle lies completely within this 3D polygon.
* @param rectangle the 3D rectangle to check
* @return true if the specified 3D rectangle lies completely within this 3D polygon
*/
//public boolean enclosesRectangle3D(GbRectangle3D rectangle)
//{
//   GbPolygon3D p = GbSystem.clipPolygon3D(this.getPoints(), rectangle.p1.x1, rectangle.p1.x2, rectangle.p2.x1, rectangle.p2.x2);
//   return(p!=null && GbSystem.equal(Math.abs(p.getArea()), rectangle.getArea()));
//}
/*
* Returns true if the specified 3D rectangle lies completely within this 3D polygon.
* @param p1 the 1st point of the rectangle to check
* @param p2 the 2nd point of the rectangle to check
* @return true if the specified 3D rectangle lies completely within this 3D polygon
*/
//public boolean enclosesRectangle3D(GbPoint3D p1, GbPoint3D p2)
//{
//   GbPolygon3D p = GbSystem.clipPolygon3D(this.getPoints(), p1.x1, p1.x2, p2.x1, p2.x2);
//   return(p!=null && GbSystem.equal(Math.abs(p.getArea()), Math.abs((p1.x1-p2.x1)*(p1.x2-p2.x2))));
//}
/*
* Returns true if the specified 3D rectangle lies completely within this 3D polygon.
* @param p1x1 the 1st x1 coordinate of the rectangle to check
* @param p1x2 the 1st x2 coordinate of the rectangle to check
* @param p2x1 the 2nd x1 coordinate of the rectangle to check
* @param p2x2 the 2nd x2 coordinate of the rectangle to check
* @return true if the specified 3D rectangle lies completely within this 3D polygon
*/
//public boolean enclosesRectangle3D(double p1x1, double p1x2, double p2x1, double p2x2)
//{
//   GbPolygon3D p = GbSystem.clipPolygon3D(this.getPoints(), p1x1, p1x2, p2x1, p2x2);
//   return(p!=null && GbSystem.equal(Math.abs(p.getArea()), Math.abs((p1x1-p2x1)*(p1x2-p2x2))));
//}

/*
* Returns true if the specified 3D rectangle is crossed by this 3D polygon.
* @param rectangle the 3D rectangle to check
* @return true if the specified 3D rectangle is crossed by this 3D polygon
*/
//public boolean crossesRectangle3D(GbRectangle3D rectangle)
//{
//   GbPolygon3D p = GbSystem.clipPolygon3D(this.getPoints(), rectangle.p1.x1, rectangle.p1.x2, rectangle.p2.x1, rectangle.p2.x2);
//   return(p!=null && GbSystem.inOpenInterval(Math.abs(p.getArea()), 0.0, rectangle.getArea()));
//}
/*
* Returns true if the specified 3D rectangle is crossed by this 3D polygon.
* @param p1 the 1st point of the rectangle to check
* @param p2 the 2nd point of the rectangle to check
* @return true if the specified 3D rectangle is crossed by this 3D polygon
*/
//public boolean crossesRectangle3D(GbPoint3D p1, GbPoint3D p2)
//{
//   GbPolygon3D p = GbSystem.clipPolygon3D(this.getPoints(), p1.x1, p1.x2, p2.x1, p2.x2);
//   return(p!=null && GbSystem.inOpenInterval(Math.abs(p.getArea()), 0.0, Math.abs((p1.x1-p2.x1)*(p1.x2-p2.x2))));
//}
/*
* Returns true if the specified 3D rectangle is crossed by this 3D polygon.
* @param p1x1 the 1st x1 coordinate of the rectangle to check
* @param p1x2 the 1st x2 coordinate of the rectangle to check
* @param p2x1 the 2nd x1 coordinate of the rectangle to check
* @param p2x2 the 2nd x2 coordinate of the rectangle to check
* @return true if the specified 3D rectangle is crossed by this 3D polygon
*/
//public boolean crossesRectangle3D(double p1x1, double p1x2, double p2x1, double p2x2)
//{
//   GbPolygon3D p = GbSystem.clipPolygon3D(this.getPoints(), p1x1, p1x2, p2x1, p2x2);
//   return(p!=null && GbSystem.inOpenInterval(Math.abs(p.getArea()), 0.0, Math.abs((p1x1-p2x1)*(p1x2-p2x2))));
//}

/*
* Returns true if the specified 3D rectangle lies (at least partly) within this 3D polygon.
* @param rectangle the 3D rectangle to check
* @return true if the specified 3D rectangle lies (at least partly) within this 3D polygon
*/
//public boolean enclosesOrCrossesRectangle3D(GbRectangle3D rectangle)
//{
//   GbPolygon3D p = GbSystem.clipPolygon3D(this.getPoints(), rectangle.p1.x1, rectangle.p1.x2, rectangle.p2.x1, rectangle.p2.x2);
//   return(p!=null && GbSystem.greater(Math.abs(p.getArea()), 0.0));
//}
/*
* Returns true if the specified 3D rectangle lies (at least partly) within this 3D polygon.
* @param p1 the 1st point of the rectangle to check
* @param p2 the 2nd point of the rectangle to check
* @return true if the specified 3D rectangle lies (at least partly) within this 3D polygon
*/
//public boolean enclosesOrCrossesRectangle3D(GbPoint3D p1, GbPoint3D p2)
//{
//   GbPolygon3D p = GbSystem.clipPolygon3D(this.getPoints(), p1.x1, p1.x2, p2.x1, p2.x2);
//   return(p!=null && GbSystem.greater(Math.abs(p.getArea()), 0.0));
//}
/*
* Returns true if the specified 3D rectangle lies (at least partly) within this 3D polygon.
* @param p1x1 the 1st x1 coordinate of the rectangle to check
* @param p1x2 the 1st x2 coordinate of the rectangle to check
* @param p2x1 the 2nd x1 coordinate of the rectangle to check
* @param p2x2 the 2nd x2 coordinate of the rectangle to check
* @return true if the specified 3D rectangle lies (at least partly) within this 3D polygon
*/
//public boolean enclosesOrCrossesRectangle3D(double p1x1, double p1x2, double p2x1, double p2x2)
//{
//   GbPolygon3D p = GbSystem.clipPolygon3D(this.getPoints(), p1x1, p1x2, p2x1, p2x2);
//   return(p!=null && GbSystem.greater(Math.abs(p.getArea()), 0.0));
//}
/*======================================================================*/

void GbPolygon3D::calculateValues()
{
   this->x1s        = 0.0;
   this->x2s        = 0.0;
   this->x3s        = 0.0;
   this->x1min      = 0.0;
   this->x1max      = 0.0;
   this->x2min      = 0.0;
   this->x2max      = 0.0;
   this->x3min      = 0.0;
   this->x3max      = 0.0;
   throw UbException(UB_EXARGS,"should be implemented");

   //this->consistent = true;

   //this->points = this->ps->getPoints();

   //int       n     = (int)this->points.size();
   //if(n < 1) return;

   //GbPoint3D p1 = (this->points)[0];
   //GbPoint3D p2 = NULL;
   //double    h1 = 0.0;
   //double    h2 = 0.0;
   //double    f=0.0;

   //this->x1s   = p1.x1;
   //this->x1min = p1.x1;
   //this->x1max = p1.x1;
   //this->x2s   = p1.x2;
   //this->x2min = p1.x2;
   //this->x2max = p1.x2;
   //this->x3s   = p1.x2;
   //this->x3min = p1.x2;
   //this->x3max = p1.x2;

   //std::cout<<"Should be implemented "<<endl;

   //for(int i=1; i<n; i++)
   //{
   //  p2         = (this->points)[i];
   //  f          = p1.x1*p2.x2 - p1.x2*p2.x1;
   //  this->area += f;
   //  h1        += f*(p1.x2 + p2.x2);
   //  h2        += f*(p1.x1 + p2.x1);
   //  p1         = p2;

   //  if(p1.x1 < this->x1min) this->x1min = p1.x1;
   //  if(p1.x1 > this->x1max) this->x1max = p1.x1;
   //  if(p1.x2 < this->x2min) this->x2min = p1.x2;
   //  if(p1.x2 > this->x2max) this->x2max = p1.x2;
   //}
   //p2         = (this->points)[0];
   //f          = p1.x1*p2.x2 - p1.x2*p2.x1;
   //this->area += f;
   //h1        += f*(p1.x2 + p2.x2);
   //h2        += f*(p1.x1 + p2.x1);

   //this->area *= 0.5;
   //h1        /= 6.0;
   //h2        /= 6.0;

   //if(n > 2)
   //{
   //   this->x1s = h2/this->area;
   //   this->x2s = h1/this->area;
   //}

   //if(n < 3 || !GbSystem::inClosedInterval(this->x1s, this->x1min, this->x1max)) this->x1s = 0.5*(this->x1min+this->x1max);
   //if(n < 3 || !GbSystem::inClosedInterval(this->x2s, this->x2min, this->x2max)) this->x2s = 0.5*(this->x2min+this->x2max);
}
/*======================================================================*/


/*======================================================================*/
// private class PointObserver implements TiObserver
// {
//    GbPolygon3D polygon = null;

//    PointObserver(GbPolygon3D polygon)
//    {
//this.polygon = polygon;
//    }

//    public void objectChanged(Object object)
//    {
//if((object instanceof GbPoint3D) && this.polygon.contains((GbPoint3D)object)>0)
//{
//   this.polygon.consistent = false;
//   this.polygon.notifyObservers();
//}
//    }
// }
/*=======================================================*/
void GbPolygon3D::write(UbFileOutput* out) 
{                                      
   throw UbException(UB_EXARGS,"not implemented");
}
/*=======================================================*/
void GbPolygon3D::read(UbFileInput* in) 
{  
   throw UbException(UB_EXARGS,"not implemented");
}
/*=======================================================*/


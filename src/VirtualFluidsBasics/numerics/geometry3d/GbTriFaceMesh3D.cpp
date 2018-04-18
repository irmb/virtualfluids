#include <numerics/geometry3d/GbTriFaceMesh3D.h>

#include <numerics/geometry3d/GbCuboid3D.h>
#include <numerics/geometry3d/GbHalfSpace3D.h>
#include <numerics/geometry3d/CoordinateTransformation3D.h>
#include <numerics/geometry3d/creator/GbTriFaceMesh3DCreator.h>
#include <basics/utilities/UbRandom.h>
#include <basics/utilities/UbTiming.h>
#include <basics/utilities/UbLogger.h>

#include <numerics/geometry3d/KdTree/KdTree.h>
#include <numerics/geometry3d/KdTree/splitalgorithms/KdSpatiallMedianSplit.h>
#include <numerics/geometry3d/KdTree/splitalgorithms/KdSAHSplit.h>
#include <numerics/geometry3d/KdTree/intersectionhandler/KdCountLineIntersectionHandler.h>
#include <numerics/geometry3d/KdTree/intersectionhandler/KdCountRayIntersectionHandler.h>

using namespace std;

GbTriFaceMesh3D::GbTriFaceMesh3D() 
   :   GbObject3D()
     , buildVertTriRelationMap(false)
     , kdTree(NULL)
     , kdTreeValid(false)
{
   this->setName("CAB_GbTriFaceMesh3D");
   this->nodes          = new vector<Vertex>;
   this->triangles      = new vector<TriFace>;
   this->consistent     = false;
   //this->kdtreeSplitAlg = KDTREE_SHAPLIT;
   this->kdtreeSplitAlg = KDTREE_SPATIALSPLIT;
}
/*=======================================================================*/
GbTriFaceMesh3D::GbTriFaceMesh3D(  string name, vector<Vertex>* nodes, vector<TriFace>* triangles, KDTREE_SPLITAGORITHM splitAlg, bool removeRedundantNodes)
   :  GbObject3D()
    , nodes(nodes)
    , triangles(triangles)
    , buildVertTriRelationMap(false)
    , consistent(false)
    , kdTree(NULL)
    , kdTreeValid(false)
    , kdtreeSplitAlg(splitAlg)
{
   if( name.empty() ) throw UbException(UB_EXARGS,"no name specified");
   if( !nodes       ) throw UbException(UB_EXARGS,"no nodes specified");
   if( !triangles   ) throw UbException(UB_EXARGS,"no triangles specified");

   this->setName(name);

   if(removeRedundantNodes)
   {
      this->deleteRedundantNodes(); //dort wird autoamtisch calculateValues() aufgerufen
   }
   else
   {
      this->calculateValues();
   }
}
/*=======================================================================*/
GbTriFaceMesh3D::~GbTriFaceMesh3D()
{
   UBLOG(logERROR,"~GbTriFaceMesh3D")
   if( nodes     ) { delete nodes;     nodes     = NULL; }
   if( triangles ) { delete triangles; triangles = NULL; }
   if( kdTree    ) { delete kdTree;    kdTree    = NULL; }
}
/*======================================================================*/
ObObjectCreator* GbTriFaceMesh3D::getCreator()
{
   return GbTriFaceMesh3DCreator::getInstance();
}
/*======================================================================*/
void GbTriFaceMesh3D::init()
{
   nodes      = NULL;
   triangles  = NULL;
   x1min      = 0.0;
   x1max      = 0.0;
   x1center   = 0.0;
   x2min      = 0.0;
   x2max      = 0.0;
   x2center   = 0.0;
   x3min      = 0.0;
   x3max      = 0.0;
   x3center   = 0.0;
   consistent = false;
   kdTreeValid = false;
}
/*======================================================================*/
//checks for doppelt nodes und fixed Dreicke die zweimal denselben Knoten haben
void GbTriFaceMesh3D::deleteRedundantNodes()
{
    UBLOG(logINFO,"GbTriFaceMesh3D::deleteRedundantNodes - Nodes before deleting redundant: "<<this->nodes->size());

    map<Vertex,size_t/*new vecIndex*/> vertexMap;
    map<Vertex,size_t/*new vecIndex*/>::iterator pos;
    
    vector<TriFace>& tris     = *this->triangles;
    vector<Vertex>&  oldNodes = *this->nodes;
    vector<Vertex>   newNodes;

    for(size_t t=0; t<tris.size(); t++)
    {
       TriFace& tri = tris[t];
       //Knoten bereits in neuem node vector?
       for(int v=0; v<=2; v++)
       {
          Vertex& vert = tri.getNode(v,oldNodes);
          pos=vertexMap.find( vert );
          if( pos!=vertexMap.end() ) tri.setNode(v, (int)pos->second);
          else
          {
             newNodes.push_back(vert);
             int index = (int)newNodes.size()-1;
             vertexMap[vert] = index;                       
             tri.setNode(v,index);
          }
       }
    }

    std::swap(*nodes,newNodes);

    UBLOG(logINFO,"GbTriFaceMesh3D::deleteRedundantNodes - Nodes after deleting redundant:"<<this->nodes->size());
    UBLOG(logINFO,"GbTriFaceMesh3D::deleteRedundantNodes - checking for triangles that have same node several times or are lines!!!");
    int counter1 = 0;
    int counter2 = 0;
    vector<TriFace> newTris;
    newTris.reserve( this->triangles->size() );
    for(size_t t=0; t<tris.size(); t++)
    {
       Vertex& v1 = tris[t].getNode(0,*nodes); 
       Vertex& v2 = tris[t].getNode(1,*nodes); 
       Vertex& v3 = tris[t].getNode(2,*nodes); 
       if( v1==v2 || v1==v3 || v2==v3 )
       {
          counter1++;
       }
       else if( tris[t].getArea(*nodes)<1.0E-8 )
       {
          counter2++;
       }
       else newTris.push_back(tris[t]);
    }
    if(counter1)
    {
       UBLOG(logWARNING,"GbTriFaceMesh3D::deleteRedundantNodes - ### Warning ###: found and removed  "<<counter1<<" triangle with double nodes!");
    }
    if(counter2)
    {
       UBLOG(logWARNING,"GbTriFaceMesh3D::deleteRedundantNodes - ### Warning ###: found and removed  "<<counter2<<" triangle that are lines!") ;
    }
    if(!counter1 && !counter2) { UBLOG(logINFO,"GbTriFaceMesh3D::deleteRedundantNodes - alles gut... nix doppelt"); }
    else swap(tris,newTris);
    
    UBLOG(logINFO,"GbTriFaceMesh3D::deleteRedundantNodes - done" );
    this->calculateValues();

    UBLOG(logINFO,"GbTriFaceMesh3D::deleteRedundantNodes !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
    UBLOG(logINFO,"GbTriFaceMesh3D::deleteRedundantNodes !!! Falls das FeTriFaceMesh verwendet wird !!!");
    UBLOG(logINFO,"GbTriFaceMesh3D::deleteRedundantNodes !!! und darin die VertexTriFaceMap         !!!");
    UBLOG(logINFO,"GbTriFaceMesh3D::deleteRedundantNodes !!! dann muss diese neu erzeugt werden     !!!");
    UBLOG(logINFO,"GbTriFaceMesh3D::deleteRedundantNodes !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
}
/*======================================================================*/
   /**
    * Returns a string representation of this triangular mesh.
    * @return a string representation of this triangular mesh
    */
string GbTriFaceMesh3D::toString()
{
	stringstream ss;
	ss<<"GbTriFaceMesh3D[";
	ss<<(int)this->triangles->size()<<"-Triangles, "<<(int)this->nodes->size()<<"-Nodes, "<<endl;
   ss<<"]";
   return(ss.str());
}
/**
 * Returns the nodes of this triangular mesh.
 * @return the nodes of this triangular mesh
 */
vector<GbTriFaceMesh3D::Vertex>* GbTriFaceMesh3D::getNodes()       {  return this->nodes;   }
/**
 * Returns the triangles of this triangular mesh.
 * @return the triangles of this triangular mesh
 */
vector<GbTriFaceMesh3D::TriFace>* GbTriFaceMesh3D::getTriangles()  { return this->triangles;  }
/*===============================================*/
double GbTriFaceMesh3D::getVolume()
{
   vector<Vertex>&  vertices = *nodes;
   vector<TriFace>& tris     = *triangles;

   double x1,x2,x3,y1,y2,y3,z1,z2,z3, G3i;
   //double rSP1 = 0.0;double rSP2 = 0.0;double rSP3 = 0.0;
   double volume = 0.0;
   for(size_t t=0; t<tris.size(); t++)
   {
      TriFace& triangle = tris[t];
      x1     = triangle.getV1x(vertices); y1 = triangle.getV1y(vertices); z1 = triangle.getV1z(vertices);
      x2     = triangle.getV2x(vertices); y2 = triangle.getV2y(vertices); z2 = triangle.getV2z(vertices);
      x3     = triangle.getV3x(vertices); y3 = triangle.getV3y(vertices); z3 = triangle.getV3z(vertices);
      G3i    = x1*(y2*z3-z2*y3)+y1*(z2*x3-x2*z3)+z1*(x2*y3-y2*x3);
      volume = volume+G3i/6.0;
   }
   return volume;
}
/*===============================================*/
UbTupleDouble3 GbTriFaceMesh3D::calculateCenterOfGravity()
{
   vector<Vertex>&  vertices = *nodes;
   vector<TriFace>& tris     = *triangles;
   
   double x1,x2,x3,y1,y2,y3,z1,z2,z3;
   double G3i;
   double rSP1 = 0.0, rSP2 = 0.0, rSP3 = 0.0, volume = 0.0;
   
   for(size_t t=0; t<tris.size(); t++)
   {
      TriFace& triangle = tris[t];
      x1     = triangle.getV1x(vertices); y1 = triangle.getV1y(vertices); z1 = triangle.getV1z(vertices);
      x2     = triangle.getV2x(vertices); y2 = triangle.getV2y(vertices); z2 = triangle.getV2z(vertices);
      x3     = triangle.getV3x(vertices); y3 = triangle.getV3y(vertices); z3 = triangle.getV3z(vertices);
      G3i    = x1*(y2*z3-z2*y3)+y1*(z2*x3-x2*z3)+z1*(x2*y3-y2*x3);
      volume = volume+G3i/6.0;
      rSP1   = rSP1+G3i*(x1+x2+x3);
      rSP2   = rSP2+G3i*(y1+y2+y3);
      rSP3   = rSP3+G3i*(z1+z2+z3);
   }
   rSP1 = rSP1/(24.0*volume);
   rSP2 = rSP2/(24.0*volume);
   rSP3 = rSP3/(24.0*volume);

   return UbTupleDouble3(rSP1, rSP2, rSP3);
}
/*===============================================*/
UbTupleDouble6 GbTriFaceMesh3D::calculateMomentOfInertia(double rhoP)
{
   vector<Vertex>& vertices = *nodes;

   double x1,x2,x3,y1,y2,y3,z1,z2,z3;
   double G3i;
   double xx,yy,zz,xy,yz,zx;
   double rSP1 = 0.0;double rSP2 = 0.0;double rSP3 = 0.0;
   double volume = 0.0;
   double top11 = 0.0;double top22 = 0.0;double top33 = 0.0;
   double top12 = 0.0;double top23 = 0.0;double top13 = 0.0;
   int size = (int)this->triangles->size();
   for(int u=0; u<size;u++)
   {
      TriFace& triangle = (*this->triangles)[u];
      x1 = triangle.getV1x(vertices); y1 = triangle.getV1y(vertices); z1 = triangle.getV1z(vertices);
      x2 = triangle.getV2x(vertices); y2 = triangle.getV2y(vertices); z2 = triangle.getV2z(vertices);
      x3 = triangle.getV3x(vertices); y3 = triangle.getV3y(vertices); z3 = triangle.getV3z(vertices);
      G3i = x1*(y2*z3-z2*y3)+y1*(z2*x3-x2*z3)+z1*(x2*y3-y2*x3);
      volume = volume+G3i/6.0;
      rSP1 = rSP1+G3i*(x1+x2+x3);
      rSP2 = rSP2+G3i*(y1+y2+y3);
      rSP3 = rSP3+G3i*(z1+z2+z3);
   }
   rSP1 = rSP1/(24.0*volume);
   rSP2 = rSP2/(24.0*volume);
   rSP3 = rSP3/(24.0*volume);

   double x1s = 0.0;//rSP1;//0.0;//
   double x2s = 0.0;//rSP2;//0.0;//
   double x3s = 0.0;//rSP3;//0.0;//

   for(int u=0; u<size;u++)
   {
      TriFace& triangle = (*this->triangles)[u];
      x1 = triangle.getV1x(vertices)-x1s;
      y1 = triangle.getV1y(vertices)-x2s;
      z1 = triangle.getV1z(vertices)-x3s;
      x2 = triangle.getV2x(vertices)-x1s;
      y2 = triangle.getV2y(vertices)-x2s;
      z2 = triangle.getV2z(vertices)-x3s;
      x3 = triangle.getV3x(vertices)-x1s;
      y3 = triangle.getV3y(vertices)-x2s;
      z3 = triangle.getV3z(vertices)-x3s;
      G3i = x1*(y2*z3-z2*y3)+y1*(z2*x3-x2*z3)+z1*(x2*y3-y2*x3);
      //rSP1 = rSP1+G3i*(x1+x2+x3)/(24.0*volume);
      //rSP2 = rSP2+G3i*(y1+y2+y3)/(24.0*volume);
      //rSP3 = rSP3+G3i*(z1+z2+z3)/(24.0*volume);
      xx = x1*x1+x2*x2+x3*x3+x1*x2+x2*x3+x3*x1;
      yy = y1*y1+y2*y2+y3*y3+y1*y2+y2*y3+y3*y1;
      zz = z1*z1+z2*z2+z3*z3+z1*z2+z2*z3+z3*z1;
      top11 = top11+(yy+zz)*rhoP*G3i/60.;
      top22 = top22+(xx+zz)*rhoP*G3i/60.;
      top33 = top33+(yy+xx)*rhoP*G3i/60.;
      xy = 2.0*(x1*y1+x2*y2+x3*y3)+x2*y3+x3*y1+x1*y2+x3*y2+x1*y3+x2*y1;
      yz = 2.0*(y1*z1+y2*z2+y3*z3)+y2*z3+y3*z1+y1*z2+y3*z2+y1*z3+y2*z1;
      zx = 2.0*(z1*x1+z2*x2+z3*x3)+z2*x3+z3*x1+z1*x2+z3*x2+z1*x3+z2*x1;
      top12 = top12-xy*rhoP*G3i/120.;
      top23 = top23-yz*rhoP*G3i/120.;
      top13 = top13-zx*rhoP*G3i/120.;
   }
   //Satz von Steiner ...
   top11 = top11-rhoP*volume*(rSP2*rSP2+rSP3+rSP3);
   top22 = top22-rhoP*volume*(rSP3*rSP3+rSP1*rSP1);
   top33 = top33-rhoP*volume*(rSP1*rSP1+rSP2*rSP2);
   top12 = top12+rhoP*volume*rSP1*rSP2;
   top23 = top23+rhoP*volume*rSP2*rSP3;
   top13 = top13+rhoP*volume*rSP3*rSP1;

   cout<<"Volume:"<<volume<<"\n Traegheitsmomente:\n";
   cout<<" top11:"<<top11<<" top22:"<<top22<<" top33:"<<top33<<endl;
   cout<<" top12:"<<top12<<" top23:"<<top23<<" top13:"<<top13<<endl;

   return UbTupleDouble6(top11,top22,top33,top12,top23,top13);
}
/*==============================================================*/
void GbTriFaceMesh3D::calculateValues()
{
   relationVertTris.clear();

   if( nodes->empty() )
   {
      x1min = x1max = x2min = x2max = x3min = x3max = 0.0;
   }
   else
   {
      Vertex& v = (*nodes)[0];
      x1min = x1max = v.x;
      x2min = x2max = v.y;
      x3min = x3max = v.z;

      for(size_t i=1; i<this->nodes->size(); i++)
      {
         Vertex& v1 = (*nodes)[i];
         x1min = UbMath::min<double>(x1min,v1.x);
         x1max = UbMath::max<double>(x1max,v1.x);
         x2min = UbMath::min<double>(x2min,v1.y);
         x2max = UbMath::max<double>(x2max,v1.y);
         x3min = UbMath::min<double>(x3min,v1.z);
         x3max = UbMath::max<double>(x3max,v1.z);
      }
      x1center = 0.5*(x1min+x1max);
      x2center = 0.5*(x2min+x2max);
      x3center = 0.5*(x3min+x3max);
      
      vector<TriFace>& tris  = *this->triangles;
      vector<Vertex>&  verts = *this->nodes;
      for(size_t i=0; i<this->triangles->size(); i++)
      {
         tris[i].calculateNormal(verts);
      }
      //relation Vertex <-> Triangle ermitteln
      if(buildVertTriRelationMap)
      {
         for(size_t t=0; t<tris.size(); t++)
         {
            TriFace& tri = tris[t];
            relationVertTris.insert( make_pair( &verts[tri.v1], &tri) );
            relationVertTris.insert( make_pair( &verts[tri.v2], &tri) );
            relationVertTris.insert( make_pair( &verts[tri.v3], &tri) );
         }
      }
   }
   this->consistent = true;

   // Need to rebuild kd-tree.
   kdTreeValid=false;
}

void GbTriFaceMesh3D::generateKdTree()
{
   //////////////////////////////////////////////////////////////////////////
   //ACHTUNG: kdTree MUSS hier erfolgen (nach this->consistent = true; da der kdTree auf Methodes des GbTriFaceMesh zurück greift)
   if (!consistent) calculateValues();

   if( !nodes->empty() )
   {
      UBLOG(logDEBUG1, "GbTriFaceMesh3D::calculateValues - build KdTree start");

      if(kdTree) delete kdTree;

      UbTimer timer; timer.start();
      if(kdtreeSplitAlg == KDTREE_SHAPLIT     ) 
      {
         UBLOG(logDEBUG1, "GbTriFaceMesh3D::calculateValues - build KdTree with SAHSplit");
         this->kdTree = new Kd::Tree<double>( *this, Kd::SAHSplit<double>()            );
      }
      else if(kdtreeSplitAlg == KDTREE_SPATIALSPLIT)
      {
         UBLOG(logDEBUG1, "GbTriFaceMesh3D::calculateValues - build KdTree with SpatialMedianSplit");
         this->kdTree = new Kd::Tree<double>( *this, Kd::SpatialMedianSplit<double>() ); 
      }
      else throw UbException(UB_EXARGS, "unknown kdtree split option)" );
      std::cout<<"GbTriFaceMesh3D::calculateValues - built kdTree in "<<timer.stop()<<"seconds"<<std::endl;
      //UBLOG(logDEBUG1, "GbTriFaceMesh3D::calculateValues - built kdTree in "<<timer.stop()<<"seconds");
   }
   kdTreeValid = true;
}
/*=========================================================================*/
std::vector<GbTriFaceMesh3D::TriFace*> GbTriFaceMesh3D::getTrianglesForVertex(Vertex* vertex)
{
   if(!buildVertTriRelationMap) { buildVertTriRelationMap=true; consistent = false; }
   if(!consistent) this->calculateValues();

   typedef std::multimap<Vertex*,TriFace*>::iterator Iterator;
   pair<Iterator,Iterator> objRange = relationVertTris.equal_range(vertex);

   std::vector<TriFace*> tmpTris;
   for(Iterator pos=objRange.first; pos!=objRange.second; ++pos) 
      tmpTris.push_back( pos->second );

   return tmpTris;
}

/*======================================================================*/
void GbTriFaceMesh3D::transform(const double matrix[4][4])
 {
    if(!this->consistent) this->calculateValues();
    vector<Vertex>& vertices = *nodes;
    float tempX = 0.f;
    float tempY = 0.f;
    float tempZ = 0.f;

    for(size_t i=0; i<vertices.size(); i++)
    {
       Vertex& v = vertices[i];
       tempX = v.x;
       tempY = v.y;
       tempZ = v.z;
       v.x = (float)(matrix[0][0] * tempX + matrix[0][1] * tempY + matrix[0][2] * tempZ + matrix[0][3] * 1.);
       v.y = (float)(matrix[1][0] * tempX + matrix[1][1] * tempY + matrix[1][2] * tempZ + matrix[1][3] * 1.);
       v.z = (float)(matrix[2][0] * tempX + matrix[2][1] * tempY + matrix[2][2] * tempZ + matrix[2][3] * 1.);
    }
    this->calculateValues();
    this->notifyObserversObjectChanged();
 }
/*======================================================================*/
void GbTriFaceMesh3D::rotate(const double& alpha, const double& beta, const double& gamma)
{
   if(!this->consistent) this->calculateValues();
   double a1 = this->getX1Centroid();
   double a2 = this->getX2Centroid();
   double a3 = this->getX3Centroid();
   CoordinateTransformation3D trafoFor(a1, a2, a3, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0);
   CoordinateTransformation3D trafoBack(a1, a2, a3, 1.0, 1.0, 1.0, alpha, beta, gamma);

   vector<Vertex>& vertices = *nodes;
   for(size_t i=0; i<vertices.size(); i++)
   {
      Vertex& v = vertices[i];
      double p1x1 = trafoFor.transformForwardToX1Coordinate(v.x, v.y, v.z);
      double p1x2 = trafoFor.transformForwardToX2Coordinate(v.x, v.y, v.z);
      double p1x3 = trafoFor.transformForwardToX3Coordinate(v.x, v.y, v.z);
      v.x = (float)trafoBack.transformBackwardToX1Coordinate(p1x1, p1x2, p1x3);
      v.y = (float)trafoBack.transformBackwardToX2Coordinate(p1x1, p1x2, p1x3);
      v.z = (float)trafoBack.transformBackwardToX3Coordinate(p1x1, p1x2, p1x3);
   }
   this->calculateValues();
   this->notifyObserversObjectChanged();
}

/*======================================================================*/
void GbTriFaceMesh3D::scale(const double& sx1, const double& sx2, const double& sx3)
{
   if(!this->consistent) this->calculateValues();
   double a1 = this->getX1Centroid();
   double a2 = this->getX2Centroid();
   double a3 = this->getX3Centroid();
   CoordinateTransformation3D trafoFor(a1, a2, a3, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0);
   CoordinateTransformation3D trafoBack(a1, a2, a3, sx1, sx2, sx3, 0.0, 0.0, 0.0);

   vector<Vertex>& vertices = *nodes;
   for(size_t i=0; i<vertices.size(); i++)
   {
      Vertex& v = vertices[i];
      double p1x1 = trafoFor.transformForwardToX1Coordinate(v.x, v.y, v.z);
      double p1x2 = trafoFor.transformForwardToX2Coordinate(v.x, v.y, v.z);
      double p1x3 = trafoFor.transformForwardToX3Coordinate(v.x, v.y, v.z);
      v.x = (float)trafoBack.transformBackwardToX1Coordinate(p1x1, p1x2, p1x3);
      v.y = (float)trafoBack.transformBackwardToX2Coordinate(p1x1, p1x2, p1x3);
      v.z = (float)trafoBack.transformBackwardToX3Coordinate(p1x1, p1x2, p1x3);
   }
   this->calculateValues();
   this->notifyObserversObjectChanged();
}

/*======================================================================*/
void GbTriFaceMesh3D::translate(const double& x1, const double& x2, const double& x3)
{
   vector<Vertex>& vertices = *nodes;
   for(size_t i=0; i<vertices.size(); i++)
   {
      Vertex& v = vertices[i];
      v.x+=static_cast<float>(x1);
      v.y+=static_cast<float>(x2);
      v.z+=static_cast<float>(x3);
   }
   this->calculateValues();
   this->notifyObserversObjectChanged();
}
/*======================================================================*/
vector<GbTriangle3D*> GbTriFaceMesh3D::getSurfaceTriangleSet()
{
   //SirAnn: eine miese Spciehrflochmethode
   //        hier werden dynmamische Objekte angelegt 
   //        mit sowas rechnet von aussen kein Mensch!!!
   vector<GbTriangle3D*> tris( triangles->size() );

   for(size_t i=0; i<this->triangles->size(); i++)
   {
      Vertex& v1 = (*nodes)[(*triangles)[i].v1];
      Vertex& v2 = (*nodes)[(*triangles)[i].v2];
      Vertex& v3 = (*nodes)[(*triangles)[i].v3];

      tris[i] = new GbTriangle3D(  new GbPoint3D(v1.x,v1.y,v1.z)
                                 , new GbPoint3D(v2.x,v2.y,v2.z)
                                 , new GbPoint3D(v3.x,v3.y,v3.z) );
   }
   return tris;
}
/*=======================================================*/
void GbTriFaceMesh3D::addSurfaceTriangleSet(vector<UbTupleFloat3>& pts, vector<UbTupleInt3>& tris)
{
   for(int i=0; i<(int)this->triangles->size(); i++)
   {
      Vertex& v1 = (*nodes)[(*triangles)[i].v1];
      Vertex& v2 = (*nodes)[(*triangles)[i].v2];
      Vertex& v3 = (*nodes)[(*triangles)[i].v3];
      pts.push_back( makeUbTuple(v1.x,v1.y,v1.z));
      pts.push_back( makeUbTuple(v2.x,v2.y,v2.z));
      pts.push_back( makeUbTuple(v3.x,v3.y,v3.z));

      tris.push_back( makeUbTuple( 3*i, 3*i+1, 3*i+2) );
   }
}
/*======================================================================*/
bool GbTriFaceMesh3D::isPointInGbObject3D(const double& x1, const double& x2, const double& x3)
{
   if (!kdTreeValid) generateKdTree();

   int iSec;
   for(int i=0; i<100; i++)
   {
      Kd::Ray<double> ray(  x1, x2, x3  //, 1, 0 ,0 );
                           , ( x1 < x1center ? UbRandom::rand(-1.0,-0.001, 10) : UbRandom::rand(0.001, 1.0, 10) )
                           , ( x2 < x2center ? UbRandom::rand(-1.0,-0.001, 10) : UbRandom::rand(0.001, 1.0, 10) )
                           , ( x3 < x3center ? UbRandom::rand(-1.0,-0.001, 10) : UbRandom::rand(0.001, 1.0, 10) ) );
        
      iSec = kdTree->intersectRay( ray, Kd::CountRayIntersectionHandler<double>() );
        
      if( iSec != Kd::Intersection::INTERSECT_EDGE ) //KEINE Kante getroffen
      {
         if(iSec == Kd::Intersection::ON_BOUNDARY )
         {
            return true;
         }
         return (iSec&1);  //ungerade anzahl an schnitten --> drinnen
      }
      UBLOG(logWARNING, "GbTriFaceMesh3D.isPointInGbObject3D.if  - an edge was hit ");
   }
   throw UbException(UB_EXARGS, "ups, nach 100 Strahlen immer noch kein Ergebnis");
}
/*======================================================================*/
bool GbTriFaceMesh3D::isPointInGbObject3D(const double& x1, const double& x2, const double& x3, bool& pointIsOnBoundary)
{
   if (!kdTreeValid) generateKdTree();
   
   int iSec;
   for(int i=0; i<100; i++)
   {
      Kd::Ray<double> ray(  x1, x2, x3 
                         , float( ( x1 < x1center ? UbRandom::rand(-1.0,-0.001, 10) : UbRandom::rand(0.001, 1.0, 10) ) )
                         , float( ( x2 < x2center ? UbRandom::rand(-1.0,-0.001, 10) : UbRandom::rand(0.001, 1.0, 10) ) )
                         , float( ( x3 < x3center ? UbRandom::rand(-1.0,-0.001, 10) : UbRandom::rand(0.001, 1.0, 10) ) ) );

      iSec = kdTree->intersectRay( ray, Kd::CountRayIntersectionHandler<double>()    );

      if( iSec != Kd::Intersection::INTERSECT_EDGE ) //KEINE Kante getroffen
      {
         if(iSec == Kd::Intersection::ON_BOUNDARY )
         {
            pointIsOnBoundary = true;
            return true;
         }
         pointIsOnBoundary = false;
         return (iSec&1);  //ungerade anzahl an schnitten --> drinnen
      }
   }

   throw UbException(UB_EXARGS, "ups, nach 100 Strahlen immer noch kein Ergebnis");
}
/*======================================================================*/
bool GbTriFaceMesh3D::intersectLine(const double& p1_x1, const double& p1_x2, const double& p1_x3, const double& p2_x1, const double& p2_x2, const double& p2_x3)
{
   if (!kdTreeValid) generateKdTree();

   int iSec = kdTree->intersectLine( UbTupleDouble3(p1_x1, p1_x2, p1_x3), UbTupleDouble3(p2_x1, p2_x2, p2_x3), Kd::CountLineIntersectionHandler<double>() );

   return (iSec != Kd::Intersection::NO_INTERSECTION);
}
/*======================================================================*/
GbLine3D* GbTriFaceMesh3D::createClippedLine3D (GbPoint3D& point1, GbPoint3D& point2)
{
   throw UbException(UB_EXARGS,"not implemented");
}
/*======================================================================*/
void GbTriFaceMesh3D::write(UbFileOutput* out)
{
   out->writeString(this->getCreator()->getTypeID());
   
   out->writeInteger((int)kdtreeSplitAlg);

   //nodes
   vector<Vertex>& vertices = *nodes;
   out->writeSize_t( nodes->size() );
   out->writeLine();
   for(size_t i=0; i<vertices.size(); i++)
   {
      Vertex& v = vertices[i];
      out->writeFloat(v.x);
      out->writeFloat(v.y);
      out->writeFloat(v.z);
      out->writeLine();
   }
   
   //triangles
   vector<TriFace>& tris = *triangles;
   out->writeSize_t( tris.size() );
   out->writeLine();
   for(size_t i=0; i<tris.size(); i++)
   {
      TriFace& t = tris[i];
      out->writeInteger(t.v1);
      out->writeInteger(t.v2);
      out->writeInteger(t.v3);
      out->writeLine();
   }
}
/*======================================================================*/
void GbTriFaceMesh3D::read(UbFileInput* in)
{
   kdtreeSplitAlg =  (KDTREE_SPLITAGORITHM)in->readInteger();

   if(!nodes) nodes = new vector<Vertex>;
   //nodes
   vector<Vertex>& vertices = *nodes;
   vertices.resize( in->readSize_t( ) );
   in->readLine();
   for(size_t i=0; i<vertices.size(); i++)
   {
      Vertex& v = vertices[i];
      v.x = in->readFloat();
      v.y = in->readFloat();
      v.z = in->readFloat();
      in->readLine();
   }

   //triangles
   if(!triangles) triangles = new vector<TriFace>;
   vector<TriFace>& tris = *triangles;
   tris.resize( in->readSize_t( ) );
   in->readLine();
   for(size_t i=0; i<tris.size(); i++)
   {
      TriFace& t = tris[i];
      t.v1 = in->readInteger();
      t.v2 = in->readInteger();
      t.v3 = in->readInteger();
      in->readLine();
   }

   this->calculateValues();
}
/*======================================================================*/
UbTuple<string, string> GbTriFaceMesh3D::writeMesh(string filename, WbWriter* writer, bool writeNormals, vector< string >* datanames, std::vector< std::vector < double > >* nodedata )
{
   UBLOG(logINFO, "GbTriFaceMesh3D::writeMesh ");

   vector<UbTupleFloat3 > triNodes(nodes->size());
   vector<UbTupleInt3 >   tris(triangles->size());

   for(size_t i=0; i<nodes->size(); i++)
      triNodes[i] = makeUbTuple( (*nodes)[i].x, (*nodes)[i].y, (*nodes)[i].z );

   for(size_t i=0; i<triangles->size(); i++)
      tris[i] = makeUbTuple( (*triangles)[i].v1, (*triangles)[i].v2, (*triangles)[i].v3 ) ;

   UbTuple<string, string> filenames("","");

   if( !datanames || datanames->empty() || !nodedata  )
   {
      val<1>(filenames) = writer->writeTriangles(filename,triNodes,tris);
   }
   else
   {
      val<1>(filenames) = writer->writeTrianglesWithNodeData(filename,triNodes,tris,*datanames,*nodedata);
   }

   if(writeNormals)
   {
      vector<UbTupleFloat3 > lineNodes(triangles->size()*2);
      vector<UbTupleInt2 >   lines(triangles->size());
      for(size_t i=0; i<triangles->size(); i++)
      {
         TriFace& triangle = (*triangles)[i];
         lineNodes[i*2  ] = makeUbTuple( triangle.getX1Centroid(*nodes)
                                        ,triangle.getX2Centroid(*nodes)
                                        ,triangle.getX3Centroid(*nodes));
         lineNodes[i*2+1] = makeUbTuple( (float)(triangle.getX1Centroid(*nodes)+1.*triangle.nx)
                                        ,(float)(triangle.getX2Centroid(*nodes)+1.*triangle.ny)
                                        ,(float)(triangle.getX3Centroid(*nodes)+1.*triangle.nz));

         lines[i] = makeUbTuple((int)i*2,(int)i*2+1);
      }
      val<2>(filenames) = writer->writeLines(filename+"_normals",lineNodes,lines);
   }

   return filenames;
}
/*======================================================================*/

#ifndef KDCOUNTRAYINTERSECTIONHANDLER_H
#define KDCOUNTRAYINTERSECTIONHANDLER_H

#include <basics/utilities/UbTuple.h>
#include <basics/utilities/UbKeys.h>

#include <numerics/geometry3d/GbTriFaceMesh3D.h>

#include <numerics/geometry3d/KdTree/KdNode.h>
//#include <numerics/geometry3d/KdTree/KdUtilities.h>
#include <numerics/geometry3d/KdTree/intersectionhandler/KdRayIntersectionHandler.h>

#include <set>

namespace Kd
{
   template< typename T >
   class CountRayIntersectionHandler : public RayIntersectionHandler<T> 
   {
   public:
      int intersectRay(const Ray<T>& ray, Node<T>& parent, Node<T>*& child1, Node<T>*& child2, std::set< UbKeys::Key3<int> >& mailbox) const
      {
         if( parent.intersectRayBoundingBox(ray)  == Intersection::INTERSECTION)
         {
            if( parent.isLeaf() ) 
            {
               return this->checkIntersectionWithTriFaces(ray, *parent.getTriFaces(), parent.getNodes(), mailbox);
            } 
            else
            {
               int sum = 0;
               if( child1 )
               {
                  int erg = child1->intersectRay(ray, *this, mailbox);
                  if(erg < 0)
                  {
                     return erg;
                  }
                  sum += erg;
               }
               if( child2 ) 
               {
                  int erg = child2->intersectRay(ray, *this, mailbox);
                  if(erg < 0)
                  {
                     return erg;
                  }
                  sum += erg;
               }
               return sum;
            }
         } 
         else 
         {
            return 0;
         }
      }
      /* ======================================================================================= */

   private:
      int checkIntersectionWithTriFaces(const Ray<T>& ray, std::vector<GbTriFaceMesh3D::TriFace>& triFaces, std::vector<GbTriFaceMesh3D::Vertex>& nodes, std::set< UbKeys::Key3<int> >& mailbox) const
      {
         T e1x,e1y,e1z,e2x,e2y,e2z,px,py,pz,a,f,sx,sy,sz,u,qx,qy,qz,v,factor;

         int counter = 0, iSec = 0;

         for( std::size_t i=0; i<triFaces.size(); i++ )
         {
            GbTriFaceMesh3D::TriFace& triFace = triFaces[i];

            if( mailbox.find( UbKeys::Key3<int>(triFace.getIndexVertex1(), triFace.getIndexVertex2(), triFace.getIndexVertex3() ) )==mailbox.end() ) 
            {
               mailbox.insert( UbKeys::Key3<int>(triFace.getIndexVertex1(), triFace.getIndexVertex2(), triFace.getIndexVertex3() ) ); //schon hier rein, ansonsten muss man es unten bei JEDEm continue und am ende des ifs machen

               GbTriFaceMesh3D::Vertex& v1 = triFace.getNode(0, nodes);
               GbTriFaceMesh3D::Vertex& v2 = triFace.getNode(1, nodes);
               GbTriFaceMesh3D::Vertex& v3 = triFace.getNode(2, nodes);

               //////////////////////////////////////////////////////////////////////////
               //Raytracing - start(  Anm.: prüft NUR in Strahlrichtung
               // Grundidee: Schnittpunkt in Baryzentrischen Koordinaten besimmten
               // t(u,v,w) = w*v0 + u*v1 + v*v2 
               // mit w = 1.0-u-v, da fuer alle Punkte (u,v,w) im Dreick gilt u+v+w = 1
               // wenn u, v oder w == 0 -> Punkt liegt auf Kante
               // wenn u, v oder w == 1 -> Punkt liegt auf Eckpunkt (-> die anderen Werte muessen 0 )
               
               //e1 = v1 - v0
               e1x = v2.x-v1.x;
               e1y = v2.y-v1.y;
               e1z = v2.z-v1.z;

               //e2 = v2 - v0
               e2x = v3.x-v1.x;
               e2y = v3.y-v1.y;
               e2z = v3.z-v1.z;

               //p = d x e2
               px = ray.directionY*e2z - ray.directionZ*e2y;
               py = ray.directionZ*e2x - ray.directionX*e2z;
               pz = ray.directionX*e2y - ray.directionY*e2x;

               //a = e1 dot p
               a = e1x*px + e1y*py + e1z*pz;
               //if( fabs(a)<1.E-10 ) continue;
               if( fabs(a) < UbMath::Epsilon<T>::val() ) continue;
               f = T(1.0/a);

               //s = o - v0
               sx = ray.originX - v1.x;
               sy = ray.originY - v1.y;
               sz = ray.originZ - v1.z;

               //u = f * ( s dot p)
               u = f * ( sx*px + sy*py + sz*pz );
               
               //u ist nur gueltig in [0;1] 
               if( ( UbMath::less(   u, T(0.0) ) ) || ( UbMath::greater(u, T(1.0) ) ) ) 
               {
                  continue;
               }

               //q = s x e1
               qx = sy*e1z - sz*e1y;
               qy = sz*e1x - sx*e1z;
               qz = sx*e1y - sy*e1x;

               //v = f*(e2 dot q)
               v = f * (ray.directionX*qx + ray.directionY*qy + ray.directionZ*qz);
               
               //v ist nur gueltig in [0;1] da aber v bereits gueltig ist -> u+v darf nicht > 1.0 sein ;-)
               if(   ( UbMath::less(v, T(0.0) ) ) || ( UbMath::greater(u+v, T(1.0) ) ) ) 
               {
                  continue;
               }

               //t = f * (e2 dot q)
               factor = f * (e2x*qx + e2y*qy + e2z*qz);
               //Raytracing - end
               //////////////////////////////////////////////////////////////////////////

               if( UbMath::zero( factor ) ) 
               {
                  return Intersection::ON_BOUNDARY; //ray.Org liegt direkt auf einem dreieck --> boundary
               }
               if( factor < 0.0 )
               {
                  continue;  //Schnittpunkt liegt in entgegengesetzter Strahlrichtung
               }

               //edge tests
               //wenn u, v oder w ==0 -> Punkt liegt auf Kante bzw. Eckpunkt
               if( UbMath::zero(u)          )  return Intersection::INTERSECT_EDGE;
               if( UbMath::zero(v)          )  return Intersection::INTERSECT_EDGE;
               if( UbMath::zero(T(1.0)-u-v) )  return Intersection::INTERSECT_EDGE;

               counter++;
            }
         }
         return counter;
      }
   };
}

#endif //KDCOUNTRAYLINEINTERSECTIONHANDLER_H

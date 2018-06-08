#ifndef KDCOUNTLINEINTERSECTIONHANDLER_H
#define KDCOUNTLINEINTERSECTIONHANDLER_H

#include <basics/utilities/UbTuple.h>
#include <basics/utilities/UbKeys.h>

#include <numerics/geometry3d/GbTriFaceMesh3D.h>

#include <numerics/geometry3d/KdTree/KdNode.h>
#include <numerics/geometry3d/KdTree/KdUtilities.h>
#include <numerics/geometry3d/KdTree/intersectionhandler/KdLineIntersectionHandler.h>

#include <set>

namespace Kd
{
   template< typename T >
   class CountLineIntersectionHandler : public LineIntersectionHandler<T> 
   {
   public:
      bool intersectLine(const UbTuple<T,T,T>& n1, const UbTuple<T,T,T>& n2, Node<T>& parent, Node<T>*& child1, Node<T>*& child2) const
      {
         if( parent.intersectLineBoundingBox(n1, n2)  == Intersection::INTERSECTION)
         {
            if( parent.isLeaf() ) 
            {
               std::vector<GbTriFaceMesh3D::TriFace>& triFaces = *parent.getTriFaces();
               std::vector<GbTriFaceMesh3D::Vertex>& nodes = parent.getNodes();

         for( std::size_t i=0; i<triFaces.size(); i++ )
         {
            GbTriFaceMesh3D::TriFace& triFace = triFaces[i];

                  if( Kd::intersectLine(n1, n2, triFace, nodes) ) return true;
               }
               return false;
               }
            else
               {
               if( child1 )
               {
                  if (child1->intersectLine(n1, n2, *this)) return true;
               }
               if( child2 ) 
               {
                  if (child2->intersectLine(n1, n2, *this)) return true;
               }
            }
         }
         return false;
      }
      /* ======================================================================================= */
   };
}

#endif //KDCOUNTLINEINTERSECTIONHANDLER_H

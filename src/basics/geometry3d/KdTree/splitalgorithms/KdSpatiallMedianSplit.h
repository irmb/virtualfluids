#ifndef SPATIALLMEDIANSPLIT_H
#define SPATIALLMEDIANSPLIT_H

#include <basics/utilities/UbMath.h>
#include <geometry3d/GbTriFaceMesh3D.h>

#include <geometry3d/KdTree/splitalgorithms/KdSplitAlgorithm.h>

namespace Kd
{
   template< typename T >
   class SpatialMedianSplit : public SplitAlgorithm<T> 
   {
      /* ======================================================================================= */
      SplitCandidate<T> findBestSplitCandidate(const int& level, const int& maxLevel, Node<T>& node )  const
      {
         if(   node.getTriFaces()->size() <= 24 //max triangles in node 
            || level >= maxLevel        )
         {
            return SplitCandidate<T>();
         }

         T dx = std::fabs(node.x[1] - node.x[0]);
         T dy = std::fabs(node.y[1] - node.y[0]);
         T dz = std::fabs(node.z[1] - node.z[0]);

         if     ( UbMath::equal(dx, UbMath::max(dx, dy, dz) ) ) return SplitCandidate<T>(Axis::X, node.x[0] + 0.5 * dx, 0, 0, 0);
         else if( UbMath::equal(dy, UbMath::max(dy, dz    ) ) ) return SplitCandidate<T>(Axis::Y, node.y[0] + 0.5 * dy, 0, 0, 0);

         return SplitCandidate<T>(Axis::Z, node.z[0] + 0.5 * dz, 0, 0, 0);

      }
      /* ======================================================================================= */
      void distributeTriFaces(const SplitCandidate<T>& candidate, std::vector<GbTriFaceMesh3D::TriFace>& primitives_child1, std::vector<GbTriFaceMesh3D::TriFace>& primitives_child2, Node<T>& node) const
      {
         if( !node.getTriFaces() )  throw UbException(UB_EXARGS, "null pointer");

         std::vector<GbTriFaceMesh3D::TriFace>& srcTriFaces = *node.getTriFaces();
         std::vector<GbTriFaceMesh3D::Vertex>&  srcNodes    =  node.getNodes();
         std::vector<T> projection;

         for(std::size_t i=0; i<srcTriFaces.size(); i++) 
         {
            GbTriFaceMesh3D::TriFace& triFace = srcTriFaces[i];
            Kd::project2Axis(triFace, srcNodes, candidate.axis, projection);

            T& min = projection[0];
            T& max = projection[2];
         
            // case 1 : object inside plane
            if( UbMath::equal(min,max) )
            {
               if( UbMath::equal(min,candidate.position) ) 
               {
                  primitives_child1.push_back(triFace);
                  primitives_child2.push_back(triFace);
               }
               else if( UbMath::less(min, candidate.position) )
               {
                  primitives_child1.push_back(triFace);
               } 
               else if( UbMath::greater(min, candidate.position) )
               {
                  primitives_child2.push_back(triFace);
               }
            } 
            // case 2 : object on left side of plane
            else if( UbMath::lessEqual(max, candidate.position) )
            {
               primitives_child1.push_back(triFace);
            } 
            // case 3 : object on right side of plane
            else if( UbMath::greaterEqual(min, candidate.position) )
            {
               primitives_child2.push_back(triFace);
            }
            // case 4 : object in both nodes
            else 
            {
               primitives_child1.push_back(triFace);
               primitives_child2.push_back(triFace);
            }
         }

         node.deleteTriFaces();
      }
   };
}

#endif //SPATIALLMEDIANSPLIT_H

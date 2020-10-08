#ifndef KDSPLITALGORITHM_H
#define KDSPLITALGORITHM_H

#include <geometry3d/GbTriFaceMesh3D.h>
//#include <geometry3d/KdTree/Node.h>
#include <geometry3d/KdTree/KdSplitCandidate.h>

#include <vector>

namespace Kd
{
   template< typename T >
   class Node;

   template< typename T >
   class SplitAlgorithm 
   {
   public:
      virtual SplitCandidate<T> findBestSplitCandidate(const int& level, const int& maxLevel, Node<T>& node ) const = 0;
      virtual void distributeTriFaces(const SplitCandidate<T>& candidate, std::vector<GbTriFaceMesh3D::TriFace>& triFacesForChild1, std::vector<GbTriFaceMesh3D::TriFace>& triFacesForChild2, Node<T>& node) const=0;
      virtual ~SplitAlgorithm() = default;
   };
}


#endif  //KDSPLITALGORITHM_H

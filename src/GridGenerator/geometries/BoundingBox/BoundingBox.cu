#include "BoundingBox.cuh"

#include "../Triangle/Triangle.cuh"
#include "../Vertex/Vertex.cuh"
#include <GridGenerator/utilities/math/CudaMath.cuh>

#include <limits>

#include "Serialization/BoundingBoxMemento.h"


template <typename T>
 BoundingBox<T>::BoundingBox(T minX, T maxX, T minY, T maxY, T minZ, T maxZ) : minX(minX), maxX(maxX), minY(minY), maxY(maxY), minZ(minZ), maxZ(maxZ) {}

 template <typename T>
 HOSTDEVICE BoundingBox<T>::BoundingBox() :
     minX(0),
     maxX(0),
     minY(0),
     maxY(0),
     minZ(0),
     maxZ(0) {}

 BoundingBox<real>::BoundingBox(const BoundingBox<real> &t) : minX(t.minX), maxX(t.maxX), minY(t.minY), maxY(t.maxY), minZ(t.minZ), maxZ(t.maxZ) {}

 BoundingBox<int>::BoundingBox(const BoundingBox<real> &box)
 {
     this->minX = (int)(box.minX);
     this->minY = (int)(box.minY);
     this->minZ = (int)(box.minZ);

     this->maxX = (int)floor(box.maxX + 1);
     this->maxY = (int)floor(box.maxY + 1);
     this->maxZ = (int)floor(box.maxZ + 1);
 }

 template <typename T>
 BoundingBox<T>::BoundingBox(const BoundingBox<int> &box)
 {
	 this->minX = (real)box.minX;
	 this->minY = (real)box.minY;
	 this->minZ = (real)box.minZ;

	 this->maxX = (real)box.maxX;
	 this->maxY = (real)box.maxY;
	 this->maxZ = (real)box.maxZ;
 }


 HOSTDEVICE BoundingBox<real> BoundingBox<real>::makeExactBox(const Triangle &t)
 {
	 BoundingBox<real> box;
     t.setMinMax(box.minX, box.maxX, box.minY, box.maxY, box.minZ, box.maxZ);
	 return box;
 }


 HOSTDEVICE BoundingBox<real> BoundingBox<real>::makeRealNodeBox(const Triangle& t, const real& delta)
{
    BoundingBox<real> box;

    real minX, maxX, minY, maxY, minZ, maxZ;
    t.setMinMax(minX, maxX, minY, maxY, minZ, maxZ);

    calculateMinMaxOnNodes(box.minX, box.maxX, minX, maxX, delta);
    calculateMinMaxOnNodes(box.minY, box.maxY, minY, maxY, delta);
    calculateMinMaxOnNodes(box.minZ, box.maxZ, minZ, maxZ, delta);

    return box;
}


template <typename T>
 HOST BoundingBox<T> BoundingBox<T>::makeInvalidMinMaxBox()
 {
     BoundingBox<T> box = BoundingBox<T>(std::numeric_limits<T>::max(),
         std::numeric_limits<T>::lowest(),
         std::numeric_limits<T>::max(),
         std::numeric_limits<T>::lowest(),
         std::numeric_limits<T>::max(),
         std::numeric_limits<T>::lowest());
     return box;
 }

 HOSTDEVICE BoundingBox<int> BoundingBox<int>::makeNodeBox(const Triangle &t)
 {
	 BoundingBox<int> box;

     real minX, maxX, minY, maxY, minZ, maxZ;
     t.setMinMax(minX, maxX, minY, maxY, minZ, maxZ);

	 calculateMinMaxOnNodes(box.minX, box.maxX, minX, maxX);
	 calculateMinMaxOnNodes(box.minY, box.maxY, minY, maxY);
	 calculateMinMaxOnNodes(box.minZ, box.maxZ, minZ, maxZ);
	 return box;
 }

 template <>
 HOSTDEVICE void BoundingBox<real>::calculateMinMaxOnNodes(real &minNode, real &maxNode, const real &minExact, const real &maxExact, const real& delta)
 {
     minNode = ceil(minExact - 1.0);
     maxNode = floor(maxExact + 1.0);
     if (minNode + 0.5 < minExact)
         minNode += 0.5;


     if (maxNode - 0.5 > maxExact)
         maxNode -= 0.5;
 }

 template <>
 HOSTDEVICE void BoundingBox<int>::calculateMinMaxOnNodes(int &minNode, int &maxNode, const real &minExact, const real &maxExact)
 {
     minNode = ceil(minExact - 1.0);
     maxNode = floor(maxExact + 1.0);
 }

 void BoundingBox<real>::setMinMax(const Triangle& t)
 {
     t.setMinMax(minX, maxX, minY, maxY, minZ, maxZ);
 }

 template <typename T>
 bool BoundingBox<T>::intersect(const Triangle &t) const
 {
	 if (isInside(t.v1) || isInside(t.v2) || isInside(t.v3))
		 return true;
	 return false;
 }

 template <typename T>
 bool BoundingBox<T>::isInside(const Triangle &t) const
 {
	 if (isInside(t.v1) && isInside(t.v2) && isInside(t.v3))
		 return true;
	 return false;
 }

 template <typename T>
 bool BoundingBox<T>::isInside(const Vertex &v) const
 {
     if (v.isXbetween(minX, maxX) && v.isYbetween(minY, maxY) && v.isZbetween(minZ, maxZ))
		 return true;
	 return false;
 }

 template <typename T>
 std::vector<std::vector<Vertex> > BoundingBox<T>::getIntersectionPoints(const BoundingBox<real> &b) const
 {
	 std::vector<std::vector<Vertex> > intersectionBox;
	 intersectionBox.resize(6);

	 int intersects = 0;
	 if (b.minX < maxX && b.maxX > maxX) { //maxX is intersect
		 intersectionBox[intersects].push_back(Vertex((real)maxX, (real)minY, (real)minZ));
		 intersectionBox[intersects].push_back(Vertex((real)maxX, (real)maxY, (real)minZ));
		 intersectionBox[intersects].push_back(Vertex((real)maxX, (real)minY, (real)maxZ));
		 intersects++;
	 }
	 if (b.minX < minX && b.maxX > minX) { //minX is intersect
		 intersectionBox[intersects].push_back(Vertex((real)minX, (real)minY, (real)minZ));
		 intersectionBox[intersects].push_back(Vertex((real)minX, (real)maxY, (real)minZ));
		 intersectionBox[intersects].push_back(Vertex((real)minX, (real)minY, (real)maxZ));
		 intersects++;
	 }
	 if (b.minY < minY && b.maxY > minY) { //minY is intersect
		 intersectionBox[intersects].push_back(Vertex((real)minX, (real)minY, (real)minZ));
		 intersectionBox[intersects].push_back(Vertex((real)maxX, (real)minY, (real)minZ));
		 intersectionBox[intersects].push_back(Vertex((real)minX, (real)minY, (real)maxZ));
		 intersects++;
	 }
	 if (b.minY < maxY && b.maxY > maxY) { //maxY is intersect
		 intersectionBox[intersects].push_back(Vertex((real)minX, (real)maxY, (real)minZ));
		 intersectionBox[intersects].push_back(Vertex((real)maxX, (real)maxY, (real)minZ));
		 intersectionBox[intersects].push_back(Vertex((real)minX, (real)maxY, (real)maxZ));
		 intersects++;
	 }
	 if (b.minZ < minZ && b.maxZ > minZ) { //minZ is intersect
		 intersectionBox[intersects].push_back(Vertex((real)minX, (real)minY, (real)minZ));
		 intersectionBox[intersects].push_back(Vertex((real)maxX, (real)minY, (real)minZ));
		 intersectionBox[intersects].push_back(Vertex((real)minX, (real)maxY, (real)minZ));
		 intersects++;
	 }
	 if (b.minZ < maxZ && b.maxZ > maxZ) { //maxZ is intersect
		 intersectionBox[intersects].push_back(Vertex((real)minX, (real)minY, (real)maxZ));
		 intersectionBox[intersects].push_back(Vertex((real)maxX, (real)minY, (real)maxZ));
		 intersectionBox[intersects].push_back(Vertex((real)minX, (real)maxY, (real)maxZ));
		 intersects++;
	 }

	 return intersectionBox;
 }


 template <typename T>
 bool BoundingBox<T>::intersect(const BoundingBox<T> &box) const
 {
	 struct Vertex v[8];
	 box.getPoints(v);

	 for (int i = 0; i < 8; i++) {
		 if (isInside(v[i]))
			 return true;
	 }
	 return false;
 }

 template <typename T>
 void BoundingBox<T>::getPoints(Vertex v[8]) const
 {
	 v[0] = Vertex(minX, minY, minZ);
	 v[1] = Vertex(maxX, minY, minZ);
	 v[2] = Vertex(minX, maxY, minZ);
	 v[3] = Vertex(maxX, maxY, minZ);

	 v[4] = Vertex(minX, minY, maxZ);
	 v[5] = Vertex(maxX, minY, maxZ);
	 v[6] = Vertex(minX, maxY, maxZ);
	 v[7] = Vertex(maxX, maxY, maxZ);
 }


 template <typename T>
 void BoundingBox<T>::print() const
 {
	 printf("min/max - x: %2.4f/ %2.4f, y: %2.4f, %2.4f, z: %2.4f, %2.4f \n", (real)minX, (real)maxX, (real)minY, (real)maxY, (real)minZ, (real)maxZ);
 }


 HOST bool BoundingBox<real>::operator==(const BoundingBox<real> &box) const
 {
     return CudaMath::equal(minX, box.minX)
         && CudaMath::equal(maxX, box.maxX)
         && CudaMath::equal(minY, box.minY)
         && CudaMath::equal(maxY, box.maxY)
         && CudaMath::equal(minZ, box.minZ)
         && CudaMath::equal(maxZ, box.maxZ);
 }


 template <typename T>
 HOST BoundingBoxMemento BoundingBox<T>::getState() const
 {
     BoundingBoxMemento memento;
     memento.minX = (real)minX;
     memento.maxX = (real)maxX;
     memento.minY = (real)minY;
     memento.maxY = (real)maxY;
     memento.minZ = (real)minZ;
     memento.maxZ = (real)maxZ;
     return memento;
 }

 template <typename T>
 HOST void BoundingBox<T>::setState(const BoundingBoxMemento &memento)
 {
     this->minX = (T)memento.minX;
     this->maxX = (T)memento.maxX;
     this->minY = (T)memento.minY;
     this->maxY = (T)memento.maxY;
     this->minZ = (T)memento.minZ;
     this->maxZ = (T)memento.maxZ;
 }


 template class VF_PUBLIC BoundingBox<int>;
 template class VF_PUBLIC BoundingBox<real>;

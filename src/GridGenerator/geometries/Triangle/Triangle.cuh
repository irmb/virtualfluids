#ifndef Triangle_h
#define Triangle_h

#include "GridGenerator/global.h"
#include "GridGenerator_EXPORT.h"

#include "../Vertex/Vertex.cuh"

#include <memory>

class TriangleMemento;


struct GridGenerator_EXPORT Triangle
{
    Vertex v1, v2, v3, normal;
    doubflo alphaAngles[3];
    
	HOSTDEVICE Triangle(Vertex &v1, Vertex &v2, Vertex &v3, Vertex &normal);
	HOSTDEVICE Triangle(Vertex &v1, Vertex &v2, Vertex &v3);
	HOSTDEVICE Triangle();

	HOSTDEVICE void set(const Vertex &v1, const Vertex &v2, const Vertex &v3);
	HOSTDEVICE void calcNormal();


	Vertex getCenterOfMass() const;
	doubflo getHalfAngleBetweenToAdjacentTriangle(const Triangle &t2) const;
	int isEqual(const Triangle &t2) const;
	bool doesNormalsShowToEachOther(const  Triangle &t2) const;
	int getCommonEdge(const Triangle &t2) const;

	bool contains(const Vertex& v)const;
	HOSTDEVICE int getNumberOfCommonEdge(const Triangle &t2) const;
	HOSTDEVICE int getTriangleIntersection(const Vertex &P, const Vertex &direction, Vertex &pointOnTri, doubflo &qVal) const;
	HOSTDEVICE void print() const;

    HOSTDEVICE char isUnderFace(const Vertex &point) const;

    HOSTDEVICE bool isUnterExtendedFace(const Vertex & point, doubflo &s, doubflo &delta) const;
    HOSTDEVICE bool isNotNextToFace(const Vertex &point) const;
    HOSTDEVICE bool isUnderAngleToNeighbors(const Vertex &point) const;
    HOSTDEVICE void getClosestPointsOnEdges(Vertex arr[], const Vertex &P) const;
    HOSTDEVICE doubflo getPerpedicularDistanceFrom(const Vertex &P) const;
    HOSTDEVICE Vertex getPerpedicularPointFrom(const Vertex &P) const;
    HOSTDEVICE bool isQNode(const Vertex & point, const doubflo &s, const doubflo &delta) const;
    HOSTDEVICE bool isStopper(const Vertex & point) const;

    HOST bool operator==(const Triangle &t) const;

    HOST TriangleMemento getState() const;
    HOST void setState(const TriangleMemento &memento);


    HOSTDEVICE void Triangle::setMinMax(doubflo &minX, doubflo &maxX, doubflo &minY, doubflo &maxY, doubflo &minZ, doubflo &maxZ) const;
};

#endif

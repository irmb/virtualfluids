#ifndef Triangle_h
#define Triangle_h

#include "GridGenerator/global.h"

#include "../Vertex/Vertex.h"

#include <memory>

class TriangleMemento;


struct VF_PUBLIC Triangle
{
    Vertex v1, v2, v3, normal;
    real alphaAngles[3];
    real layerThickness;
    
	HOSTDEVICE Triangle(Vertex &v1, Vertex &v2, Vertex &v3, Vertex &normal);
	HOSTDEVICE Triangle(Vertex &v1, Vertex &v2, Vertex &v3);
	HOSTDEVICE Triangle();

	HOSTDEVICE void set(const Vertex &v1, const Vertex &v2, const Vertex &v3);
    HOSTDEVICE void set(int index, Vertex value);
    HOSTDEVICE Vertex get(int index);
	HOSTDEVICE void calcNormal();

    HOSTDEVICE void initalLayerThickness(real delta);


	Vertex getCenterOfMass() const;
	real getHalfAngleBetweenToAdjacentTriangle(const Triangle &t2) const;
	int isEqual(const Triangle &t2) const;
	bool doesNormalsShowToEachOther(const  Triangle &t2) const;
	int getCommonEdge(const Triangle &t2) const;

	bool contains(const Vertex& v)const;
	HOSTDEVICE int getNumberOfCommonEdge(const Triangle &t2) const;
	HOSTDEVICE int getTriangleIntersection(const Vertex &P, const Vertex &direction, Vertex &pointOnTri, real &qVal) const;
	HOSTDEVICE void print() const;

    HOSTDEVICE char isUnderFace(const Vertex &point) const;

    HOSTDEVICE bool isUnterExtendedFace(const Vertex & point, real &s) const;
    HOSTDEVICE bool isNotNextToFace(const Vertex &point) const;
    HOSTDEVICE bool isUnderAngleToNeighbors(const Vertex &point) const;
    HOSTDEVICE void getClosestPointsOnEdges(Vertex arr[], const Vertex &P) const;
    HOSTDEVICE real getPerpedicularDistanceFrom(const Vertex &P) const;
    HOSTDEVICE Vertex getPerpedicularPointFrom(const Vertex &P) const;
    HOSTDEVICE bool isQNode(const Vertex & point, const real &s) const;
    HOSTDEVICE bool isNegativeDirectionBorder(const Vertex & point) const;

    HOST bool operator==(const Triangle &t) const;

    HOST TriangleMemento getState() const;
    HOST void setState(const TriangleMemento &memento);


    HOSTDEVICE void setMinMax(real &minX, real &maxX, real &minY, real &maxY, real &minZ, real &maxZ) const;
};

#endif

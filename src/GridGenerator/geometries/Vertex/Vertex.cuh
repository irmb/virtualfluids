#ifndef VERTEX_H
#define VERTEX_H

#include "GridGenerator/global.h"
#include "GridGenerator_EXPORT.h"

#include <stdio.h>
#include <cuda_runtime.h>
#include <memory>
#include <ostream>

class VertexMemento;

struct GridGenerator_EXPORT Vertex 
{
public:
    doubflo x, y, z;

	HOSTDEVICE Vertex(doubflo x, doubflo y, doubflo z);
	HOSTDEVICE Vertex();
	HOSTDEVICE Vertex(const Vertex& v);
	HOSTDEVICE ~Vertex() {}

	HOSTDEVICE doubflo getEuclideanDistanceTo(const Vertex &w) const;
	HOSTDEVICE Vertex operator-(const Vertex &v) const;
	HOSTDEVICE Vertex operator+(const Vertex &v) const;
	HOSTDEVICE Vertex operator*(const doubflo value) const;
	HOSTDEVICE doubflo operator*(const Vertex &w) const;
	HOSTDEVICE struct Vertex crossProduct(const Vertex &w) const;
	HOSTDEVICE doubflo length() const;
	HOSTDEVICE void normalize();
	HOSTDEVICE doubflo getMagnitude() const;
	HOSTDEVICE int isEqual(const Vertex &w) const;
	HOSTDEVICE doubflo getInnerAngle(const Vertex &w) const;

	HOSTDEVICE bool operator==(const Vertex &v) const;

    HOST VertexMemento getState() const;
    HOST void setState(const VertexMemento &memento);

    HOST bool isXbetween(doubflo min, doubflo max) const;
    HOST bool isYbetween(doubflo min, doubflo max) const;
    HOST bool isZbetween(doubflo min, doubflo max) const;

    HOSTDEVICE static void setMinMax(doubflo &minX, doubflo &maxX, doubflo &minY, doubflo &maxY, doubflo &minZ, doubflo &maxZ, const Vertex &v1, const Vertex &v2, const Vertex &v3); 
    HOSTDEVICE static void calculateMinMax(const doubflo &value1, const doubflo &value2, const doubflo &value3, doubflo &min, doubflo &max);

    HOSTDEVICE void print() const;
    HOST void print(std::ostream &ost) const;
    HOST void printFormatted(std::ostream &ost) const;

};



#endif

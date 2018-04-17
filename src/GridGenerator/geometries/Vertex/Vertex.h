#ifndef VERTEX_H
#define VERTEX_H

#include "GridGenerator/global.h"


#include <stdio.h>
#include <cuda_runtime.h>
#include <memory>
#include <ostream>

class VertexMemento;

struct VF_PUBLIC Vertex 
{
public:
    real x, y, z;

	HOSTDEVICE Vertex(real x, real y, real z);
	HOSTDEVICE Vertex();
	HOSTDEVICE Vertex(const Vertex& v);
	HOSTDEVICE ~Vertex() {}

	HOSTDEVICE real getEuclideanDistanceTo(const Vertex &w) const;
	HOSTDEVICE Vertex operator-(const Vertex &v) const;
	HOSTDEVICE Vertex operator+(const Vertex &v) const;
	HOSTDEVICE Vertex operator*(const real value) const;
	HOSTDEVICE real operator*(const Vertex &w) const;
	HOSTDEVICE struct Vertex crossProduct(const Vertex &w) const;
	HOSTDEVICE real length() const;
	HOSTDEVICE void normalize();
	HOSTDEVICE real getMagnitude() const;
	HOSTDEVICE int isEqual(const Vertex &w) const;
	HOSTDEVICE real getInnerAngle(const Vertex &w) const;

	HOSTDEVICE bool operator==(const Vertex &v) const;

    HOST VertexMemento getState() const;
    HOST void setState(const VertexMemento &memento);

    HOST bool isXbetween(real min, real max) const;
    HOST bool isYbetween(real min, real max) const;
    HOST bool isZbetween(real min, real max) const;

    HOSTDEVICE static void setMinMax(real &minX, real &maxX, real &minY, real &maxY, real &minZ, real &maxZ, const Vertex &v1, const Vertex &v2, const Vertex &v3); 
    HOSTDEVICE static void calculateMinMax(const real &value1, const real &value2, const real &value3, real &min, real &max);

    HOSTDEVICE void print() const;
    HOST void print(std::ostream &ost) const;
    HOST void printFormatted(std::ostream &ost) const;

};



#endif

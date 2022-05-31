#ifndef VERTEX_H
#define VERTEX_H

#include <stdio.h>
#include <cuda_runtime.h>
#include <memory>
#include <ostream>

#include "gpu/GridGenerator/global.h"

class VertexMemento;

struct GRIDGENERATOR_EXPORT Vertex
{
public:
    real x, y, z;

	HOSTDEVICE Vertex(real x, real y, real z);
	HOSTDEVICE Vertex();
	HOSTDEVICE Vertex(const Vertex& v);
    HOSTDEVICE Vertex& operator=(const Vertex&);
	HOSTDEVICE ~Vertex() {}

	HOSTDEVICE real getEuclideanDistanceTo(const Vertex &w) const;
	HOSTDEVICE Vertex operator-(const Vertex &v) const;
	HOSTDEVICE Vertex operator+(const Vertex &v) const;
	HOSTDEVICE Vertex operator*(const real& value) const;
    HOSTDEVICE Vertex operator/(const real& value) const;

	HOSTDEVICE real operator*(const Vertex &w) const;
	HOSTDEVICE struct Vertex crossProduct(const Vertex &w) const;
	HOSTDEVICE real length() const;
	HOSTDEVICE void normalize();
	HOSTDEVICE real getMagnitude() const;
	HOSTDEVICE int isEqual(const Vertex &w) const;
	HOSTDEVICE real getInnerAngle(const Vertex &w) const;

	HOSTDEVICE bool operator==(const Vertex &v) const;

    CUDA_HOST VertexMemento getState() const;
    CUDA_HOST void setState(const VertexMemento &memento);

    CUDA_HOST bool isXbetween(real min, real max) const;
    CUDA_HOST bool isYbetween(real min, real max) const;
    CUDA_HOST bool isZbetween(real min, real max) const;

    HOSTDEVICE static void setMinMax(real &minX, real &maxX, real &minY, real &maxY, real &minZ, real &maxZ, const Vertex &v1, const Vertex &v2, const Vertex &v3); 
    HOSTDEVICE static void calculateMinMax(const real &value1, const real &value2, const real &value3, real &min, real &max);

    HOSTDEVICE void print() const;
    CUDA_HOST void print(std::ostream &ost) const;
    CUDA_HOST void printFormatted(std::ostream &ost) const;

};



#endif

#ifndef TriangleNeighborFinder_h
#define TriangleNeighborFinder_h

#include "GridGenerator/global.h"
#include "GridGenerator_EXPORT.h"

#include <vector>

struct Triangle;
struct Vertex;
struct Geometry;

struct IntegerPtr2D {
    int *ptr;
    int size;
    int DIM;
};

class TriangleNeighborFinder
{
public:
    GridGenerator_EXPORT TriangleNeighborFinder(Triangle *triangles_ptr, int size);
    GridGenerator_EXPORT ~TriangleNeighborFinder();
    
    void GridGenerator_EXPORT fillWithNeighborIndices(IntegerPtr2D *indices, Triangle *triangles);
	void GridGenerator_EXPORT fillWithNeighborAngles(Geometry *geom) const;

private:
    void initalSortedInSpaceWithCoords(Triangle *triangles_ptr, int size);
    void fillSortedInSpaceWithFirstVertexAndCoordinateIDs(int numberOfRows);
    void fillVectorWithIndicesOfTriangleNeighbors();
    unsigned int findTriangleID(unsigned int uniqueCoordID);
    Vertex getCoordinatesIDfromTriangle(int triangleID);
    bool isNeighborNotAlreadyInside(unsigned int iTriangle, unsigned int jTriangle);
    std::vector<std::vector<unsigned int> > indicesOfTriangleNeighbors;

    bool isTriangleNeighborOfParentTriangle(Vertex, Vertex);

    int numberOfRows;
    doubflo **sortedToTriangles;
    doubflo **sortedInSpace;
};

#endif

#ifndef TriangleNeighborFinder_h
#define TriangleNeighborFinder_h

#include "GridGenerator/global.h"

#include <vector>

struct IDS {
    enum IDs { vertexID = 0, firstVertexID = 1, coordinateID = 2, uniqueCoordID = 3, x = 3, y = 4, z = 5 };
};

struct Triangle;
struct Vertex;
class TriangularMesh;

struct IntegerPtr2D {
    int *ptr;
    int size;
    int DIM;
};

class TriangleNeighborFinder
{
public:
    VIRTUALFLUIDS_GPU_EXPORT TriangleNeighborFinder(Triangle *triangles, int size);
    VIRTUALFLUIDS_GPU_EXPORT ~TriangleNeighborFinder();
    
    std::vector<int> getTriangleIDsWithCommonVertex(int vertexID) const;
    std::vector< std::vector<Triangle> > getTrianglesPerVertex() const;

    void VIRTUALFLUIDS_GPU_EXPORT fillWithNeighborIndices(IntegerPtr2D *indices, Triangle *triangles);
	void VIRTUALFLUIDS_GPU_EXPORT fillWithNeighborAngles(TriangularMesh *geom) const;

    void printSortedToTriangles() const;
    void printSortedInSpace() const;

private:
    void initalSortedInSpaceWithCoords(Triangle *triangles_ptr, int size);
    void fillSortedInSpaceWithFirstVertexAndCoordinateIDs(int numberOfRows);
    void fillVectorWithIndicesOfTriangleNeighbors();
    unsigned int findTriangleID(unsigned int uniqueCoordID);
    Vertex getCoordinatesIDfromTriangle(int triangleID);
    bool isNeighborNotAlreadyInside(unsigned int iTriangle, unsigned int jTriangle);

    bool isTriangleNeighborOfParentTriangle(Vertex, Vertex);

    int numberOfRows;
    Triangle *triangles;

public:
    std::vector<std::vector<uint> > indicesOfTriangleNeighbors;
    real **sortedToTriangles;
    real **sortedInSpace;
};

#endif

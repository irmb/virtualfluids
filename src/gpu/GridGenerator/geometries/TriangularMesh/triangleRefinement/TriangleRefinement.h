#ifndef TriangleRefinement_h
#define TriangleRefinement_h

#include <vector>

#include "global.h"

#include "geometries/TriangularMesh/triangleNeighborFinder/TriangleNeighborFinder.h"

struct Triangle;
struct Vertex;
struct IntegerPtr2D;

class  TriangleRefinement
{
public:
    VIRTUALFLUIDS_GPU_EXPORT TriangleRefinement(std::vector<Triangle> *triangles);
    VIRTUALFLUIDS_GPU_EXPORT ~TriangleRefinement();

    void VIRTUALFLUIDS_GPU_EXPORT refine(int iTriangle);
    static void VIRTUALFLUIDS_GPU_EXPORT refine(Triangle t, Triangle &firstNewTriangle, Triangle &secondNewTriangle);
    
    void VIRTUALFLUIDS_GPU_EXPORT refineUntilMinDistance(double d_min);
    void VIRTUALFLUIDS_GPU_EXPORT refineUntilcountTriangle(int countTri);
    void VIRTUALFLUIDS_GPU_EXPORT redoubleTriangles();

    static VIRTUALFLUIDS_GPU_EXPORT Vertex getHalfVertex(const Vertex &v, const Vertex &w);

private:
    std::vector<Triangle> *triangles;
    IntegerPtr2D indices;

    int findIndexFromTriangleWithLongestEdge(double *d);
    std::vector<Vertex> getVertexArray(int iTriangle);
    Vertex getNewhalfVertexFromTrianglesEdge(std::vector<Vertex> v, int edge);
    void sortNeighborIndices();
    void createTwoTriangles(std::vector<Vertex> v, int edge, Vertex newEdge);
    void eraseOldTrianglesFromVector(int iTriangle, int indexNeighbor);
    int findCommonEdgeFromTriangles(int indexNeighbor, int iTriangle);

public:
    static int VIRTUALFLUIDS_GPU_EXPORT getEdgeWithLongestDistance(Triangle &t);
    static real VIRTUALFLUIDS_GPU_EXPORT getLongestEdgeDistance(Triangle &t);
};


#endif

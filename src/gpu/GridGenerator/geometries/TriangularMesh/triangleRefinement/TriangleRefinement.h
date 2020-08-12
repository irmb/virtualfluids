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
    GRIDGENERATOR_EXPORT TriangleRefinement(std::vector<Triangle> *triangles);
    GRIDGENERATOR_EXPORT ~TriangleRefinement();

    void GRIDGENERATOR_EXPORT refine(int iTriangle);
    static void GRIDGENERATOR_EXPORT refine(Triangle t, Triangle &firstNewTriangle, Triangle &secondNewTriangle);
    
    void GRIDGENERATOR_EXPORT refineUntilMinDistance(double d_min);
    void GRIDGENERATOR_EXPORT refineUntilcountTriangle(int countTri);
    void GRIDGENERATOR_EXPORT redoubleTriangles();

    static GRIDGENERATOR_EXPORT Vertex getHalfVertex(const Vertex &v, const Vertex &w);

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
    static int GRIDGENERATOR_EXPORT getEdgeWithLongestDistance(Triangle &t);
    static real GRIDGENERATOR_EXPORT getLongestEdgeDistance(Triangle &t);
};


#endif

#ifndef TriangleRefinement_h
#define TriangleRefinement_h

#include "GridGenerator_EXPORT.h"
#include <vector>
#include "GridGenerator/global.h"

#include <GridGenerator/utilities/triangleNeighborFinder/TriangleNeighborFinder.h>

struct Triangle;
struct Vertex;
struct IntegerPtr2D;

class  TriangleRefinement
{
public:
    GridGenerator_EXPORT TriangleRefinement(std::vector<Triangle> *triangles);
    GridGenerator_EXPORT ~TriangleRefinement();

    void GridGenerator_EXPORT refine(int iTriangle);
    static void GridGenerator_EXPORT refine(Triangle t, Triangle &firstNewTriangle, Triangle &secondNewTriangle);
    
    void GridGenerator_EXPORT refineUntilMinDistance(double d_min);
    void GridGenerator_EXPORT refineUntilcountTriangle(int countTri);
    void GridGenerator_EXPORT redoubleTriangles();

    static GridGenerator_EXPORT Vertex getHalfVertex(const Vertex &v, const Vertex &w);

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
    static int GridGenerator_EXPORT getEdgeWithLongestDistance(Triangle &t);
    static doubflo GridGenerator_EXPORT getLongestEdgeDistance(Triangle &t);
};


#endif

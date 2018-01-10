#ifndef TriangleRefinement_h
#define TriangleRefinement_h


#include <vector>
#include "GridGenerator/global.h"

#include <GridGenerator/utilities/triangleNeighborFinder/TriangleNeighborFinder.h>

struct Triangle;
struct Vertex;
struct IntegerPtr2D;

class  TriangleRefinement
{
public:
    VF_PUBLIC TriangleRefinement(std::vector<Triangle> *triangles);
    VF_PUBLIC ~TriangleRefinement();

    void VF_PUBLIC refine(int iTriangle);
    static void VF_PUBLIC refine(Triangle t, Triangle &firstNewTriangle, Triangle &secondNewTriangle);
    
    void VF_PUBLIC refineUntilMinDistance(double d_min);
    void VF_PUBLIC refineUntilcountTriangle(int countTri);
    void VF_PUBLIC redoubleTriangles();

    static VF_PUBLIC Vertex getHalfVertex(const Vertex &v, const Vertex &w);

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
    static int VF_PUBLIC getEdgeWithLongestDistance(Triangle &t);
    static real VF_PUBLIC getLongestEdgeDistance(Triangle &t);
};


#endif

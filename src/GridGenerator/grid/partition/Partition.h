#ifndef Partition_H
#define Partition_H

#include "GridGenerator/global.h"


#include <vector>
#include <string>
#include <memory>

template <class T>
class BoundingBox;
struct Triangle;
struct Vertex;
class Grid;
class Transformator;

class VF_PUBLIC Partition
{
public:
    static void partitionGridMesh(const Grid &grid);
    static void partitionGrid(const Grid &grid);

	static std::vector<BoundingBox<int> > getProcessBoxes(const int numberOfProcesses, const int globalNx, const int globalNy, const int globalNz);
    static std::vector<std::vector<int> > partitionBoxes(std::vector<std::vector<BoundingBox<int>> >, int processes, std::vector< std::shared_ptr<Transformator> > transformators);

    static std::vector<std::vector<Triangle> > getTrianglesPerProcess(std::vector<BoundingBox<int>> &boxes, Triangle *triangles, int nTriangles);

    static std::vector<std::vector<Triangle> >  splitTriangles(std::vector<Triangle> &triangleVec, std::vector<BoundingBox<int>> boxes);
    
private:
    Partition(){};
    ~Partition(){};

    static int calcEdgesFromGraph(const Grid &grid);

    static void splitAllBoxesInThreePieces(std::vector<BoundingBox<int>> &boxes, bool splitX, bool splitY, bool splitZ);
    static void splitAllBoxesInTwoPieces(std::vector<BoundingBox<int>> &boxes, bool splitX, bool splitY, bool splitZ);
    static void findMaxBoxSize(std::vector<BoundingBox<int>> &boxes, bool &splitX, bool &splitY, bool &splitZ);

    static void splitAndPushTriangle(BoundingBox<int> &box, std::vector<Triangle> &trianglesPerProcess, std::vector<Triangle> &triangleVec, int indexTriangle);
    static void sliceTriangle(std::vector<Triangle> &out, std::vector<Triangle>& triangleVec, int index, const Vertex& a, const Vertex& b, const Vertex& c, real d1, real d2, real d3);

};

#endif

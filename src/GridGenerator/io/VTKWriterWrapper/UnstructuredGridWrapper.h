//#ifndef UnstructuredGridWrapper_H
//#define UnstructuredGridWrapper_H
//
//
//#include <string>
//#include <vector>
//
//class UnstructuredGridWriter;
//struct Triangle;
//
//template <class T>
//class BoundingBox;
//struct Vertex;
//struct Grid;
//
//class VF_PUBLIC UnstructuredGridWrapper
//{
//public:
//    UnstructuredGridWrapper();
//    ~UnstructuredGridWrapper();
//
//    void addBoundingBox(double boundingBox[6], int type);
//	template <typename T>
//    void addBoundingBox(BoundingBox<T> b, int type);
//    void addVoxel(double node[3], double nodeNX[3], double nodeNY[3], double nodeNZ[3],
//        double nodeNYNX[3], double nodeNZNY[3], double nodeNZNX[3], double nodeNZNYNX[3], int nodeTypes[8]);
//
//    void addTriangle(double v1[3], double v2[3], double v3[3], int type[3]);
//    void addTriangle(Triangle t, int type[3]);
//    void addTriangles(std::vector<Triangle> triangles);
//
//    void addGridPoints(double grid[], int nx, int ny, int nz);
//
//    void addGridPoint(int x, int y, int z, int type);
//    void addGridPoint(Vertex v, int type);
//
//    void writeUnstructuredGridToFile(std::string filename);
//	
//    void displayUnstructuredGrid();
//    void displayUnstructuredGrid(std::string filename);
//
//	template <typename T>
//	static void writeBoxesToFile(std::vector<BoundingBox<T>>, std::string name);
//
//private:
//#ifndef __unix__
//    UnstructuredGridWriter *writer;
//#endif
//
//};
//
//
//#endif

#ifndef STLReader_H
#define STLReader_H

#include "GridGenerator_EXPORT.h"
#include <vector>
#include <fstream>
#include <string>
#include <memory>

class Transformator;
struct Triangle;
struct Vertex;
template<class T>
class BoundingBox;
struct Geometry;

class GridGenerator_EXPORT STLReader
{
public:
	static std::vector<Triangle> readSTL(const std::string& name, const Transformator& trans);
	static std::vector<Triangle> readSTL(const BoundingBox<int> &box, const std::string name, const Transformator& trans);

	static Geometry getGeometry(const std::string& name, const Transformator& trans);
	static Geometry getGeometryFromBinarySTL(const std::string& name, const Transformator& trans);
	static Geometry getGeometryFromASCIISTL(const std::string& name, const Transformator& trans);

    static std::vector<Triangle> readBinarySTL(const std::string& name, const Transformator& trans);
    static std::vector<Triangle> readASCIISTL(const std::string& name, const Transformator& trans);
	static std::vector<Triangle> readBinarySTL(const BoundingBox<int> &box, const std::string name, const Transformator& trans);
	static std::vector<Triangle> readASCIISTL(const BoundingBox<int> &box, const std::string name, const Transformator& trans);

private:
    STLReader(){};
    ~STLReader(){};

    static int countLinesInFile(std::string name);
    static Vertex parseLineToCoordinates(std::ifstream& file, std::string format);
    static Vertex getVertexFromChar(const char* facet);
};

#endif

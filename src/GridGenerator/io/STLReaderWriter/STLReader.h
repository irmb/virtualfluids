#ifndef STLReader_H
#define STLReader_H


#include <vector>
#include <fstream>
#include <string>
#include <memory>

#include <VirtualFluidsDefinitions.h>

class Transformator;
struct Triangle;
struct Vertex;
template<class T>
class BoundingBox;
struct TriangularMesh;

class VF_PUBLIC STLReader
{
public:
    static std::vector<Triangle> readSTL(const std::string& name);
	static std::vector<Triangle> readSTL(const std::string& name, const Transformator& trans);
	static std::vector<Triangle> readSTL(const BoundingBox<int> &box, const std::string name, const Transformator& trans);

	static TriangularMesh* makeTriangularMesh(const std::string& name);
	static TriangularMesh* getGeometryFromBinarySTL(const std::string& name);
	static TriangularMesh* getGeometryFromASCIISTL(const std::string& name);

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

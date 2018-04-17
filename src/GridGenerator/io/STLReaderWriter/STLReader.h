#ifndef STLReader_H
#define STLReader_H


#include <vector>
#include <string>

#include <VirtualFluidsDefinitions.h>

struct Triangle;
struct Vertex;
template<class T>
class BoundingBox;

class VF_PUBLIC STLReader
{
public:
    static std::vector<Triangle> readSTL(const std::string& name);
	static std::vector<Triangle> readSTL(const BoundingBox<int> &box, const std::string& name);

    static std::vector<Triangle> readBinarySTL(const std::string& name);
    static std::vector<Triangle> readASCIISTL(const std::string& name);
	static std::vector<Triangle> readBinarySTL(const BoundingBox<int> &box, const std::string& name);
	static std::vector<Triangle> readASCIISTL(const BoundingBox<int> &box, const std::string& name);

private:
    STLReader(){};
    ~STLReader(){};

    static int countLinesInFile(std::string name);
    static Vertex parseLineToCoordinates(std::ifstream& file, std::string format);
    static Vertex getVertexFromChar(const char* facet);
};

#endif

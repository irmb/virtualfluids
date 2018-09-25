#ifndef STLReader_H
#define STLReader_H


#include <vector>
#include <string>

#include <VirtualFluidsDefinitions.h>

struct Triangle;
struct Vertex;
class BoundingBox;

class VF_PUBLIC STLReader
{
public:

    enum FileType { ascii, binary };

    static std::vector<Triangle> readSTL(const std::string& name);
    static std::vector<Triangle> readSTL(const std::string& name, FileType fileType );
	static std::vector<Triangle> readSTL(const BoundingBox &box, const std::string& name);

    static std::vector<Triangle> readBinarySTL(const std::string& name);
    static std::vector<Triangle> readASCIISTL(const std::string& name);
    static std::vector<Triangle> readASCIISTLWithPatches(const std::string& name);
	static std::vector<Triangle> readBinarySTL(const BoundingBox &box, const std::string& name);
	static std::vector<Triangle> readASCIISTL(const BoundingBox &box, const std::string& name);

private:
    STLReader(){};
    ~STLReader(){};

    static int countLinesInFile(std::string name);
    static Vertex parseLineToCoordinates(std::ifstream& file, std::string format);
    static Vertex parseLineToCoordinates(const std::string& file, const std::string format);
    static Vertex getVertexFromChar(const char* facet);
};

#endif

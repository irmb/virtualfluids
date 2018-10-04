#ifndef STLReader_H
#define STLReader_H


#include <vector>
#include <string>

#include <VirtualFluidsDefinitions.h>

#include "core/DataTypes.h"

struct Triangle;
struct Vertex;
class BoundingBox;

class VF_PUBLIC STLReader
{
public:

    enum FileType { ascii, binary };

    static std::vector<Triangle> readSTL(const std::string& name);
    static std::vector<Triangle> readSTL(const std::string& name, FileType fileType, const std::vector<uint> ignorePatches = std::vector<uint>() );
	static std::vector<Triangle> readSTL(const BoundingBox &box, const std::string& name);

    static std::vector<Triangle> readBinarySTL(const std::string& name);
    static std::vector<Triangle> readASCIISTL(const std::string& name);
    static std::vector<Triangle> readASCIISTLWithPatches(const std::string& name, const std::vector<uint> ignorePatches);
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

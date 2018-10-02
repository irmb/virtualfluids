#define _CRT_SECURE_NO_DEPRECATE
#include "STLWriter.h"
#include <fstream>
#include <sstream>

#include <GridGenerator/geometries/Vertex/Vertex.h>
#include <GridGenerator/geometries/Triangle/Triangle.h>

#include <utilities/logger/Logger.h>


void STLWriter::writeSTL(std::vector<Triangle> &vec, const std::string &name, bool writeBinary)
{
    const int size = vec.size();
    *logging::out << logging::Logger::INFO_INTERMEDIATE << "Write " << size << " Triangles to STL : " + name + "\n";

    std::ios_base::openmode mode = std::ios::out;
    if (writeBinary)
        mode = std::ios::out | std::ios::binary;

    std::ofstream ofstream(name, mode);

    if (!ofstream.is_open()) {
        *logging::out << logging::Logger::INFO_HIGH << " Output file not open - exit function\n";
        return;
    }

    if (writeBinary)
        writeBinarySTL(ofstream, vec);
    else
        writeAsciiSTL(ofstream, vec);

    ofstream.close();
    *logging::out << logging::Logger::INFO_INTERMEDIATE << "Output file closed\n";
}


void STLWriter::writeAsciiSTL(std::ofstream &ofstream, std::vector<Triangle> &vec)
{
    ofstream << "solid ascii\n";
    for (size_t i = 0; i < vec.size(); i++) 
    {
        Triangle t = vec[i];

        ofstream << "facet normal ";
        t.normal.printFormatted(ofstream);
        ofstream << "\n";

        ofstream << "outer loop\n";

        ofstream << "vertex ";
        t.v1.printFormatted(ofstream);
        ofstream << "\n";
        ofstream << "vertex ";
        t.v2.printFormatted(ofstream);
        ofstream << "\n";
        ofstream << "vertex ";
        t.v3.printFormatted(ofstream);
        ofstream << "\n";

        ofstream << "endloop\n";
        ofstream << "endfacet\n";
    }
    ofstream << "endsolid\n";
}

void STLWriter::writeBinarySTL(std::ofstream &ofstream, std::vector<Triangle> &vec)
{
    char header_info[80] = "GridGeneration-File iRMB";
    unsigned long nTriLong = (unsigned long)vec.size();
    ofstream.write(header_info, sizeof(header_info));
    ofstream.write((char*)&nTriLong, 4);
    char attribute[2] = "0";
    for (size_t i = 0; i < vec.size(); i++)
    {
        Triangle t = vec[i];

        t.normal.print(ofstream);

        t.v1.print(ofstream);
        t.v2.print(ofstream);
        t.v3.print(ofstream);

        ofstream.write(attribute, 2);
    }
}

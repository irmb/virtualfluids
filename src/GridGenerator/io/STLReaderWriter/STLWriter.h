#ifndef STLWriter_H
#define STLWriter_H



#include <vector>
#include <string>
#include <memory>
#include <fstream>

#include <VirtualFluidsDefinitions.h>

class Transformator;
struct Triangle;

class VF_PUBLIC STLWriter
{
public:
    static void writeSTL(std::vector<Triangle> &vec, const std::string &name, bool writeBinary = false);

private:
    STLWriter() {};
    STLWriter(const STLWriter &) {};
    virtual ~STLWriter() {};

    static void writeAsciiSTL(std::ofstream &ofstream, std::vector<Triangle> &vec);
    static void writeBinarySTL(std::ofstream &ofstream, std::vector<Triangle> &vec);
};

#endif

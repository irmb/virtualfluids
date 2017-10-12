#ifndef STLWriter_H
#define STLWriter_H

#include "GridGenerator_EXPORT.h"

#include <vector>
#include <string>
#include <memory>
#include <fstream>

class Transformator;
struct Triangle;

class GridGenerator_EXPORT STLWriter
{
public:
    static void writeSTL(std::vector<Triangle> &vec, const std::string &name, std::shared_ptr<const Transformator> trans, bool writeBinary = false);

private:
    STLWriter() {};
    STLWriter(const STLWriter &) {};
    virtual ~STLWriter() {};

    static void writeAsciiSTL(std::ofstream &ofstream, std::vector<Triangle> &vec, std::shared_ptr<const Transformator> trans);
    static void writeBinarySTL(std::ofstream &ofstream, std::vector<Triangle> &vec, std::shared_ptr<const Transformator> trans);
};

#endif

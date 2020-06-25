#ifndef QLINEWRITER_H
#define QLINEWRITER_H

#include <vector>
#include <string>
#include <fstream>

#include "global.h"

class GeometryBoundaryCondition;
class Grid;
struct Vertex;


class QLineWriter 
{
public:
    static void writeArrows(std::string fileName, SPtr<GeometryBoundaryCondition> geometryBoundaryCondition, SPtr<Grid> grid);

private:
    static Vertex getVertex(int matrixIndex, SPtr<Grid> grid);
};


#endif

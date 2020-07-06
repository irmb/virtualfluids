#ifndef VTKGEOMETRYREADER_H
#define VTKGEOMETRYREADER_H

#include <string>

/**
* @file
* @author  Kostyantyn Kucher <k.kucher@gmx.net>
* @version 1.0
*
* @section LICENSE
*
* @section DESCRIPTION
*
* Read geometry from VTK Lagacy File 
* geometry is integer values for node types difinition
* respectively
* SOLID = 
* FLUID = 
* 
*/

class VtkGeometryReader
{
public:
   VtkGeometryReader(){}
   ~VtkGeometryReader(){}
   static void readFile(const std::string& fileName, unsigned int* geoMat);
protected:
private:
};

#endif
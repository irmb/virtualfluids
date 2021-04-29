#ifndef READ_GEOMETRY_H
#define READ_GEOMETRY_H

#include "Parameter/Parameter.h"
#include "Communication/Communicator.h"
#include "Input/VtkGeometryReader.h"

#include <iostream>

extern "C" void readGeometry(Parameter* para, VF::GPU::Communicator* comm, int lev, std::string geometryFile);

#endif

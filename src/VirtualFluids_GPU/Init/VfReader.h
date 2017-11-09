#ifndef VF_READER_H
#define VF_READER_H

#include "Parameter/Parameter.h"
#include "Communication/Communicator.h"
#include "Input/kFullReader.h"
#include "Input/PositionReader.h"
//#include "GPU/GPU_Interface.h"

#include <iostream>

extern "C" void readVFkFull(Parameter* para, const std::string geometryFile);

extern "C" void readVFgeoFull(Parameter* para, const std::string geometryFile);

extern "C" void readVecSP(Parameter* para);

extern "C" void readInterfaceCF(Parameter* para);

extern "C" void readInterfaceFC(Parameter* para);

extern "C" void readInterfaceOffCF(Parameter* para, const std::string geometryFile);

extern "C" void readInterfaceOffFC(Parameter* para, const std::string geometryFile);

extern "C" void readNoSlipBc(Parameter* para);

extern "C" void readSlipBc(Parameter* para);

extern "C" void readPressBc(Parameter* para);

extern "C" void readPropellerCylinder(Parameter* para);

extern "C" void readMeasurePoints(Parameter* para);

#endif

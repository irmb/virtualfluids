#ifndef WRITE_DATA_H
#define WRITE_DATA_H

#include "Parameter/Parameter.h"
#include "Utilities/StringUtil.hpp"
#include "Communication/Communicator.h"

#include <iostream>

extern "C" void writeInit(Parameter* para);
extern "C" void writeTimestep(Parameter* para, unsigned int t);
extern "C" void writeParticle(Parameter* para, unsigned int t);

#endif

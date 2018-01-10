#ifndef WRITE_DATA_H
#define WRITE_DATA_H

#include <core/PointerDefinitions.h>

class Parameter;

void writeInit(SPtr<Parameter> para);
void writeTimestep(Parameter* para, unsigned int t);
void writeParticle(Parameter* para, unsigned int t);

#endif

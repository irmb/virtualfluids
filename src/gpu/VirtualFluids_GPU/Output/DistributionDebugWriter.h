#ifndef DISTRIBUTION_DEBUG_WRITER
#define DISTRIBUTION_DEBUG_WRITER

#include <basics/DataTypes.h>

class Parameter;

class DistributionDebugWriter
{
public:
    static void writeDistributions(const Parameter* para, uint timestep);
    static void writeDistributionsForLevel(const Parameter* para, uint level, uint timestep);
};

#endif

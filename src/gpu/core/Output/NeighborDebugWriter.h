#ifndef NEIGHBORDEBUG_HPP
#define NEIGHBORDEBUG_HPP

#include <string>
class Parameter;
struct LBMSimulationParameter;
class WbWriter;

class NeighborDebugWriter
{
public:
    static void writeNeighborLinkLinesDebug(Parameter *para);

protected:
    static void writeNeighborLinkLines(LBMSimulationParameter *parH, int direction, const std::string &name,
                                       WbWriter *writer);
};

#endif

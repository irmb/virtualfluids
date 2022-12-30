#ifndef NEIGHBORDEBUG_HPP
#define NEIGHBORDEBUG_HPP

#include "LBM/LB.h"
#include "Logger.h"
#include "Parameter/Parameter.h"
#include "basics/utilities/UbSystem.h"
#include "grid/NodeValues.h"
#include "lbm/constants/D3Q27.h"
#include <basics/writer/WbWriterVtkXmlBinary.h>

#include "Utilities/FindNeighbors.h"
#include "VirtualFluids_GPU/Communication/Communicator.h"
#include "Core/StringUtilities/StringUtil.h"

namespace NeighborDebugWriter
{

inline void writeNeighborLinkLines(Parameter *para, const int level, const unsigned long long numberOfNodes, const int direction,
                                   const std::string &name)
{
    VF_LOG_INFO("Write node links in direction {}.", direction);
    std::vector<UbTupleFloat3> nodes(numberOfNodes * 2);
    std::vector<UbTupleInt2> cells(numberOfNodes);

    for (size_t position = 0; position < numberOfNodes; position++) {
        if (para->getParH(level)->typeOfGridNode[position] != GEO_FLUID)
            continue;

        const double x1 = para->getParH(level)->coordinateX[position];
        const double x2 = para->getParH(level)->coordinateY[position];
        const double x3 = para->getParH(level)->coordinateZ[position];

        const uint positionNeighbor = getNeighborIndex(para->getParH(level).get(), position, direction);

        const double x1Neighbor = para->getParH(level)->coordinateX[positionNeighbor];
        const double x2Neighbor = para->getParH(level)->coordinateY[positionNeighbor];
        const double x3Neighbor = para->getParH(level)->coordinateZ[positionNeighbor];

        nodes.emplace_back(float(x1), float(x2), float(x3));
        nodes.emplace_back(float(x1Neighbor), float(x2Neighbor), float(x3Neighbor));

        cells.emplace_back((int)nodes.size() - 2, (int)nodes.size() - 1);
    }
    WbWriterVtkXmlBinary::getInstance()->writeLines(name, nodes, cells);
}

inline void writeNeighborLinkLinesDebug(Parameter *para)
{
    for (int level = 0; level <= para->getMaxLevel(); level++) {
        for (int direction = vf::lbm::dir::STARTDIR; direction <= vf::lbm::dir::ENDDIR; direction++) {
            const std::string fileName = para->getFName() + "_" + StringUtil::toString<int>(level) + "_Link_" +
                                         std::to_string(direction) + "_Debug.vtk";
            writeNeighborLinkLines(para, level, para->getParH(level)->numberOfNodes, direction, fileName);
        }
    }
}

} // namespace NeighborDebugWriter

#endif

#ifndef EDGENODEDEBUG_HPP
#define EDGENODEDEBUG_HPP

#include <fstream>
#include <sstream>
#include <stdio.h>
// #include <math.h>
#include "Core/StringUtilities/StringUtil.h"
#include "LBM/D3Q27.h"
#include "LBM/LB.h"
#include "Parameter/Parameter.h"
#include "basics/utilities/UbSystem.h"
#include <basics/writer/WbWriterVtkXmlBinary.h>
#include <cmath>

#include "VirtualFluids_GPU/Communication/Communicator.h"

namespace EdgeNodeDebugWriter
{

void writeEdgeNodesXZ_Send(Parameter *para)
{
    std::vector<UbTupleFloat3> nodesVec;

    // nodedata
    std::vector<std::string> datanames = { "SparseIndex", "ProcessNeighbor", "IndexInSendVector", "AfterFtoC" };
    std::vector<std::vector<double>> nodedata;

    int numberOfNodes = 0;
    for (int level = 0; level < para->getMaxLevel(); level++){
        numberOfNodes += (int) para->getParH(level)->edgeNodesXtoZ.size();
    }

    nodesVec.resize(numberOfNodes);
    nodedata.resize(datanames.size(), std::vector<double>(numberOfNodes));

    int nodeCount = 0;
    for (int level = 0; level < para->getMaxLevel(); level++) {
        for (int u = 0; u < numberOfNodes; u++) {
            // node data section
            int indexOfProcessNeighborSend = para->getParH(level)->edgeNodesXtoZ[u].indexOfProcessNeighborSend;
            int indexInSendBuffer = para->getParH(level)->edgeNodesXtoZ[u].indexInSendBuffer;
            int sparseIndex = para->getParH(level)->sendProcessNeighborZ[indexOfProcessNeighborSend].index[indexInSendBuffer];
            nodedata[0][nodeCount] = sparseIndex;
            nodedata[1][nodeCount] = indexOfProcessNeighborSend;
            nodedata[2][nodeCount] = indexInSendBuffer;
            nodedata[3][nodeCount] = indexInSendBuffer < para->getParH(level)->sendProcessNeighborsAfterFtoCZ[indexOfProcessNeighborSend].numberOfNodes;

            // coordinate section
            double x1           = para->getParH(level)->coordX_SP[sparseIndex];
            double x2           = para->getParH(level)->coordY_SP[sparseIndex];
            double x3           = para->getParH(level)->coordZ_SP[sparseIndex];
            nodesVec[nodeCount] = (makeUbTuple((float)(x1), (float)(x2), (float)(x3)));

            nodeCount++;
        }
        std::string filenameVec = para->getFName() + "_writeEdgeNodesXZ_Send_PID_" +
                                  std::to_string(vf::gpu::Communicator::getInstanz()->getPID()) + "_" +
                                  StringUtil::toString<int>(level);

        WbWriterVtkXmlBinary::getInstance()->writeNodesWithNodeData(filenameVec, nodesVec, datanames, nodedata);
    }
}
} // namespace EdgeNodeDebugWriter

#endif

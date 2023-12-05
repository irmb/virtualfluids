#ifndef EDGENODEDEBUG_HPP
#define EDGENODEDEBUG_HPP

#include <fstream>
#include <sstream>
#include <cstdio>
// #include <math.h>
#include "StringUtilities/StringUtil.h"
#include "lbm/constants/D3Q27.h"
#include "Calculation/Calculation.h"
#include "Parameter/Parameter.h"
#include "basics/utilities/UbSystem.h"
#include <basics/writer/WbWriterVtkXmlBinary.h>
#include <cmath>


namespace EdgeNodeDebugWriter
{

void addCoordinatesToNodeVector(SPtr<LBMSimulationParameter> parH, std::vector<UbTupleFloat3> &nodesVec, int indexInNodesVector, int sparseIndexOfNode){
            double x1           = parH->coordinateX[sparseIndexOfNode];
            double x2           = parH->coordinateY[sparseIndexOfNode];
            double x3           = parH->coordinateZ[sparseIndexOfNode];
            nodesVec[indexInNodesVector] = (makeUbTuple((float)(x1), (float)(x2), (float)(x3)));
}

void writeEdgeNodesXZ_Send(SPtr<Parameter> para, int processID = 0)
{
    std::vector<UbTupleFloat3> nodesVec;
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
            int indexOfProcessNeighborSend = para->getParH(level)->edgeNodesXtoZ[u].indexOfProcessNeighborSend;
            int indexInSendBuffer = para->getParH(level)->edgeNodesXtoZ[u].indexInSendBuffer;
            int sparseIndex = para->getParH(level)->sendProcessNeighborZ[indexOfProcessNeighborSend].index[indexInSendBuffer];
            nodedata[0][nodeCount] = sparseIndex;
            nodedata[1][nodeCount] = indexOfProcessNeighborSend;
            nodedata[2][nodeCount] = indexInSendBuffer;
            nodedata[3][nodeCount] = indexInSendBuffer < para->getParH(level)->sendProcessNeighborsAfterFtoCZ[indexOfProcessNeighborSend].numberOfNodes;

            addCoordinatesToNodeVector(para->getParH(level), nodesVec, nodeCount, sparseIndex);

            nodeCount++;
        }
        std::string filenameVec = para->getFName() + "_writeEdgeNodesXZ_Send_PID_" +
                                  std::to_string(processID) + "_" +
                                  StringUtil::toString<int>(level);

        WbWriterVtkXmlBinary::getInstance()->writeNodesWithNodeData(filenameVec, nodesVec, datanames, nodedata);
    }
}

void writeEdgeNodesXZ_Recv(SPtr<Parameter> para, int processID = 0)
{
    std::vector<UbTupleFloat3> nodesVec;
    std::vector<std::string> datanames = { "SparseIndex", "ProcessNeighbor", "IndexInRecvVector", "AfterFtoC" };
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
            int indexOfProcessNeighborRecv = para->getParH(level)->edgeNodesXtoZ[u].indexOfProcessNeighborRecv;
            int indexInRecvBuffer = para->getParH(level)->edgeNodesXtoZ[u].indexInRecvBuffer;
            int sparseIndex = para->getParH(level)->recvProcessNeighborX[indexOfProcessNeighborRecv].index[indexInRecvBuffer];
            nodedata[0][nodeCount] = sparseIndex;
            nodedata[1][nodeCount] = indexOfProcessNeighborRecv;
            nodedata[2][nodeCount] = indexInRecvBuffer;
            nodedata[3][nodeCount] = indexInRecvBuffer < para->getParH(level)->recvProcessNeighborX[indexOfProcessNeighborRecv].numberOfNodes;

            addCoordinatesToNodeVector(para->getParH(level), nodesVec, nodeCount, sparseIndex);

            nodeCount++;
        }
        std::string filenameVec = para->getFName() + "_writeEdgeNodesXZ_Recv_PID_" +
                                  std::to_string(processID) + "_" +
                                  StringUtil::toString<int>(level);

        WbWriterVtkXmlBinary::getInstance()->writeNodesWithNodeData(filenameVec, nodesVec, datanames, nodedata);
    }
}
} // namespace EdgeNodeDebugWriter

#endif

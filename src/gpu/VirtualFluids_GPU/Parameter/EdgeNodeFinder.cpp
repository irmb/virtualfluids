#include "EdgeNodeFinder.h"
#include "Parameter.h"

namespace vf::gpu
{
void findEdgeNodesCommMultiGPU(SPtr<Parameter> parameter)
{
    for (int level = 0; level <= parameter->getFine(); level++) {
        findEdgeNodesXY(level, parameter);
        findEdgeNodesXZ(level, parameter);
        findEdgeNodesYZ(level, parameter);
    }
}
} // namespace vf::gpu

namespace
{
void findEdgeNodesXY(int level, SPtr<Parameter> parameter)
{
    int indexOfProcessNeighborSend;
    int indexInSendBuffer;
    for (uint i = 0; i < (unsigned int)(parameter->getNumberOfProcessNeighborsX(level, "recv")); i++) {
        for (int j = 0; j < parameter->getParH(level)->recvProcessNeighborX[i].numberOfNodes; j++) {
            int nodeIndex = parameter->getParH(level)->recvProcessNeighborX[i].index[j];
            bool foundIndex =
                findIndexInSendNodes(nodeIndex, parameter->getParH(level)->sendProcessNeighborY, indexOfProcessNeighborSend,  indexInSendBuffer);
            if (foundIndex) {
                parameter->getParH(level)->edgeNodesXtoY.emplace_back(i, j, indexOfProcessNeighborSend,
                                                                      indexInSendBuffer);
            }
        }
    }
}

void findEdgeNodesXZ(int level, SPtr<Parameter> parameter)
{
    int indexOfProcessNeighborSend;
    int indexInSendBuffer;
    for (uint i = 0; i < (unsigned int)(parameter->getNumberOfProcessNeighborsX(level, "recv")); i++) {
        for (int j = 0; j < parameter->getParH(level)->recvProcessNeighborX[i].numberOfNodes; j++) {
            int nodeIndex = parameter->getParH(level)->recvProcessNeighborX[i].index[j];
            bool foundIndex =
                findIndexInSendNodes(nodeIndex, parameter->getParH(level)->sendProcessNeighborZ, indexOfProcessNeighborSend,  indexInSendBuffer);
            if (foundIndex) {
                parameter->getParH(level)->edgeNodesXtoZ.emplace_back(i, j, indexOfProcessNeighborSend,
                                                                      indexInSendBuffer);
            }
        }
    }
}

void findEdgeNodesYZ(int level, SPtr<Parameter> parameter)
{
    int indexOfProcessNeighborSend;
    int indexInSendBuffer;
    for (uint i = 0; i < (unsigned int)(parameter->getNumberOfProcessNeighborsY(level, "recv")); i++) {
        for (int j = 0; j < parameter->getParH(level)->recvProcessNeighborY[i].numberOfNodes; j++) {
            int nodeIndex = parameter->getParH(level)->recvProcessNeighborY[i].index[j];
            bool foundIndex =
                findIndexInSendNodes(nodeIndex, parameter->getParH(level)->sendProcessNeighborZ,indexOfProcessNeighborSend,  indexInSendBuffer);
            if (foundIndex) {
                parameter->getParH(level)->edgeNodesYtoZ.emplace_back(i, j, indexOfProcessNeighborSend,
                                                                      indexInSendBuffer);
            }
        }
    }
}

bool findIndexInSendNodes(int nodeIndex, const std::vector<ProcessNeighbor27>& sendProcessNeighbor, int &indexOfProcessNeighborSend, int &indexInSendBuffer)
{
    for (uint neighbor = 0; neighbor < (unsigned int)sendProcessNeighbor.size(); neighbor++) {
        for (int node = 0; node < sendProcessNeighbor[neighbor].numberOfNodes; node++) {
            if (sendProcessNeighbor[neighbor].index[node] == nodeIndex) {
                indexOfProcessNeighborSend = neighbor;
                indexInSendBuffer = node;
                return true;
            }
        }
    }
    return false;
}

} // namespace
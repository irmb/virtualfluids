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
            int index = parameter->getParH(level)->recvProcessNeighborX[i].index[j];
            bool foundIndex =
                findIndexInSendNodesXY(level, index, indexOfProcessNeighborSend, indexInSendBuffer, parameter);
            if (foundIndex) {
                parameter->getParH(level)->edgeNodesXtoY.emplace_back(i, j, indexOfProcessNeighborSend,
                                                                      indexInSendBuffer);
            }
        }
    }
}

bool findIndexInSendNodesXY(int level, int index, int &indexOfProcessNeighborSend, int &indexInSendBuffer,
                            SPtr<Parameter> parameter)
{
    for (uint k = 0; k < (unsigned int)(parameter->getNumberOfProcessNeighborsY(level, "send")); k++) {
        for (int l = 0; l < parameter->getParH(level)->sendProcessNeighborY[k].numberOfNodes; l++) {
            if (parameter->getParH(level)->sendProcessNeighborY[k].index[l] == index) {
                indexOfProcessNeighborSend = k;
                indexInSendBuffer = l;
                return true;
            }
        }
    }
    return false;
}

void findEdgeNodesXZ(int level, SPtr<Parameter> parameter)
{
    int indexOfProcessNeighborSend;
    int indexInSendBuffer;
    for (uint i = 0; i < (unsigned int)(parameter->getNumberOfProcessNeighborsX(level, "recv")); i++) {
        for (int j = 0; j < parameter->getParH(level)->recvProcessNeighborX[i].numberOfNodes; j++) {
            int index = parameter->getParH(level)->recvProcessNeighborX[i].index[j];
            bool foundIndex =
                findIndexInSendNodesXZ(level, index, indexOfProcessNeighborSend, indexInSendBuffer, parameter);
            if (foundIndex) {
                parameter->getParH(level)->edgeNodesXtoZ.emplace_back(i, j, indexOfProcessNeighborSend,
                                                                      indexInSendBuffer);
            }
        }
    }
}

bool findIndexInSendNodesXZ(int level, int index, int &indexOfProcessNeighborSend, int &indexInSendBuffer,
                            SPtr<Parameter> parameter)
{
    for (uint k = 0; k < (unsigned int)(parameter->getNumberOfProcessNeighborsZ(level, "send")); k++) {
        for (int l = 0; l < parameter->getParH(level)->sendProcessNeighborZ[k].numberOfNodes; l++) {
            if (parameter->getParH(level)->sendProcessNeighborZ[k].index[l] == index) {
                indexOfProcessNeighborSend = k;
                indexInSendBuffer = l;
                return true;
            }
        }
    }
    return false;
}

void findEdgeNodesYZ(int level, SPtr<Parameter> parameter)
{
    int indexOfProcessNeighborSend;
    int indexInSendBuffer;
    for (uint i = 0; i < (unsigned int)(parameter->getNumberOfProcessNeighborsY(level, "recv")); i++) {
        for (int j = 0; j < parameter->getParH(level)->recvProcessNeighborY[i].numberOfNodes; j++) {
            int index = parameter->getParH(level)->recvProcessNeighborY[i].index[j];
            bool foundIndex =
                findIndexInSendNodesYZ(level, index, indexOfProcessNeighborSend, indexInSendBuffer, parameter);
            if (foundIndex) {
                parameter->getParH(level)->edgeNodesYtoZ.emplace_back(i, j, indexOfProcessNeighborSend,
                                                                      indexInSendBuffer);
            }
        }
    }
}

bool findIndexInSendNodesYZ(int level, int index, int &indexOfProcessNeighborSend, int &indexInSendBuffer,
                            SPtr<Parameter> parameter)
{
    for (uint k = 0; k < (unsigned int)(parameter->getNumberOfProcessNeighborsZ(level, "send")); k++) {
        for (int l = 0; l < parameter->getParH(level)->sendProcessNeighborZ[k].numberOfNodes; l++) {
            if (parameter->getParH(level)->sendProcessNeighborZ[k].index[l] == index) {
                indexOfProcessNeighborSend = k;
                indexInSendBuffer = l;
                return true;
            }
        }
    }
    return false;
}
} // namespace
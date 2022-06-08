#include <vector>

#include "EdgeNodeFinder.h"
#include "Parameter.h"

namespace vf::gpu
{
//! \brief Find nodes that are both received in the x-direction and sent in the y-direction
void findEdgeNodesXY(int level, SPtr<Parameter> parameter);
//! \brief Find nodes that are both received in the x-direction and sent in the z-direction
void findEdgeNodesXZ(int level, SPtr<Parameter> parameter);
//! \brief Find nodes that are both received in the y-direction and sent in the z-direction
void findEdgeNodesYZ(int level, SPtr<Parameter> parameter);
void findEdgeNodes(const std::vector<ProcessNeighbor27> &recvProcessNeighbor,
                   const std::vector<ProcessNeighbor27> &sendProcessNeighbor,
                   std::vector<LBMSimulationParameter::EdgeNodePositions> &edgeNodes);
bool findIndexInSendNodes(int nodeIndex, const std::vector<ProcessNeighbor27> &sendProcessNeighbor,
                          int &indexOfProcessNeighborSend, int &indexInSendBuffer);

void findEdgeNodesCommMultiGPU(SPtr<Parameter> parameter)
{
    for (int level = 0; level <= parameter->getFine(); level++) {
        findEdgeNodesXY(level, parameter);
        findEdgeNodesXZ(level, parameter);
        findEdgeNodesYZ(level, parameter);
    }
}

void findEdgeNodesXY(int level, SPtr<Parameter> parameter)
{
    findEdgeNodes(parameter->getParH(level)->recvProcessNeighborX, parameter->getParH(level)->sendProcessNeighborY,
                  parameter->getParH(level)->edgeNodesXtoY);
}

void findEdgeNodesXZ(int level, SPtr<Parameter> parameter)
{
    findEdgeNodes(parameter->getParH(level)->recvProcessNeighborX, parameter->getParH(level)->sendProcessNeighborZ,
                  parameter->getParH(level)->edgeNodesXtoZ);
}

void findEdgeNodesYZ(int level, SPtr<Parameter> parameter)
{
    findEdgeNodes(parameter->getParH(level)->recvProcessNeighborY, parameter->getParH(level)->sendProcessNeighborZ,
                  parameter->getParH(level)->edgeNodesYtoZ);
}

void findEdgeNodes(const std::vector<ProcessNeighbor27> &recvProcessNeighbor,
                   const std::vector<ProcessNeighbor27> &sendProcessNeighbor,
                   std::vector<LBMSimulationParameter::EdgeNodePositions> &edgeNodes)
{
    int indexOfProcessNeighborSend;
    int indexInSendBuffer;
    for (uint i = 0; i < (unsigned int)(recvProcessNeighbor.size()); i++) {
        for (int j = 0; j < recvProcessNeighbor[i].numberOfNodes; j++) {
            int nodeIndex = recvProcessNeighbor[i].index[j];
            bool foundIndex =
                findIndexInSendNodes(nodeIndex, sendProcessNeighbor, indexOfProcessNeighborSend, indexInSendBuffer);
            if (foundIndex) {
                edgeNodes.emplace_back(i, j, indexOfProcessNeighborSend, indexInSendBuffer);
            }
        }
    }
}

bool findIndexInSendNodes(int nodeIndex, const std::vector<ProcessNeighbor27> &sendProcessNeighbor,
                          int &indexOfProcessNeighborSend, int &indexInSendBuffer)
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

} // namespace vf::gpu

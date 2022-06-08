#include <vector>
#include <optional>

#include "EdgeNodeFinder.h"
#include "Parameter.h"

namespace vf::gpu
{
//! \brief Find nodes that are both received in the x-direction and sent in the y-direction
void findEdgeNodesXY(const int level, LBMSimulationParameter& parameterLB);
//! \brief Find nodes that are both received in the x-direction and sent in the z-direction
void findEdgeNodesXZ(const int level, LBMSimulationParameter& parameterLB);
//! \brief Find nodes that are both received in the y-direction and sent in the z-direction
void findEdgeNodesYZ(const int level, LBMSimulationParameter& parameterLB);
void findEdgeNodes(const std::vector<ProcessNeighbor27> &recvProcessNeighbor,
                   const std::vector<ProcessNeighbor27> &sendProcessNeighbor,
                   std::vector<LBMSimulationParameter::EdgeNodePositions> &edgeNodes);
std::optional<std::pair<int, int>> findIndexInSendNodes(const int nodeIndex,
                                                        const std::vector<ProcessNeighbor27> &sendProcessNeighbor);

void findEdgeNodesCommMultiGPU(Parameter& parameter)
{
    for (int level = 0; level <= parameter.getFine(); level++) {
        findEdgeNodesXY(level, *parameter.getParH(level));
        findEdgeNodesXZ(level, *parameter.getParH(level));
        findEdgeNodesYZ(level, *parameter.getParH(level));
    }
}

void findEdgeNodesXY(const int level, LBMSimulationParameter& parameterLB)
{
    findEdgeNodes(parameterLB.recvProcessNeighborX, parameterLB.sendProcessNeighborY,
                  parameterLB.edgeNodesXtoY);
}

void findEdgeNodesXZ(const int level, LBMSimulationParameter& parameterLB)
{
    findEdgeNodes(parameterLB.recvProcessNeighborX, parameterLB.sendProcessNeighborZ,
                  parameterLB.edgeNodesXtoZ);
}

void findEdgeNodesYZ(const int level, LBMSimulationParameter& parameterLB)
{
    findEdgeNodes(parameterLB.recvProcessNeighborY, parameterLB.sendProcessNeighborZ,
                  parameterLB.edgeNodesYtoZ);
}

void findEdgeNodes(const std::vector<ProcessNeighbor27> &recvProcessNeighbor,
                   const std::vector<ProcessNeighbor27> &sendProcessNeighbor,
                   std::vector<LBMSimulationParameter::EdgeNodePositions> &edgeNodes)
{
    for (uint neighbor = 0; neighbor < (unsigned int)(recvProcessNeighbor.size()); neighbor++) {
        for (int index = 0; index < recvProcessNeighbor[neighbor].numberOfNodes; index++) {
            if (auto sendIndices = findIndexInSendNodes(recvProcessNeighbor[neighbor].index[index], sendProcessNeighbor)) {
                edgeNodes.emplace_back(neighbor, index, sendIndices->first, sendIndices->second);
            }
        }
    }
}

std::optional<std::pair<int, int>> findIndexInSendNodes(const int nodeIndex,
                                                        const std::vector<ProcessNeighbor27> &sendProcessNeighbor)
{
    for (uint neighbor = 0; neighbor < (unsigned int)sendProcessNeighbor.size(); neighbor++) {
        for (int node = 0; node < sendProcessNeighbor[neighbor].numberOfNodes; node++) {
            if (sendProcessNeighbor[neighbor].index[node] == nodeIndex) {
                return std::pair<int, int>(neighbor, node);
            }
        }
    }
    return std::nullopt;
}

} // namespace vf::gpu

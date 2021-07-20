#include "NodeIndicesMultiGPU.h"

NodeIndicesMultiGPU::NodeIndicesMultiGPU(const std::vector<uint> *geoFluidSizes,
                                         const std::vector<const uint*> *geoFluidNodesIndices)
{
    this->geoFluidSizes = geoFluidSizes;
    this->geoFluidNodeIndices = geoFluidNodesIndices;
}

uint NodeIndicesMultiGPU::getGeoFluidSize(uint level) 
{ 
    return (*this->geoFluidSizes)[level];
}

const uint* NodeIndicesMultiGPU::getGeoFluidNodeIndices(uint level)
{
    return (*this->geoFluidNodeIndices)[level];
}

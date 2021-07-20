#include "NodeIndicesMultiGPU.h"

NodeIndicesMultiGPU::NodeIndicesMultiGPU(const std::vector<uint> *geoFluidSizes,
                                         const std::vector<const uint*> *geoFluidNodesIndices)
{
    this->geoFluidSizes = geoFluidSizes;
    this->geoFluidNodeIndices = geoFluidNodesIndices;
}

uint NodeIndicesMultiGPU::getGeoFluidSize(uint gridNumber) 
{ 
    return (*this->geoFluidSizes)[gridNumber];
}

const uint* NodeIndicesMultiGPU::getGeoFluidNodeIndices(uint gridNumber)
{
    return (*this->geoFluidNodeIndices)[gridNumber];
}

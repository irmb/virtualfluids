#include "NodeIndicesMultiGPU.h"

NodeIndicesMultiGPU::NodeIndicesMultiGPU(const std::vector<uint> *geoFluidSizes,
                                         const std::vector<const std::vector<uint> *> *geoFluidNodesIndices)
{
    this->geoFluidSizes = geoFluidSizes;
    this->geoFluidNodeIndices = geoFluidNodesIndices;
}

uint NodeIndicesMultiGPU::getGeoFluidSize(uint gridNumber) 
{ 
    return (*this->geoFluidSizes)[gridNumber];
}

const std::vector<uint>* NodeIndicesMultiGPU::getGeoFluidNodeIndices(uint gridNumber)
{
    return (*this->geoFluidNodeIndices)[gridNumber];
}

#include "NodeIndicesMultiGPU.h"

NodeIndicesMultiGPU::NodeIndicesMultiGPU(SPtr<MultipleGridBuilder> gridBuilder)
{
    std::vector<std::shared_ptr<Grid>> grids = gridBuilder->getGrids();
    for (uint i = 0; i < grids.size(); i++) {
        geoFluidSize.push_back(grids[i]->getGeoFluidSize());
        geoFluidNodeIndices.push_back(grids[i]->getGeoFluidNodes());
    }
}

uint NodeIndicesMultiGPU::getGeoFluidSize(uint gridNumber) 
{ 
    return this->geoFluidSize[gridNumber]; 
}

const std::vector<uint>* NodeIndicesMultiGPU::getGeoFluidNodeIndices(uint gridNumber)
{
    return this->geoFluidNodeIndices[gridNumber];
}

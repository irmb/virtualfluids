#include <vector>
#include <memory>

#include "GridGenerator/grid/GridBuilder/MultipleGridBuilder.h"

class NodeIndicesMultiGPU
{
    std::vector<uint> geoFluidSize;
    std::vector<const std::vector<uint>*> geoFluidNodeIndices; 

public:
    NodeIndicesMultiGPU(SPtr<MultipleGridBuilder> gridBuilder);

    uint getGeoFluidSize(uint gridNumber);
    const std::vector<uint>* getGeoFluidNodeIndices(uint gridNumber);
};
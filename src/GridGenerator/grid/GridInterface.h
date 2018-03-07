#ifndef GRID_INTERFACE_H
#define GRID_INTERFACE_H

#include "core/DataTypes.h"
#include "VirtualFluidsDefinitions.h"
#include "utilities/cuda/cudaDefines.h"



class GridImp;

class GridInterface
{
public:
    HOSTDEVICE VF_PUBLIC GridInterface();
    HOSTDEVICE VF_PUBLIC ~GridInterface();

    HOSTDEVICE void VF_PUBLIC findInterfaceCF(const uint& indexOnCoarseGrid, GridImp* coarseGrid, GridImp* fineGrid);
    HOSTDEVICE void VF_PUBLIC findInterfaceFC(const uint& indexOnCoarseGrid, GridImp* coarseGrid, GridImp* fineGrid);
    HOSTDEVICE void VF_PUBLIC findOverlapStopper(const uint& indexOnCoarseGrid, GridImp* coarseGrid, GridImp* fineGrid);

    HOSTDEVICE void VF_PUBLIC print() const;

    struct Interface
    {
        uint *fine, *coarse;
        uint numberOfEntries = 0;
    } fc, cf;


private:
    HOSTDEVICE int getCoarseToFineIndexOnFineGrid(const uint& indexOnCoarseGrid, const GridImp* coarseGrid, const GridImp* fineGrid);
    HOSTDEVICE bool isNeighborFineFluid(real x, real y, real z, const GridImp* coarseGrid, const GridImp* fineGrid);

    HOSTDEVICE int getFineToCoarseIndexOnFineGrid(const uint& indexOnCoarseGrid, const GridImp* coarseGrid, const GridImp* fineGrid);


};


#endif
#ifndef GRID_INTERFACE_H
#define GRID_INTERFACE_H

#include "core/DataTypes.h"
#include "VirtualFluidsDefinitions.h"
#include "utilities/cuda/cudaDefines.h"

struct Grid;

class GridInterface
{
public:
    HOSTDEVICE VF_PUBLIC GridInterface(const Grid* finerGrid);

    HOSTDEVICE VF_PUBLIC ~GridInterface();

    HOSTDEVICE void VF_PUBLIC findCF(const uint& index, const Grid* coarseGrid, const Grid* fineGrid);
    HOSTDEVICE void VF_PUBLIC findFC(const uint& index, const Grid* coarseGrid, const Grid* fineGrid);

    HOSTDEVICE void VF_PUBLIC print() const;

    struct Interface
    {
        uint *fine, *coarse;
        uint numberOfEntries = 0;
        real startCoarseX;
        real startCoarseY;
        real startCoarseZ;

        real endCoarseX;
        real endCoarseY;
        real endCoarseZ;

        char coarseEntry, fineEntry;
    } fc, cf;

private:
    HOSTDEVICE void initalCoarseToFine(uint sizeCF, const Grid* fineGrid);
    HOSTDEVICE void initalFineToCoarse(uint sizeCF, const Grid* fineGrid);

    HOSTDEVICE static void findInterface(Interface& interface, const int& factor, const uint& index, const Grid* coarseGrid, const Grid* fineGrid);
    HOSTDEVICE static bool isOnInterface(Interface& interface, const real& x, const real& y, const real& z);
    HOSTDEVICE static uint getIndexOnFinerGrid(const real& factor, const Grid* fineGrid, const real& x, const real& y, const real& z);
};


#endif
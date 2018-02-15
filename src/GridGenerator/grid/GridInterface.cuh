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

    HOSTDEVICE void initalGridInterface(const GridImp* fineGrid);

    HOSTDEVICE void VF_PUBLIC findCF(const uint& index, const GridImp* coarseGrid, const GridImp* fineGrid);
    HOSTDEVICE void VF_PUBLIC findFC(const uint& index, const GridImp* coarseGrid, const GridImp* fineGrid);

    HOSTDEVICE void VF_PUBLIC print() const;

    struct Interface
    {
        uint *fine, *coarse;
        uint numberOfEntries = 0;

        real startOffset;
        real endOffset;

        char coarseEntry, fineEntry;
    } fc, cf;


private:
    HOSTDEVICE void initalCoarseToFine(const GridImp* fineGrid);
    HOSTDEVICE void initalFineToCoarse(const GridImp* fineGrid);

    HOSTDEVICE static void findInterface(Interface& interface, const int& factor, const uint& index, const GridImp* coarseGrid, const GridImp* fineGrid);
    HOSTDEVICE static uint getIndexOnFinerGrid(const real& factor, const GridImp* fineGrid, const real& x, const real& y, const real& z);
};


#endif
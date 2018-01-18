#ifndef GRID_INTERFACE_H
#define GRID_INTERFACE_H

#include "core/DataTypes.h"
#include "VirtualFluidsDefinitions.h"

struct Grid;

struct GridInterface
{
    VF_PUBLIC GridInterface(const Grid* finerGrid);
    VF_PUBLIC GridInterface();
    void VF_PUBLIC findCF(uint index, const Grid* coarseGrid, const Grid* fineGrid);

private:
    real startCFCx;
    real startCFCy;
    real startCFCz;

    real endCFCx;
    real endCFCy;
    real endCFCz;

    uint *cfc, *cff, *fcf, *fcc;

    int numberOfEntriesInCF = 0;
};


#endif
#ifndef GRID_INTERFACE_H
#define GRID_INTERFACE_H

#include "Grid.cuh"

struct GridInterface
{
    GridInterface(const Grid& finerGrid);

    real startCFCx;
    real startCFCy;
    real startCFCz;

    real endCFCx;
    real endCFCy;
    real endCFCz;

    uint *cfc, *cff, *fcf, *fcc;
};


#endif
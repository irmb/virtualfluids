#ifndef FINDINTERFACE_H
#define FINDINTERFACE_H

#include "LBM/LB.h"
#include "lbm/constants/D3Q27.h"


void interpolation(InterpolationCellCoarseToFine &intCF, InterpolationCellFineToCoarse &intFC, 
                               unsigned int LxCoarse, unsigned int LyCoarse, unsigned int LzCoarse, 
                               unsigned int LxFine, unsigned int LyFine, unsigned int LzFine, 
                               unsigned int dNx, unsigned int dNy, unsigned int dNz, 
                               unsigned int *kCoarse, unsigned int *kFine, bool* needInterface,
                               InterpolationCellNeighbor &offCF, InterpolationCellNeighbor &offFC);

#endif

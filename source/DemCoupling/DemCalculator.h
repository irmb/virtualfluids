/*
*  Author: S. Peters
*  mail: peters@irmb.tu-bs.de
*/
#ifndef DEM_CALCULATOR_H
#define DEM_CALCULATOR_H

#include <memory>

#include "BasicCalculator.h"

class Grid3D;
class Synchronizer;
class DemCoProcessor;
class CalculationManager;

class DemCalculator : public BasicCalculator
{
public:
    DemCalculator(SPtr<Grid3D> grid, SPtr<UbScheduler> additionalGhostLayerUpdateScheduler, int numberOfTimeSteps);
    virtual ~DemCalculator() {}

    void calculate() override;

};

#endif

/*
*  Author: S. Peters
*  mail: peters@irmb.tu-bs.de
*/
#ifndef DEM_CALCULATOR_H
#define DEM_CALCULATOR_H

#include <memory>

#include <VirtualFluids/Grid/Calculator.h>

class Grid3D;
class Synchronizer;
class DemCoProcessor;
class CalculationManager;

class DemCalculator
    : public Calculator
{
public:
    DemCalculator(std::shared_ptr<Grid3D> grid, std::shared_ptr<Synchronizer> synchronizer, bool mainThread = true);
    virtual ~DemCalculator() {}

    void calculate(const double& endTime, std::shared_ptr<CalculationManager> cm, boost::exception_ptr& error) override;

};

#endif

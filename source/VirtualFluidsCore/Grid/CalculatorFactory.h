/*
*  Author: S. Peters
*  mail: peters@irmb.tu-bs.de
*/
#ifndef CalculatorFactory_H
#define CalculatorFactory_H

#include <memory>

class Calculator;
class UbScheduler;
class Synchronizer;
class Grid3D;

enum class CalculatorType
{
    HYBRID, MPI, OMP, PREPOSTBC, FETOL, DEM
};

class CalculatorFactory
{
public:
    explicit CalculatorFactory(std::shared_ptr<UbScheduler> visScheduler) : visScheduler(visScheduler) { }
    virtual ~CalculatorFactory() {}

    virtual std::shared_ptr<Calculator> makeCalculator(std::shared_ptr<Grid3D> grid, std::shared_ptr<Synchronizer> sync, bool isMainThread, CalculatorType type) = 0;

protected:
    std::shared_ptr<UbScheduler> visScheduler;
}; 

#endif 


#ifndef AverageValuesSimulationObserver_H
#define AverageValuesSimulationObserver_H

#include <PointerDefinitions.h>
#include <string>
#include <vector>

#include "SimulationObserver.h"
#include "LBMSystem.h"
#include "UbTuple.h"

class UbScheduler;
class WbWriter;
class Grid3D;
class Block3D;

//! \brief  Computes the time averaged mean velocity and RMS values and writes to parallel .vtk
//! \details writes at given time intervals specified in scheduler (s), does averaging according to scheduler (Avs) and
//! resets according to scheduler (rs).  <br>
//!  Computes  the time averaged mean velocity  \f$ u_{mean}=\frac{1}{N}\sum\limits_{i=1}^n u_{i} \f$  and RMS of
//!  fluctuations. You need to calculate a square root before plotting RMS. <br>
//
//! \author  Sonja Uphoff, Kostyantyn Kucher
// \f$ u_{mean}=\frac{1}{N}\sum\limits_{i=1}^n u_{i} \f$
class AverageValuesSimulationObserver : public SimulationObserver
{
public:
    AverageValuesSimulationObserver();
    AverageValuesSimulationObserver(SPtr<Grid3D> grid, const std::string &path, WbWriter *const writer, SPtr<UbScheduler> s,
                             SPtr<UbScheduler> Avs, SPtr<UbScheduler> rsMeans, SPtr<UbScheduler> rsRMS, bool restart);
    //! Make update
    void update(real step) override;
    //! Resets averaged velocity and RMS-values according to ResetSceduler
    void reset(real step);

protected:
    //! Prepare data and write in .vtk file
    void collectData(real step);
    //! Reset data
    void resetDataRMS(real step);
    void resetDataMeans(real step);
    //! prepare data
    void addData(const SPtr<Block3D> block);
    void clearData();
    //! Computes average and RMS values of macroscopic quantities
    void calculateAverageValues(real timeStep);
    ////! write .txt file spatial intergrated averaged value, fluctuation, porous features
    // void collectPlotDataZ(double step);
    ////! create txt file and write head line
    // void initPlotDataZ(double step);

private:
    std::vector<UbTupleFloat3> nodes;
    std::vector<UbTupleUInt8> cells;
    std::vector<std::string> datanames;
    std::vector<std::vector<real>> data;
    std::vector<std::vector<SPtr<Block3D>>> blockVector;
    int minInitLevel; // min init level
    int maxInitLevel;
    int gridRank;
    int resetStepRMS;
    int resetStepMeans;
    real averageInterval;
    std::string path;
    WbWriter *writer;
    bool restart, compressible;
    SPtr<UbScheduler> averageScheduler;    // additional scheduler to averaging after a given interval
    SPtr<UbScheduler> resetSchedulerRMS;   // additional scheduler to restart averaging after a given interval
    SPtr<UbScheduler> resetSchedulerMeans; // additional scheduler to restart averaging after a given interval
    // labels for the different components, e.g. AvVxx for time averaged RMS: 1/n SUM((U-Umean)^2)
    // you need to calculate a square root before plotting RMS
    enum Values {
        AvVx   = 0,
        AvVy   = 1,
        AvVz   = 2,
        AvVxx  = 3,
        AvVyy  = 4,
        AvVzz  = 5,
        AvVxy  = 6,
        AvVxz  = 7,
        AvVyz  = 8,
        AvP    = 9,
        AvPrms = 10
    };

    using CalcMacrosFct = void (*)(const real *const &, real &, real &, real &, real &);
    CalcMacrosFct calcMacros;
};
#endif

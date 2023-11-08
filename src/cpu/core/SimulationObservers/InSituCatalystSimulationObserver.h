#ifdef VF_CATALYST

#ifndef InSituCatalystSimulationObserver_h__
#define InSituCatalystSimulationObserver_h__

#include <SimulationObserver.h>
#include <Grid3D.h>
#include <LBMUnitConverter.h>
#include "lbm/constants/D3Q27.h"

#include <string>

#include <vtkCPDataDescription.h>
#include <vtkCPInputDataDescription.h>
#include <vtkCPProcessor.h>
#include <vtkCPPythonScriptPipeline.h>
#include <vtkDoubleArray.h>
#include <vtkNew.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>

class InSituCatalystSimulationObserver : public SimulationObserver
{
public:
    InSituCatalystSimulationObserver();
    InSituCatalystSimulationObserver(SPtr<Grid3D> grid, SPtr<UbScheduler> s, std::string script);
    virtual ~InSituCatalystSimulationObserver();
    void update(real step);

protected:
    void collectData(real step);
    void addData(SPtr<Block3D> block);
    void buildVTKGrid();
    void addVTKGridData(SPtr<Block3D> block);

private:
    std::vector<std::vector<SPtr<Block3D>>> blockVector;
    int minInitLevel;
    int maxInitLevel;
    int gridRank;
    vtkSmartPointer<vtkCPProcessor> Processor;
    vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid;
    vtkSmartPointer<vtkPoints> points;
    vtkSmartPointer<vtkDoubleArray> arrays[4];
    std::vector<real> vx1Array;
    std::vector<real> vx2Array;
    std::vector<real> vx3Array;
    std::vector<real> rhoArray;
    int index;
    int numOfPoints;
    typedef void (*CalcMacrosFct)(const real *const & /*feq[27]*/, real & /*(d)rho*/, real & /*vx1*/,
                                  real & /*vx2*/, real & /*vx3*/);
    CalcMacrosFct calcMacros;
};
#endif // InSituCatalystSimulationObserver_h__

#endif

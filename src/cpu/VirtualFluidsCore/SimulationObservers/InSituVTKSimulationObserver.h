#ifdef VF_VTK

#ifndef InSituVTKSimulationObserver_h__
#define InSituVTKSimulationObserver_h__

#include <SimulationObserver.h>
#include <Grid3D.h>
#include <LBMUnitConverter.h>

#include <string>

// VTK headers
#include <vtkDoubleArray.h>
#include <vtkSmartPointer.h>
#include <vtkSocketCommunicator.h>
#include <vtkSocketController.h>
#include <vtkUnstructuredGrid.h>

class InSituVTKSimulationObserver : public SimulationObserver
{
public:
    InSituVTKSimulationObserver();
    InSituVTKSimulationObserver(SPtr<Grid3D> grid, SPtr<UbScheduler> s, const std::string &configFile,
                         SPtr<LBMUnitConverter> conv);
    virtual ~InSituVTKSimulationObserver();
    void update(real step);

protected:
    void collectData(real step);
    void addData(SPtr<Block3D> block);
    void readConfigFile(const std::string &configFile);

    // void clearData();
private:
    std::string path;
    SPtr<LBMUnitConverter> conv;
    std::vector<std::vector<SPtr<Block3D>>> blockVector;
    int minInitLevel;
    int maxInitLevel;
    int gridRank;
    vtkSmartPointer<vtkSocketCommunicator> comm;
    vtkSmartPointer<vtkSocketController> contr;
    vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid;
    vtkSmartPointer<vtkPoints> points;
    vtkSmartPointer<vtkDoubleArray> arrays[5];
    int wPort;
    std::string wHostname;
    std::string wIP;
};

#endif // InSituVTKSimulationObserver_h__

#endif

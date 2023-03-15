#ifndef D3Q27ShearStressSimulationObserver_H
#define D3Q27ShearStressSimulationObserver_H

#include <PointerDefinitions.h>
#include <string>
#include <vector>

#include <basics/utilities/UbTuple.h>

#include "SimulationObserver.h"

class Block3D;
class Grid3D;
class UbScheduler;
class D3Q27Interactor;
class BCArray3D;
class WbWriter;

//! \brief  Computes the shear stress and y plus values and writes to parallel .vtk
//! \details writes at given time intervals specified in scheduler (s) and resets according to scheduler (rs).
//!          Take root to obtain  during post processing (paraview).
//! \author  K. Kucher, S. Uphoff, M. Geier, E. Goraki Fard

class ShearStressSimulationObserver : public SimulationObserver
{
public:
    //! Default constructor
    ShearStressSimulationObserver() = default;
    //! Constructor
    ShearStressSimulationObserver(SPtr<Grid3D> grid, const std::string &path, WbWriter *const writer, SPtr<UbScheduler> s,
                           SPtr<UbScheduler> rs);
    ~ShearStressSimulationObserver() override;

    void process(real step) override;

    void addInteractor(SPtr<D3Q27Interactor> interactor);

protected:
    //! Computes average and shear stress values of macroscopic quantities
    void calculateShearStress(real timeStep);
    //! Prepare data and write in .vtk file
    void collectData(real step);
    //! Reset data
    void resetData(real step);
    //! prepare data
    void addData();
    void clearData();
    void reset(real step);
    void findPlane(int ix1, int ix2, int ix3, SPtr<Grid3D> grid, SPtr<Block3D> block, real &A, real &B, real &C,
                   real &D, real &ii);
    bool checkUndefindedNodes(SPtr<BCArray3D> bcArray, int ix1, int ix2, int ix3);
    void initDistance();

private:
    std::vector<UbTupleFloat3> nodes;
    std::vector<std::string> datanames;
    std::vector<std::vector<real>> data;
    std::string path;
    std::vector<SPtr<D3Q27Interactor>> interactors;
    std::vector<real> normals;
    int gridRank;
    WbWriter *writer;
    SPtr<UbScheduler> Resetscheduler; // additional scheduler to restart averaging after a given interval
    int minInitLevel;                 // min init level
    int maxInitLevel;
    std::vector<std::vector<SPtr<Block3D>>> blockVector;
    enum Values {
        AvVx          = 0,
        AvVy          = 1,
        AvVz          = 2,
        AvSxx         = 3,
        AvSyy         = 4,
        AvSzz         = 5,
        AvSxy         = 6,
        AvSyz         = 7,
        AvSxz         = 8,
        normalX1      = 9,
        normalX2      = 10,
        normalX3      = 11,
        normalq       = 12,
        numberOfPoint = 13
    };
};

#endif /* D3Q27ShearStressSimulationObserver_H */

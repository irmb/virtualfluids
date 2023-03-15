//! \file QCriterionSimulationObserver.h
//!  \brief Created on: 25.08.2013
//!  \author: Sonja Uphoff

#ifndef QCriterionSimulationObserver_H
#define QCriterionSimulationObserver_H

#include <PointerDefinitions.h>
#include <string>
#include <vector>

#include "SimulationObserver.h"
#include "LBMSystem.h"
#include "UbTuple.h"

namespace vf::mpi {class Communicator;}
class Grid3D;
class UbScheduler;
class WbWriter;
class Block3D;

//! \brief  Computes the value Q with which vortices can be visualized as isocontours to Q=0, writes to .vtk, For
//! uniform, serial setups only! \details writes at given time intervals specified in scheduler (s)
//!          Processing: paraview, take isolines of entry for Q-criterion vortex detection
//!			 Q-Criterion: Visualize Vorteces as regions where Vorticity is larger than strain rate (Hunt, 1988)
//! \author  Sonja Uphoff

class QCriterionSimulationObserver : public SimulationObserver
{
public:
    QCriterionSimulationObserver(SPtr<Grid3D> grid, const std::string &path, WbWriter *const writer, SPtr<UbScheduler> s,
                          std::shared_ptr<vf::mpi::Communicator> comm);
    //! Make update if timestep is write-timestep specified in SPtr<UbScheduler> s
    void process(real step) override;

protected:
    //! Prepare data and write in .vtk file
    void collectData(real step);
    //! Q is computed for all points in a block. Data for writing is added to data and cell vectors.
    void addData(const SPtr<Block3D> block);
    //! After writing to .vtk-file, all vectors are reset
    void clearData();
    //! Computes macroscopic velocities
    void computeVelocity(real *f, real *v, bool compressible);
    //! Computes average and RMS values of macroscopic quantities
    void getNeighborVelocities(int offx, int offy, int offz, int ix1, int ix2, int ix3, const SPtr<Block3D> block,
                               real *vE, real *vW);

private:
    void init();
    std::vector<UbTupleFloat3> nodes;
    std::vector<UbTupleUInt8> cells;
    std::vector<std::string> datanames; // only one entry for QKrit-SimulationObserver: Q
    std::vector<std::vector<real>> data;
    std::vector<std::vector<SPtr<Block3D>>> blockVector;
    int minInitLevel; // go through all levels for block vector of current process from minInitLevel to maxInitLevel
    int maxInitLevel;
    int gridRank; // comm-Rank des aktuellen prozesses
    std::string path;
    WbWriter *writer;
    std::shared_ptr<vf::mpi::Communicator> comm;
    enum Values { xdir = 0, ydir = 1, zdir = 2 }; // labels for the different components
};

#endif

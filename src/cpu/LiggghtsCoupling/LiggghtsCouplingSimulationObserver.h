//=======================================================================================
// ____          ____    __    ______     __________   __      __       __        __
// \    \       |    |  |  |  |   _   \  |___    ___| |  |    |  |     /  \      |  |
//  \    \      |    |  |  |  |  |_)   |     |  |     |  |    |  |    /    \     |  |
//   \    \     |    |  |  |  |   _   /      |  |     |  |    |  |   /  /\  \    |  |
//    \    \    |    |  |  |  |  | \  \      |  |     |   \__/   |  /  ____  \   |  |____
//     \    \   |    |  |__|  |__|  \__\     |__|      \________/  /__/    \__\  |_______|
//      \    \  |    |   ________________________________________________________________
//       \    \ |    |  |  ______________________________________________________________|
//        \    \|    |  |  |         __          __     __     __     ______      _______
//         \         |  |  |_____   |  |        |  |   |  |   |  |   |   _  \    /  _____)
//          \        |  |   _____|  |  |        |  |   |  |   |  |   |  | \  \   \_______
//           \       |  |  |        |  |_____   |   \_/   |   |  |   |  |_/  /    _____  |
//            \ _____|  |__|        |________|   \_______/    |__|   |______/    (_______/
//
//  This file is part of VirtualFluids. VirtualFluids is free software: you can
//  redistribute it and/or modify it under the terms of the GNU General Public
//  License as published by the Free Software Foundation, either version 3 of
//  the License, or (at your option) any later version.
//
//  VirtualFluids is distributed in the hope that it will be useful, but WITHOUT
//  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
//  for more details.
//
//  You should have received a copy of the GNU General Public License along
//  with VirtualFluids (see COPYING.txt). If not, see <http://www.gnu.org/licenses/>.
//
//! \file LiggghtsCouplingSimulationObserver.h
//! \ingroup LiggghtsCoupling
//! \author Konstantin Kutscher
//=======================================================================================

#ifndef LiggghtsCouplingSimulationObserver_h
#define LiggghtsCouplingSimulationObserver_h

#include "SimulationObserver.h"

#include "lammps.h"
#include "input.h"
#include "atom.h"
#include "modify.h"

#include <memory>
#include <vector>


class SimulationObserver;
namespace vf::mpi {class Communicator;}
class LiggghtsCouplingWrapper;
class Grid3D;
class Block3D;
struct IBdynamicsParticleData;
class LBMUnitConverter;

struct ParticleData {
    typedef typename std::vector<std::array<double, 3>> ParticleDataArrayVector;
    typedef typename std::vector<double> ParticleDataScalarVector;
};

class LiggghtsCouplingSimulationObserver : public SimulationObserver
{
public:
    LiggghtsCouplingSimulationObserver(SPtr<Grid3D> grid, SPtr<UbScheduler> s, SPtr<vf::mpi::Communicator> comm,
                                LiggghtsCouplingWrapper &wrapper, int demSteps, SPtr<LBMUnitConverter> units);
    virtual ~LiggghtsCouplingSimulationObserver();

    void process(double actualTimeStep) override;

    
protected:
    void setSpheresOnLattice();
    
    void setSingleSphere3D(double *x, double *v, double *omega, double r, int id /*, bool initVelFlag*/);
    
    double calcSolidFraction(double const dx_, double const dy_, double const dz_, double const r_);
    
    void setValues(IBdynamicsParticleData &p, int const id, double const sf, double const *v, double const dx, double const dy, double const dz, double const *omega);
    
    void setToZero(IBdynamicsParticleData &p);
    
    void getForcesFromLattice();
    
    void SumForceTorque3D(ParticleData::ParticleDataArrayVector &x, double *force, double *torque);

    void addForce(int const partId, int const coord, double const value, double *force);

    void addTorque(int const partId, int const coord, double const value, double *torque);

private:
    SPtr<vf::mpi::Communicator> comm;
    LiggghtsCouplingWrapper &wrapper;
    SPtr<LBMUnitConverter> units;
    int demSteps;
    //std::vector<std::vector<SPtr<Block3D>>> blockVector;
    //int minInitLevel;
    //int maxInitLevel;
    //int gridRank;

    double *force, *torque;
};

#endif


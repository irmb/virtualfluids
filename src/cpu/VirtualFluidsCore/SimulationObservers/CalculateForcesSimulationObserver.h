/*
 *  D3Q27ForcesSimulationObserver.h
 *
 *  Created on: 29.09.2012
 *  Author: K. Kucher
 */

#ifndef D3Q27ForcesSimulationObserver_H
#define D3Q27ForcesSimulationObserver_H

#include <PointerDefinitions.h>
#include <string>
#include <vector>

#include "SimulationObserver.h"
#include "UbTuple.h"
#include "lbm/constants/D3Q27.h"

class ForceCalculator;
namespace vf::mpi {class Communicator;}
class Grid3D;
class UbScheduler;
class D3Q27Interactor;
class DistributionArray3D;
class BoundaryConditions;

class CalculateForcesSimulationObserver : public SimulationObserver
{
public:
    //! Constructor
    //! \param v - velocity of fluid in LB units
    //! \param a - area of object in LB units
    CalculateForcesSimulationObserver(SPtr<Grid3D> grid, SPtr<UbScheduler> s, const std::string &path, std::shared_ptr<vf::mpi::Communicator> comm,
                               real v, real a);
    ~CalculateForcesSimulationObserver() override;
    void update(real step) override;
    void addInteractor(SPtr<D3Q27Interactor> interactor);

protected:
    void collectData(real step);
    void calculateForces();
    UbTupleDouble3 getForces(int x1, int x2, int x3, SPtr<DistributionArray3D> distributions,
                             SPtr<BoundaryConditions> bc);
    void calculateCoefficients();
    void write(std::ofstream *fileObject, real value, char *separator);

private:
    std::string path;
    std::shared_ptr<vf::mpi::Communicator> comm;
    std::vector<SPtr<D3Q27Interactor>> interactors;
    real forceX1global;
    real forceX2global;
    real forceX3global;
    real v; //!< is the speed of the object relative to the fluid
    real a; //!< is the reference area
    real C1;
    real C2;
    real C3;
};

#endif /* D3Q27ForcesSimulationObserver_H */

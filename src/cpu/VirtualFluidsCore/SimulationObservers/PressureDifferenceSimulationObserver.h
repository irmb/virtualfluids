/*
 *  D3Q27PressureDifferenceSimulationObserver.h
 *
 *  Created on: 28.12.2010
 *  Author: kucher
 */

#ifndef D3Q27PRESSUREDIFFERENCESimulationObserver_H
#define D3Q27PRESSUREDIFFERENCESimulationObserver_H

#include <PointerDefinitions.h>
#include <string>

#include "SimulationObserver.h"
#include "LBMSystem.h"

namespace vf::mpi {class Communicator;}
class Grid3D;
class UbScheduler;
class LBMUnitConverter;
class IntegrateValuesHelper;

class PressureDifferenceSimulationObserver : public SimulationObserver
{
public:
    PressureDifferenceSimulationObserver(SPtr<Grid3D> grid, SPtr<UbScheduler> s, const std::string &path,
                                  SPtr<IntegrateValuesHelper> h1, SPtr<IntegrateValuesHelper> h2, real rhoReal,
                                  real uReal, real uLB,
                                  /*const SPtr<LBMUnitConverter> conv,*/ std::shared_ptr<vf::mpi::Communicator> comm);
    ~PressureDifferenceSimulationObserver() override;

    void update(real step) override;

protected:
    SPtr<IntegrateValuesHelper> h1, h2;
    std::string path;
    SPtr<LBMUnitConverter> conv;
    void collectData(real step);
    std::shared_ptr<vf::mpi::Communicator> comm;
    real factor1; //= (1/3)*rhoReal*(uReal/uLB)^2 for calculation pReal = rhoLB * (1/3)*rhoReal*(uReal/uLB)^2,
                     //rhoReal and uReal in SI
    real factor2; //= rhoReal*(uReal/uLB)^2       for calculation pReal = press * rhoReal*(uReal/uLB)^2, rhoReal and
                     //uReal in SI
};

#endif /* D3Q27RHODIFFERENCESimulationObserver_H_ */

/*
 * D3Q27RhoSimulationObserver.cpp
 *
 *  Created on: 28.12.2010
 *      Author: kucher
 */

#include "PressureDifferenceSimulationObserver.h"

#include <fstream>

#include <mpi/Communicator.h>
#include "Grid3D.h"
#include "IntegrateValuesHelper.h"
#include "LBMUnitConverter.h"
#include "UbScheduler.h"

PressureDifferenceSimulationObserver::PressureDifferenceSimulationObserver(SPtr<Grid3D> grid, SPtr<UbScheduler> s,
                                                             const std::string &path, SPtr<IntegrateValuesHelper> h1,
                                                             SPtr<IntegrateValuesHelper> h2, real rhoReal,
                                                             real uReal, real uLB, std::shared_ptr<vf::mpi::Communicator> comm)

    : SimulationObserver(grid, s), path(path), h1(h1), h2(h2), comm(comm)
{
    if (comm->getProcessID() == comm->getRoot()) {
        std::ofstream ostr;
        std::string fname = path;
        ostr.open(fname.c_str(), std::ios_base::out);
        if (!ostr) {
            ostr.clear();
            std::string file_path = UbSystem::getPathFromString(fname);
            if (file_path.size() > 0) {
                UbSystem::makeDirectory(file_path);
                ostr.open(fname.c_str(), std::ios_base::out);
            }
            if (!ostr)
                throw UbException(UB_EXARGS, "couldn't open file " + fname);
        }
        ostr << "step"
             << "\t"
             << "nodes1"
             << "\t"
             << "nodes2"
             << "\t";
        ostr << "sRho1"
             << "\t"
             << "p1_1"
             << "\t"
             << "sRho2"
             << "\t"
             << "p1_2"
             << "\t"
             << "deltaP1"
             << "\t";
        ostr << "sPress1"
             << "\t"
             << "p2_1"
             << "\t"
             << "sPress2"
             << "\t"
             << "p2_2"
             << "\t"
             << "deltaP2";
        ostr << std::endl;
        ostr.close();

        factor1 = (1.0 / 3.0) * rhoReal * (uReal / uLB) * (uReal / uLB);
        factor2 = rhoReal * (uReal / uLB) * (uReal / uLB);
    }
}
//////////////////////////////////////////////////////////////////////////
PressureDifferenceSimulationObserver::~PressureDifferenceSimulationObserver() = default;
//////////////////////////////////////////////////////////////////////////
void PressureDifferenceSimulationObserver::update(real step)
{
    if (scheduler->isDue(step))
        collectData(step);
}
//////////////////////////////////////////////////////////////////////////
void PressureDifferenceSimulationObserver::collectData(real step)
{
    h1->calculateMQ();
    h2->calculateMQ();

    if (comm->getProcessID() == comm->getRoot()) {
        int istep = static_cast<int>(step);
        std::ofstream ostr;
        real nn1  = h1->getNumberOfFluidsNodes();
        real nn2  = h2->getNumberOfFluidsNodes();
        real rho1 = h1->getRho();
        real rho2 = h2->getRho();
        real p1_1 = (rho1 / nn1) * factor1;
        real p1_2 = (rho2 / nn2) * factor1;
        real dp1  = p1_1 - p1_2;

        // double press1 = h1->getPress();
        // double press2 = h2->getPress();
        // double p2_1 = (press1/nn1) * factor2;
        // double p2_2 = (press2/nn2) * factor2;
        // double dp2 = p2_1 - p2_2;

        std::string fname = path;
        ostr.open(fname.c_str(), std::ios_base::out | std::ios_base::app);
        if (!ostr) {
            ostr.clear();
            std::string path = UbSystem::getPathFromString(fname);
            if (path.size() > 0) {
                UbSystem::makeDirectory(path);
                ostr.open(fname.c_str(), std::ios_base::out | std::ios_base::app);
            }
            if (!ostr)
                throw UbException(UB_EXARGS, "couldn't open file " + fname);
        }

        ostr << istep << "\t" << nn1 << "\t" << nn2 << "\t";
        ostr << rho1 << "\t" << p1_1 << "\t" << rho2 << "\t" << p1_2 << "\t" << dp1 << "\t";
        // ostr << press1 << "\t" << p2_1 << "\t" << press2 << "\t" << p2_2 << "\t" << dp2;
        ostr << std::endl;
        ostr.close();
    }
}

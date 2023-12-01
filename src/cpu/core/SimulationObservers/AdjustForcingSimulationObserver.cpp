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
//! \file AdjustForcingSimulationObserver.cpp
//! \ingroup SimulationObservers
//! \author Konstantin Kutscher
//=======================================================================================
#include "AdjustForcingSimulationObserver.h"

#include <fstream>

#include <parallel/Communicator.h>
#include "Grid3D.h"
#include "IntegrateValuesHelper.h"
#include "UbScheduler.h"
#include <SetForcingBlockVisitor.h>

AdjustForcingSimulationObserver::AdjustForcingSimulationObserver(SPtr<Grid3D> grid, SPtr<UbScheduler> s, const std::string &path,
                                                   SPtr<IntegrateValuesHelper> integrateValues, real vTarged,
                                                   std::shared_ptr<vf::parallel::Communicator> comm)

    : SimulationObserver(grid, s), path(path), integrateValues(integrateValues), comm(comm), vx1Targed(vTarged)
{
    using namespace  vf::basics::constant;

    // cnodes = integrateValues->getCNodes();
    root = comm->isRoot();

    Ta = scheduler->getMaxStep();

    Kpcrit = c3o1 / Ta; // 0.3;
    Tcrit  = c3o1 * Ta; // 30.0;
    Tn     = c1o2 * Tcrit;
    Tv     = c12o1 / c100o1 * Tcrit;

    Kp = c6o1 / c10o1 * Kpcrit;
    Ki = Kp / Tn;
    Kd = Kp * Tv;

    y       = c0o1;
    e       = c0o1;
    esum    = c0o1;
    eold    = c0o1;
    forcing = c0o1;

    if (root) {
        std::string fname = path + "/forcing/forcing.csv";
        std::ofstream ostr;
        ostr.open(fname.c_str(), std::ios_base::out | std::ios_base::app);
        if (!ostr) {
            ostr.clear();
            std::string file_path = UbSystem::getPathFromString(fname);
            if (file_path.size() > 0) {
                UbSystem::makeDirectory(file_path);
                ostr.open(fname.c_str(), std::ios_base::out | std::ios_base::app);
            }
            if (!ostr)
                throw UbException(UB_EXARGS, "couldn't open file " + fname);
        }
        ostr << "step;volume;vx1average;forcing\n";
        ostr.close();

        //////////////////////////////////////////////////////////////////////////////////////////////////
        // temporary solution
        std::string fNameCfg = path + "/forcing/forcing.cfg";
        std::ifstream istr2;
        istr2.open(fNameCfg.c_str(), std::ios_base::in);
        if (istr2) {
            istr2 >> forcing;
            // istr2 >> esum;
            // istr2 >> eold;
        }
        istr2.close();
    }
    ////////////////////////////////////////////////////////////////////////////////////////////////////////
}

//////////////////////////////////////////////////////////////////////////
void AdjustForcingSimulationObserver::update(real step)
{
    if (scheduler->isDue(step))
        collectData(step);
}
//////////////////////////////////////////////////////////////////////////
void AdjustForcingSimulationObserver::collectData(real step)
{
    //////////////////////////////////////////////////////////////////////////////////////////////////
    // temporary solution
    if (root) {
        std::string fNameCfg = path + "/forcing/forcing.cfg";
        std::ofstream ostr2;
        ostr2.open(fNameCfg.c_str(), std::ios_base::out);
        if (!ostr2) {
            ostr2.clear();
            std::string path = UbSystem::getPathFromString(fNameCfg);
            if (path.size() > 0) {
                UbSystem::makeDirectory(path);
                ostr2.open(fNameCfg.c_str(), std::ios_base::out);
            }
            if (!ostr2)
                throw UbException(UB_EXARGS, "couldn't open file " + fNameCfg);
        }
        ostr2 << forcing << " " << esum << " " << eold;
        ostr2.close();
    }
    ////////////////////////////////////////////////////////////////////////////////////////////////////////

    integrateValues->calculateMQ();

    if (root) {
        cellsVolume = integrateValues->getCellsVolume();
        real vx1  = integrateValues->getVx1();
        vx1Average  = (vx1 / cellsVolume);

        //////////////////////////////////////////////////////////////////////////
        // PID-Controller (PID-Regler)
        e    = vx1Targed - vx1Average;
        esum = esum + e;
        y    = Kp * e + Ki * Ta * esum + Kd * (e - eold) / Ta;
        eold = e;

        forcing = forcing + y;
        //////////////////////////////////////////////////////////////////////////
    }
    //////////////////////////////////////////////////////////////////////////
    comm->broadcast(forcing);

    mu::Parser fctForcingX1, fctForcingX2, fctForcingX3;
    fctForcingX1.SetExpr("Fx1");
    fctForcingX1.DefineConst("Fx1", forcing);
    fctForcingX2.SetExpr("0.0");
    fctForcingX3.SetExpr("0.0");
    SetForcingBlockVisitor forcingVisitor(fctForcingX1, fctForcingX2, fctForcingX3);
    grid->accept(forcingVisitor);

    // for(CalcNodes cn : cnodes)
    //{
    //   LBMKernel3DPtr kernel = cn.block->getKernel();
    //   if (kernel)
    //   {
    //      kernel->setForcingX1(fctForcingX1);
    //      kernel->setWithForcing(true);
    //   }
    //
    //}

    if (root) {
        // UBLOG(logINFO, "D3Q27AdjustForcingSimulationObserver step: " << static_cast<int>(step));
        // UBLOG(logINFO, "new forcing is: " << forcing);
        std::string fname = path + "/forcing/forcing.csv";
        // std::string fname = path + "/forcing/forcing_"+UbSystem::toString(comm->getProcessID())+".csv";
        std::ofstream ostr;
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
        int istep = static_cast<int>(step);

        // ostr << istep << ";" << cellsVolume << ";" << vx1Average << "; " << forcing << "\n";
        ostr << istep << ";" << cellsVolume << ";" << vx1Average << "; " << forcing << "; " << e << "; " << esum << "; "
             << y << "\n";
        ostr.close();
    }
}

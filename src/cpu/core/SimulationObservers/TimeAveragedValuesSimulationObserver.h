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
//  FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
//  for more details.
//
//  SPDX-License-Identifier: GPL-3.0-or-later
//  SPDX-FileCopyrightText: Copyright © VirtualFluids Project contributors, see AUTHORS.md in root folder
//
//! \addtogroup cpu_SimulationObservers SimulationObservers
//! \ingroup cpu_core core
//! \{
//! \author Konstantin Kutscher
//=======================================================================================
#ifndef TimeAveragedValuesSimulationObserver_H
#define TimeAveragedValuesSimulationObserver_H

#include <PointerDefinitions.h>
#include <string>
#include <vector>

#include "SimulationObserver.h"
#include "IntegrateValuesHelper.h"
#include "LBMSystem.h"

namespace vf::parallel {class Communicator;}
class Grid3D;
class UbScheduler;
class WbWriter;
class Block3D;

//! \brief  Computes the time averaged mean velocity and RMS values and writes to parallel .vtk
//! \details writes at given time intervals specified in scheduler (s), does averaging according to scheduler (Avs) and
//! resets according to scheduler (rs).  <br>
//!  Computes  the time averaged mean velocity  \f$ u_{mean}=\frac{1}{N}\sum\limits_{i=1}^n u_{i} \f$  and RMS of
//!  fluctuations. You need to calculate a square root before plotting RMS. <br>
//
//! \author  Konstantin Kutscher
// \f$ u_{mean}=\frac{1}{N}\sum\limits_{i=1}^n u_{i} \f$

class TimeAveragedValuesSimulationObserver : public SimulationObserver
{
public:
    enum Options {
        Density            = 1,
        Velocity           = 2,
        Fluctuations       = 4,
        Triplecorrelations = 8,

        // Velocity           = 1,
        // Fluctuations       = 2,
        // Triplecorrelations = 4,
    };

public:
    TimeAveragedValuesSimulationObserver();
    TimeAveragedValuesSimulationObserver(SPtr<Grid3D> grid, const std::string &path, WbWriter *const writer,
                                  SPtr<UbScheduler> s, std::shared_ptr<vf::parallel::Communicator> comm, int options);
    TimeAveragedValuesSimulationObserver(SPtr<Grid3D> grid, const std::string &path, WbWriter *const writer,
                                  SPtr<UbScheduler> s, std::shared_ptr<vf::parallel::Communicator> comm, int options, std::vector<int> levels,
                                  std::vector<real> &levelCoords, std::vector<real> &bounds,
                                  bool timeAveraging = true);
    //! Make update
    void update(real step) override;
    //! Computes subtotal of velocity , fluctuations and triple correlations
    void calculateSubtotal(real step);
    void addLevelCoordinate(real c);
    void reset();
    void setWithGhostLayer(bool val);
    bool getWithGhostLayer();

protected:
    //! Prepare data and write in .vtk file
    void collectData(real step);
    //! prepare data
    void addData(const SPtr<Block3D> block);
    void clearData();
    //! Computes average values of velocity , fluctuations and triple correlations
    void calculateAverageValues(real timeStep);

    void init();
    void initData();
    void planarAverage(real step);
    void calculateAverageValuesForPlane(std::vector<IntegrateValuesHelper::CalcNodes> &cnodes);

private:
    std::shared_ptr<vf::parallel::Communicator> comm;
    std::vector<UbTupleFloat3> nodes;
    std::vector<UbTupleUInt8> cells;
    std::vector<std::string> datanames;
    std::vector<std::vector<double>> data;
    std::vector<std::vector<SPtr<Block3D>>> blockVector;
    bool root;
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
    enum Density { Rho, RhoF };
    enum Velocity { Vx, Vy, Vz };
    enum Fluctuations { Vxx, Vyy, Vzz, Vxy, Vxz, Vyz };
    enum Triplecorrelations { Vxxx, Vxxy, Vxxz, Vyyy, Vyyx, Vyyz, Vzzz, Vzzx, Vzzy, Vxyz };

    real saRho, saRhoF;
    real saVx, saVy, saVz;
    real saVxx, saVyy, saVzz, saVxy, saVxz, saVyz;
    real saVxxx, saVxxy, saVxxz, saVyyy, saVyyx, saVyyz, saVzzz, saVzzx, saVzzy, saVxyz;

    int options;
    real numberOfSteps;
    real minStep;
    real maxStep;

    int iMinX1, iMinX2, iMinX3;
    // int iMaxX1, iMaxX2, iMaxX3;
    int iMinC;
    int iMaxC;

    using CalcMacrosFct = void (*)(const real *const &, real &, real &, real &, real &);
    CalcMacrosFct calcMacros;

    bool planarAveraging;
    bool timeAveraging;
    std::vector<real> levelCoords;
    std::vector<int> levels;
    std::vector<real> bounds;

    bool withGhostLayer;
};
#endif

//! \}

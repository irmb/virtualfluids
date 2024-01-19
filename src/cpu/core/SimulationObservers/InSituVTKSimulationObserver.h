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

//! \}

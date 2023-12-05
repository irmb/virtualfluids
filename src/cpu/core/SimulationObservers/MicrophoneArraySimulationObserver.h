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
//! \file MicrophoneArraySimulationObserver.h
//! \ingroup SimulationObservers
//! \author Konstantin Kutscher
//=======================================================================================
#ifndef MicrophoneArraySimulationObserver_h__
#define MicrophoneArraySimulationObserver_h__

#include "SimulationObserver.h"
#include "LBMSystem.h"
#include "UbTuple.h"
#include <array>
#include <string>
#include <vector>

namespace vf::parallel {class Communicator;}
class Grid3D;
class UbScheduler;
class Vector3D;
class DistributionArray3D;

//! \brief     Class implements microphone array.
//! \details   It samples pressure (LBM rho value) and stores to .csv file.
//! \author    Konstantin Kutscher
//! \date      February 2019

class MicrophoneArraySimulationObserver : public SimulationObserver
{
public:
    MicrophoneArraySimulationObserver(SPtr<Grid3D> grid, SPtr<UbScheduler> s, const std::string &path,
                               std::shared_ptr<vf::parallel::Communicator> comm);
    ~MicrophoneArraySimulationObserver() override;

    //! calls collectData.
    void update(real step) override;

    //! add microphone
    bool addMicrophone(Vector3D coords);

protected:
    void collectData(real step);
    void writeFile(real step);

private:
    std::string path;
    std::shared_ptr<vf::parallel::Communicator> comm;

    struct Mic {
        unsigned int id;
        SPtr<DistributionArray3D> distridution;
        UbTupleInt3 nodeIndexes;
    };
    std::vector<SPtr<Mic>> microphones;

    std::vector<SPtr<std::stringstream>> strVector;

    int count;
    int micID;

    using CalcMacrosFct = void (*)(const real *const &, real &, real &, real &, real &);
    CalcMacrosFct calcMacros;
};

#endif // MicrophoneArraySimulationObserver_h__

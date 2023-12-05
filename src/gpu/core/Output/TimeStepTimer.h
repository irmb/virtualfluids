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
//! \author Martin Schoenherr
//=======================================================================================
#ifndef TIMESTEPTIMER_H
#define TIMESTEPTIMER_H

#include <string>

#include <basics/DataTypes.h>

#include <basics/Timer/Timer.h>

class TimeStepTimer
{
public:
    TimeStepTimer(std::string _name, uint _tActivate) : name(_name), tActivate(_tActivate) {};

    void startTotalTimer            (uint t);
    void stopTotalTimer             (uint t);
    void startCollisionTimer        (uint t);
    void stopCollisionTimer         (uint t);
    void startPostCollisionBCTimer  (uint t);
    void stopPostCollisionBCTimer   (uint t);
    void startPreCollisionBCTimer   (uint t);
    void stopPreCollisionBCTimer    (uint t);
    void startEddyViscosityTimer    (uint t);
    void stopEddyViscosityTimer     (uint t);
    void startActuatorTimer         (uint t);
    void stopActuatorTimer          (uint t);
    void startProbeTimer            (uint t);
    void stopProbeTimer             (uint t);
    void startExchangeTimer         (uint t);
    void stopExchangeTimer          (uint t);

    void outputPerformance(uint t);

private:
    vf::basics::Timer totalTimer { "total" };
    vf::basics::Timer collisionTimer { "collision" };
    vf::basics::Timer postCollisionBCTimer { "postCollisionBC" };
    vf::basics::Timer preCollisionBCTimer { "preCollisionBC" };
    vf::basics::Timer eddyViscosityTimer { "eddyViscosity" };
    vf::basics::Timer actuatorTimer { "actuator" };
    vf::basics::Timer probeTimer { "probes" };
    vf::basics::Timer exchangeTimer { "exchange" };

    std::string name;
    uint tActivate;
};

#endif 

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
//! \addtogroup gpu_Output Output
//! \ingroup gpu_core core
//! \{
//! \author Martin Schoenherr
//=======================================================================================
#include "TimeStepTimer.h"

#include <basics/DataTypes.h>

#include <logger/Logger.h>

void TimeStepTimer::startTotalTimer             (uint t){ if(t%this->tActivate==0) this->totalTimer.start();              }
void TimeStepTimer::stopTotalTimer              (uint t){ if(t%this->tActivate==0) this->totalTimer.end();               }
void TimeStepTimer::startCollisionTimer         (uint t){ if(t%this->tActivate==0) this->collisionTimer.start();          }
void TimeStepTimer::stopCollisionTimer          (uint t){ if(t%this->tActivate==0) this->collisionTimer.end();           }
void TimeStepTimer::startPostCollisionBCTimer   (uint t){ if(t%this->tActivate==0) this->postCollisionBCTimer.start();    }
void TimeStepTimer::stopPostCollisionBCTimer    (uint t){ if(t%this->tActivate==0) this->postCollisionBCTimer.end();     }
void TimeStepTimer::startPreCollisionBCTimer    (uint t){ if(t%this->tActivate==0) this->preCollisionBCTimer.start();     }
void TimeStepTimer::stopPreCollisionBCTimer     (uint t){ if(t%this->tActivate==0) this->preCollisionBCTimer.end();      }
void TimeStepTimer::startEddyViscosityTimer     (uint t){ if(t%this->tActivate==0) this->eddyViscosityTimer.start();      }
void TimeStepTimer::stopEddyViscosityTimer      (uint t){ if(t%this->tActivate==0) this->eddyViscosityTimer.end();       }
void TimeStepTimer::startActuatorTimer          (uint t){ if(t%this->tActivate==0) this->actuatorTimer.start();           }
void TimeStepTimer::stopActuatorTimer           (uint t){ if(t%this->tActivate==0) this->actuatorTimer.end();            }
void TimeStepTimer::startProbeTimer             (uint t){ if(t%this->tActivate==0) this->probeTimer.start();              }
void TimeStepTimer::stopProbeTimer              (uint t){ if(t%this->tActivate==0) this->probeTimer.end();               }
void TimeStepTimer::startExchangeTimer          (uint t){ if(t%this->tActivate==0) this->exchangeTimer.start();           }
void TimeStepTimer::stopExchangeTimer           (uint t){ if(t%this->tActivate==0) this->exchangeTimer.end();            }

void TimeStepTimer::outputPerformance(uint t)
{
    if (t % this->tActivate == 0) {
        float tCollision         = this->collisionTimer.getTimeInSeconds();
        float tPostCollisionBC   = this->postCollisionBCTimer.getTimeInSeconds();
        float tPreCollisionBC    = this->preCollisionBCTimer.getTimeInSeconds();
        float tEddyViscosity     = this->eddyViscosityTimer.getTimeInSeconds();
        float tAcutator          = this->actuatorTimer.getTimeInSeconds();
        float tProbe             = this->probeTimer.getTimeInSeconds();
        float tExchange          = this->exchangeTimer.getTimeInSeconds();
        float tTotal             = tCollision+tPostCollisionBC+tPreCollisionBC+tEddyViscosity+tAcutator+tProbe+tExchange;
        
        VF_LOG_INFO(" --- Collision \t {}%",        (tCollision/tTotal)*100 );
        VF_LOG_INFO(" --- PostCollisionBCs \t {}%", (tPostCollisionBC/tTotal)*100 );
        VF_LOG_INFO(" --- PreCollisionBCs \t {}%",  (tPreCollisionBC/tTotal)*100 );
        VF_LOG_INFO(" --- Eddy viscosity \t {}%",   (tEddyViscosity/tTotal)*100 );
        VF_LOG_INFO(" --- Actuators \t {}%",        (tAcutator/tTotal)*100 );
        VF_LOG_INFO(" --- Probes \t\t {}%",           (tProbe/tTotal)*100 );
        VF_LOG_INFO(" --- Data exchange \t {}%",    (tExchange/tTotal)*100 );
    }
}

//! \}

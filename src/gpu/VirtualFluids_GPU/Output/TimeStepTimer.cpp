#include "helper_cuda.h"
#include <cuda_runtime.h>
#include "Core/DataTypes.h"
#include "UbScheduler.h"
#include "Parameter/Parameter.h"

#include "Timer.h"
#include "TimeStepTimer.h"

void TimeStepTimer::startTotalTimer             (uint t){ if(t%this->tActivate==0) this->totalTimer->startTimer();              }
void TimeStepTimer::stopTotalTimer              (uint t){ if(t%this->tActivate==0) this->totalTimer->stopTimer();               }
void TimeStepTimer::startCollisionTimer         (uint t){ if(t%this->tActivate==0) this->collisionTimer->startTimer();          }
void TimeStepTimer::stopCollisionTimer          (uint t){ if(t%this->tActivate==0) this->collisionTimer->stopTimer();           }
void TimeStepTimer::startPostCollisionBCTimer   (uint t){ if(t%this->tActivate==0) this->postCollisionBCTimer->startTimer();    }
void TimeStepTimer::stopPostCollisionBCTimer    (uint t){ if(t%this->tActivate==0) this->postCollisionBCTimer->stopTimer();     }
void TimeStepTimer::startPreCollisionBCTimer    (uint t){ if(t%this->tActivate==0) this->preCollisionBCTimer->startTimer();     }
void TimeStepTimer::stopPreCollisionBCTimer     (uint t){ if(t%this->tActivate==0) this->preCollisionBCTimer->stopTimer();      }
void TimeStepTimer::startEddyViscosityTimer     (uint t){ if(t%this->tActivate==0) this->eddyViscosityTimer->startTimer();      }
void TimeStepTimer::stopEddyViscosityTimer      (uint t){ if(t%this->tActivate==0) this->eddyViscosityTimer->stopTimer();       }
void TimeStepTimer::startActuatorTimer          (uint t){ if(t%this->tActivate==0) this->actuatorTimer->startTimer();           }
void TimeStepTimer::stopActuatorTimer           (uint t){ if(t%this->tActivate==0) this->actuatorTimer->stopTimer();            }
void TimeStepTimer::startProbeTimer             (uint t){ if(t%this->tActivate==0) this->probeTimer->startTimer();              }
void TimeStepTimer::stopProbeTimer              (uint t){ if(t%this->tActivate==0) this->probeTimer->stopTimer();               }
void TimeStepTimer::startExchangeTimer          (uint t){ if(t%this->tActivate==0) this->exchangeTimer->startTimer();           }
void TimeStepTimer::stopExchangeTimer           (uint t){ if(t%this->tActivate==0) this->exchangeTimer->stopTimer();            }


void TimeStepTimer::resetTimers(uint t)
{
    if(t%this->tActivate==0)
    {
        this->totalTimer->resetTimer();
        this->collisionTimer->resetTimer();
        this->postCollisionBCTimer->resetTimer();
        this->preCollisionBCTimer->resetTimer();
        this->eddyViscosityTimer->resetTimer();
        this->actuatorTimer->resetTimer();
        this->probeTimer->resetTimer();
    }
}

void TimeStepTimer::outputPerformance(uint t, Parameter* para)
{
    if(t%this->tActivate==0)
    {
        
        real tCollision         = this->collisionTimer->getTotalElapsedTime();
        real tPostCollisionBC   = this->postCollisionBCTimer->getTotalElapsedTime();
        real tPreCollisionBC    = this->preCollisionBCTimer->getTotalElapsedTime();
        real tEddyViscosity     = this->eddyViscosityTimer->getTotalElapsedTime();
        real tAcutator          = this->actuatorTimer->getTotalElapsedTime();
        real tProbe             = this->probeTimer->getTotalElapsedTime();
        real tExchange          = this->exchangeTimer->getTotalElapsedTime();
        real tTotal             = tCollision+tPostCollisionBC+tPreCollisionBC+tEddyViscosity+tAcutator+tProbe+tExchange;
        
        VF_LOG_INFO(" --- Collision \t {}%",        (tCollision/tTotal)*100 );
        VF_LOG_INFO(" --- PostCollisionBCs \t {}%", (tPostCollisionBC/tTotal)*100 );
        VF_LOG_INFO(" --- PreCollisionBCs \t {}%",  (tPreCollisionBC/tTotal)*100 );
        VF_LOG_INFO(" --- Eddy viscosity \t {}%",   (tEddyViscosity/tTotal)*100 );
        VF_LOG_INFO(" --- Actuators \t {}%",        (tAcutator/tTotal)*100 );
        VF_LOG_INFO(" --- Probes \t\t {}%",           (tProbe/tTotal)*100 );
        VF_LOG_INFO(" --- Data exchange \t {}%",    (tExchange/tTotal)*100 );
    }
}
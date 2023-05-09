#ifndef TIMESTEPTIMER_H
#define TIMESTEPTIMER_H

#include "helper_cuda.h"
#include <cuda_runtime.h>
#include "DataTypes.h"
#include "UbScheduler.h"
#include "Parameter/Parameter.h"

#include "Timer.h"

class TimeStepTimer
{
    public:
    TimeStepTimer(std::string _name, uint _tActivate): name(_name), tActivate(_tActivate)
    {
        
    };
    
    ~TimeStepTimer(){};

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

    void resetTimers(uint t);
    void outputPerformance(uint t, Parameter* para);

    private:
    
    Timer* totalTimer           = new Timer("total");
    Timer* collisionTimer       = new Timer("collision");
    Timer* postCollisionBCTimer = new Timer("postCollisionBC");
    Timer* preCollisionBCTimer  = new Timer("preCollisionBC");
    Timer* eddyViscosityTimer   = new Timer("eddyViscosity");
    Timer* actuatorTimer        = new Timer("actuator");
    Timer* probeTimer           = new Timer("probes");
    Timer* exchangeTimer        = new Timer("exchange");
    
    std::string name;
    uint tActivate;
};



#endif 
#ifndef SIMULATION_OBSERVER_H
#define SIMULATION_OBSERVER_H

class SimulationObserver
{
public:
    virtual ~SimulationObserver() = default;
    virtual void update() = 0;

};
#endif 
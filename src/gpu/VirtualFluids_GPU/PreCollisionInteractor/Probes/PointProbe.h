#ifndef PointProbe_H
#define PointProbe_H

#include "Probe.h"

class PointProbe: public Probe
{
public:
    PointProbe(
        const std::string _probeName,
        uint _tStartAvg,
        uint _tStartOut,
        uint _tOut
    ): Probe(_probeName, 
             _tStartAvg, 
             _tStartOut, 
             _tOut)
    {}

    void addProbePointsFromList(std::vector<real>& _pointCoordsX, std::vector<real>& _pointCoordsY, std::vector<real>& _pointCoordsZ);
    void addProbePointsFromXNormalPlane(real pos_x, real pos0_y, real pos0_z, real pos1_y, real pos1_z, real delta_y, real delta_z);
    
private:
    void findPoints(Parameter* para, GridProvider* gridProvider, std::vector<int>& probeIndices_level,
                    std::vector<real>& distX_level, std::vector<real>& distY_level, std::vector<real>& distZ_level,      
                    std::vector<real>& pointCoordsX_level, std::vector<real>& pointCoordsY_level, std::vector<real>& pointCoordsZ_level,
                    int level) override;

    void calculateQuantities(SPtr<ProbeStruct> probeStruct, Parameter* para, int level) override;

private:
    uint nProbePoints;
    std::vector<real> pointCoordsX, pointCoordsY, pointCoordsZ; 

};

#endif
#ifndef PlaneProbe_H
#define PlaneProbe_H

#include "Probe.h"

class PlaneProbe : public Probe
{
public: 
    PlaneProbe(
        const std::string _probeName,
        const std::string _outputPath,
        uint _tStartAvg,
        uint _tAvg,
        uint _tStartOut,
        uint _tOut
    ): Probe(_probeName, 
             _outputPath,
             _tStartAvg, 
             _tAvg,
             _tStartOut, 
             _tOut,
             true)
    {}

    void setProbePlane(real _posX, real _posY, real _posZ, real _deltaX, real _deltaY, real _deltaZ)
    {
        this->posX = _posX; 
        this->posY = _posY; 
        this->posZ = _posZ;         
        this->deltaX = _deltaX; 
        this->deltaY = _deltaY; 
        this->deltaZ = _deltaZ; 
    }

private:
    bool isAvailablePostProcessingVariable(PostProcessingVariable _variable) override;

    void findPoints(Parameter* para, GridProvider* gridProvider, std::vector<int>& probeIndices_level,
                    std::vector<real>& distX_level, std::vector<real>& distY_level, std::vector<real>& distZ_level,      
                    std::vector<real>& pointCoordsX_level, std::vector<real>& pointCoordsY_level, std::vector<real>& pointCoordsZ_level,
                    int level) override;
    void calculateQuantities(SPtr<ProbeStruct> probeStruct, Parameter* para, int level) override;

private:
    real posX, posY, posZ;
    real deltaX, deltaY, deltaZ;
};

#endif
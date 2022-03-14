#ifndef PlanarAverageProbe_H
#define PlanarAverageProbe_H

#include "Probe.h"

__global__ void moveIndicesInNegNormalDir( uint* pointIndices, uint nPoints, uint* neighborWSB, uint* neighborInplane1, uint* neighborInplane2, real* coordsX, real* coordsY, real* coordsZ ); 

__global__ void moveIndicesInPosNormalDir( uint* pointIndices, uint nPoints, uint* neighborNormal, real* coordsX, real* coordsY, real* coordsZ );

///////////////////////////////////////////////////////////////////////////////////

class PlanarAverageProbe : public Probe
{
public: 
    PlanarAverageProbe(
        const std::string _probeName,
        const std::string _outputPath,
        uint _tStartAvg,
        uint _tAvg,
        uint _tStartOut,
        uint _tOut,
        char _planeNormal
    ):  Probe(_probeName, 
             _outputPath,
             _tStartAvg, 
             _tAvg,
             _tStartOut, 
             _tOut,
             false),
        planeNormal(_planeNormal)
    {
        assert(_planeNormal == 'x' || _planeNormal == 'y' || _planeNormal == 'z');
        // this->hasDeviceQuantityArray = false;
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
    char planeNormal;
    bool isEvenTAvg = true;
};

#endif
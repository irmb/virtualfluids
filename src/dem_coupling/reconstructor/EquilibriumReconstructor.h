/*
*  Author: S. Peters
*  mail: peters@irmb.tu-bs.de
*/
#ifndef EQUILIBRIUM_RECONSTRUCTOR_H
#define EQUILIBRIUM_RECONSTRUCTOR_H

#include "boost/shared_ptr.hpp"
#include "UbTuple.h"

#include "Reconstructor.h"

class ILBMKernel;
class PhysicsEngineGeometryAdapter;

class EquilibriumReconstructor : public Reconstructor
{
public:
    virtual ~EquilibriumReconstructor() {}

    void reconstructNode(const int &x1, const int &x2, const int &x3, const Vector3D& worldCoordinates, std::shared_ptr<PhysicsEngineGeometryAdapter> physicsEngineGeometry, std::shared_ptr<ILBMKernel> kernel) const override;

private:
    double getLocalAverageDensity(const int &x1, const int &x2, const int &x3, std::shared_ptr<ILBMKernel> kernel) const;
};



#endif


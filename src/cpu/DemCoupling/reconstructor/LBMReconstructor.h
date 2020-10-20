/*
 *  Author: S. Peters
 *  mail: peters@irmb.tu-bs.de
 */
#ifndef LBM_RECONSTRUCTOR_H
#define LBM_RECONSTRUCTOR_H

#include "UbTuple.h"

#include "Reconstructor.h"

#include "LBMSystem.h"

class ILBMKernel;
class PhysicsEngineGeometryAdapter;

class LBMReconstructor : public Reconstructor
{
public:
    LBMReconstructor(bool compressible);
    virtual ~LBMReconstructor() {}

    void reconstructNode(const int &x1, const int &x2, const int &x3, const Vector3D &worldCoordinates,
                         std::shared_ptr<PhysicsEngineGeometryAdapter> physicsEngineGeometry,
                         std::shared_ptr<ILBMKernel> kernel) const override;

private:
    static void collide(LBMReal *f, LBMReal collFactor);

    typedef void (*CalcMacrosFct)(const LBMReal *const & /*f[27]*/, LBMReal & /*rho*/, LBMReal & /*vx1*/,
                                  LBMReal & /*vx2*/, LBMReal & /*vx3*/);
    CalcMacrosFct calcMacrosFct;
};

#endif

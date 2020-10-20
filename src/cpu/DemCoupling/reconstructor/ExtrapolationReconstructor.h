/*
 *  Author: S. Peters
 *  mail: peters@irmb.tu-bs.de
 */
#ifndef EXTRAPOLATION_RECONSTRUCTOR_H
#define EXTRAPOLATION_RECONSTRUCTOR_H

#include <memory>

#include "UbTuple.h"

#include "Reconstructor.h"

class ILBMKernel;
class PhysicsEngineGeometryAdapter;

class ExtrapolationReconstructor : public Reconstructor
{
public:
    ExtrapolationReconstructor(std::shared_ptr<Reconstructor> alternativeReconstructor);
    virtual ~ExtrapolationReconstructor() {}

    void reconstructNode(const int &x1, const int &x2, const int &x3, const Vector3D &worldCoordinates,
                         std::shared_ptr<PhysicsEngineGeometryAdapter> physicsEngineGeometry,
                         std::shared_ptr<ILBMKernel> kernel) const override;

    void setAlternativeReconstructor(std::shared_ptr<Reconstructor> alternativeReconstructor);

private:
    int getNumberOfExtrapolationCells(const int x1, const int x2, const int x3, const UbTupleInt3 &ubTuple,
                                      std::shared_ptr<ILBMKernel> kernel) const;
    UbTupleInt3 getSphereDirection(const Vector3D &worldCoordinates,
                                   std::shared_ptr<PhysicsEngineGeometryAdapter> physicsEngineGeometry) const;
    UbTupleInt3 getCorrespondingLatticeDirection(const Vector3D &direction) const;
    void extrapolatePdFs(const int x1, const int x2, const int x3, const UbTupleInt3 &ubTuple,
                         int numberOfCellsForExtrapolation, std::shared_ptr<ILBMKernel> kernel) const;

    std::shared_ptr<Reconstructor> alternativeReconstructor;
};

#endif

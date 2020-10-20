#ifndef ThinWallBCProcessor_H
#define ThinWallBCProcessor_H

#include <PointerDefinitions.h>

#include "BCProcessor.h"

class ILBMKernel;

class ThinWallBCProcessor : public BCProcessor
{
public:
    ThinWallBCProcessor() = default;
    explicit ThinWallBCProcessor(SPtr<ILBMKernel> kernel);
    SPtr<BCProcessor> clone(SPtr<ILBMKernel> kernel) override;
    void applyPostCollisionBC(); // FIXME: should the base method virtual??
protected:
private:
};

#endif

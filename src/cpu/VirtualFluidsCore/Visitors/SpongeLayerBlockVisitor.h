#ifndef SpongeLayerBlockVisitor_h__
#define SpongeLayerBlockVisitor_h__

#include "Block3DVisitor.h"
#include "D3Q27System.h"

class Grid3D;
class Block3D;
class GbCuboid3D;
class LBMKernel;

//! \brief Set sponge layer for all blocks inside boundingBox
//! \details This visitor sets viscosity gradient inside bounding box.
//! \author K. Kutscher
class SpongeLayerBlockVisitor : public Block3DVisitor
{
public:
    SpongeLayerBlockVisitor(SPtr<GbCuboid3D> boundingBox, SPtr<LBMKernel> kernel, double nue, int dir);
    ~SpongeLayerBlockVisitor() override;

    void visit(SPtr<Grid3D> grid, SPtr<Block3D> block) override;

private:
    SPtr<GbCuboid3D> boundingBox;
    SPtr<LBMKernel> kernel;
    double nue;
    int dir;
};

#endif // SetSpongeLayerBlockVisitor_h__

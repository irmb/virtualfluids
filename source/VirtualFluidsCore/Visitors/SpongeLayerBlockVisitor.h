#ifndef SpongeLayerBlockVisitor_h__
#define SpongeLayerBlockVisitor_h__

#include "Block3DVisitor.h"
#include "D3Q27System.h"

class Grid3D;
class Block3D;
class GbCuboid3D;
class LBMKernel;

//! \brief Set sponge layer for all blocks inside boundingBox
//! \details This visitor is useful if you need to set or reset sponge layer in kernels (e.g. after restart because sponge layer is not serializable). 
//! \author K. Kutscher
class SpongeLayerBlockVisitor : public Block3DVisitor
{
public:
   SpongeLayerBlockVisitor(SPtr<GbCuboid3D> boundingBox, LBMReal collFactor);
   virtual ~SpongeLayerBlockVisitor();

   void visit(SPtr<Grid3D> grid, SPtr<Block3D> block) override;

   void setKernel(SPtr<LBMKernel> k);
   void setViscosity(LBMReal v);

private:
    SPtr<GbCuboid3D> boundingBox;
    LBMReal collFactor;
    SPtr<LBMKernel> kernel;
    LBMReal viscosity;
};

#endif // SetSpongeLayerBlockVisitor_h__

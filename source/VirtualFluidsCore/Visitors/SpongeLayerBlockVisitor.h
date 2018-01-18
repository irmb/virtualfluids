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
   SpongeLayerBlockVisitor();
   virtual ~SpongeLayerBlockVisitor();

   void visit(SPtr<Grid3D> grid, SPtr<Block3D> block) override;

   void setBoundingBox(SPtr<GbCuboid3D> boundingBox);
   void setKernel(SPtr<LBMKernel> k);

private:
    SPtr<GbCuboid3D> boundingBox;
    SPtr<LBMKernel> kernel;
    LBMReal viscosity;
};

#endif // SetSpongeLayerBlockVisitor_h__

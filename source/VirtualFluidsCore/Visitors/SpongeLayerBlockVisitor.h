#ifndef SpongeLayerBlockVisitor_h__
#define SpongeLayerBlockVisitor_h__

#include "Block3DVisitor.h"
#include "D3Q27System.h"

class Grid3D;
class Block3D;
class GbCuboid3D;

//! \brief Set sponge layer for all blocks inside boundingBox
//! \details This visitor is useful if you need to set or reset sponge layer in kernels (e.g. after restart because sponge layer is not serializable). 
//! \author K. Kutscher
class SpongeLayerBlockVisitor : public Block3DVisitor
{
public:
   SpongeLayerBlockVisitor(std::shared_ptr<GbCuboid3D> boundingBox, LBMReal collFactor);
   virtual ~SpongeLayerBlockVisitor();

   void visit(std::shared_ptr<Grid3D> grid, std::shared_ptr<Block3D> block) override;

private:
    std::shared_ptr<GbCuboid3D> boundingBox;
    LBMReal collFactor;
};

#endif // SetSpongeLayerBlockVisitor_h__

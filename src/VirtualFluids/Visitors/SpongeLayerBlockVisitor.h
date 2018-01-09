#ifndef SpongeLayerBlockVisitor_h__
#define SpongeLayerBlockVisitor_h__

#include "Block3DVisitor.h"

class Grid3D;
class Block3D;
class GbCuboid3D;

//! \brief Set sponge layer for all kernels of grid
//! \details This visitor is useful if you need to set or reset sponge layer in kernels (e.g. after restart because sponge layer is not serializable). 
//! \author K. Kutscher
class SpongeLayerBlockVisitor : public Block3DVisitor
{
public:
   SpongeLayerBlockVisitor(std::shared_ptr<GbCuboid3D> boundingBox);
   virtual ~SpongeLayerBlockVisitor();

   void visit(std::shared_ptr<Grid3D> grid, std::shared_ptr<Block3D> block) override;

private:
    std::shared_ptr<GbCuboid3D> boundingBox;
};

#endif // SetSpongeLayerBlockVisitor_h__

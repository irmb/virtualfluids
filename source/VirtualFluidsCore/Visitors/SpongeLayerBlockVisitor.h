#ifndef SpongeLayerBlockVisitor_h__
#define SpongeLayerBlockVisitor_h__

#include "Block3DVisitor.h"
#include "LBMKernel3D.h"
#include <numerics/geometry3d/GbCuboid3D.h>

//! \brief Set sponge layer for all kernels of grid
//! \details This visitor is useful if you need to set or reset sponge layer in kernels (e.g. after restart because sponge layer is not serializable). 
//! \author K. Kutscher

class SpongeLayerBlockVisitor : public Block3DVisitor
{
public:
   SpongeLayerBlockVisitor(GbCuboid3DPtr boundingBox);
   virtual ~SpongeLayerBlockVisitor();
   virtual void visit(Grid3DPtr grid, Block3DPtr block);
protected:
private:
   GbCuboid3DPtr boundingBox;
};

#endif // SetSpongeLayerBlockVisitor_h__

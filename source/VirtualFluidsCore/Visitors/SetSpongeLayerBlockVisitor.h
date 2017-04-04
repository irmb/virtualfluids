#ifndef SetSpongeLayerBlockVisitor_h__
#define SetSpongeLayerBlockVisitor_h__

#include "Block3DVisitor.h"
#include "LBMKernel.h"

//! \brief Set sponge layer for all kernels of grid
//! \details This visitor is useful if you need to set or reset sponge layer in kernels (e.g. after restart because sponge layer is not serializable). 
//! \author K. Kucher

class SetSpongeLayerBlockVisitor : public Block3DVisitor
{
public:
   SetSpongeLayerBlockVisitor(const mu::Parser& spongeLayer);
   virtual ~SetSpongeLayerBlockVisitor();
   virtual void visit(Grid3DPtr grid, Block3DPtr block);
protected:
private:
   mu::Parser spongeLayer;
};

#endif // SetSpongeLayerBlockVisitor_h__

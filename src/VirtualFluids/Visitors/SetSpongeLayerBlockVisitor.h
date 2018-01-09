#ifndef SetSpongeLayerBlockVisitor_h__
#define SetSpongeLayerBlockVisitor_h__

#include <memory>

#include <MuParser/include/muParser.h>

#include "Block3DVisitor.h"

class Grid3D;
class Block3D;

//! \brief Set sponge layer for all kernels of grid
//! \details This visitor is useful if you need to set or reset sponge layer in kernels (e.g. after restart because sponge layer is not serializable). 
//! \author K. Kucher
class SetSpongeLayerBlockVisitor : public Block3DVisitor
{
public:
   SetSpongeLayerBlockVisitor(const mu::Parser& spongeLayer);
   virtual ~SetSpongeLayerBlockVisitor();

   void visit(std::shared_ptr<Grid3D> grid, std::shared_ptr<Block3D> block) override;
protected:
private:
   mu::Parser spongeLayer;
};

#endif // SetSpongeLayerBlockVisitor_h__

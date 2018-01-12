#ifndef SetForcingBlockVisitor_h
#define SetForcingBlockVisitor_h

#include "Block3DVisitor.h"
#include "LBMKernel.h"

class Block3D;
class Grid3D;

//! \brief Set forcing for all kernels of grid
//! \details This visitor is useful if you need to set or reset forcing in kernels (e.g. after restart because forcing is not serializable). 
//! \author K. Kucher
class SetForcingBlockVisitor : public Block3DVisitor
{
public:
   SetForcingBlockVisitor(LBMReal forcingX1, LBMReal forcingX2, LBMReal forcingX3);
   
   SetForcingBlockVisitor(const mu::Parser& muForcingX1, const mu::Parser& muForcingX2, const mu::Parser& muForcingX3);

   SetForcingBlockVisitor(const std::string& sForcingX1, const std::string& sForcingX2, const std::string& sForcingX3);

   virtual ~SetForcingBlockVisitor() {}

   virtual void visit(SPtr<Grid3D> grid, SPtr<Block3D> block) override;

private:
   int ftype;
   LBMReal forcingX1;
   LBMReal forcingX2;
   LBMReal forcingX3;
   mu::Parser muForcingX1;
   mu::Parser muForcingX2;
   mu::Parser muForcingX3;
   std::string sForcingX1;
   std::string sForcingX2;
   std::string sForcingX3;
};

#endif

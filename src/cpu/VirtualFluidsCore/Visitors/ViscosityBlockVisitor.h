#ifndef ViscosityBlockVisitor_h
#define ViscosityBlockVisitor_h

#include <PointerDefinitions.h>

#include "Block3DVisitor.h"
#include "LBMSystem.h"

class Grid3D;
class Block3D;

class ViscosityBlockVisitor : public Block3DVisitor
{
public:
   ViscosityBlockVisitor(LBMReal nu);

   virtual ~ViscosityBlockVisitor() {}

   void visit(SPtr<Grid3D> grid, SPtr<Block3D> block) override;

private:
   LBMReal nu;
};

#endif

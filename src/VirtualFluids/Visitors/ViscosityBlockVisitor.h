#ifndef ViscosityBlockVisitor_h
#define ViscosityBlockVisitor_h

#include <memory>

#include "Block3DVisitor.h"
#include "LBMSystem.h"

class Grid3D;
class Block3D;

class ViscosityBlockVisitor : public Block3DVisitor
{
public:
   ViscosityBlockVisitor(LBMReal nu);

   virtual ~ViscosityBlockVisitor() {}

   void visit(std::shared_ptr<Grid3D> grid, std::shared_ptr<Block3D> block) override;

private:
   LBMReal nu;
};

#endif

#ifndef ViscosityBlockVisitor_h
#define ViscosityBlockVisitor_h

#include "Block3DVisitor.h"
#include "LBMKernel3D.h"


class ViscosityBlockVisitor : public Block3DVisitor
{
public:
   ViscosityBlockVisitor(LBMReal nu);

   virtual ~ViscosityBlockVisitor() {}

   virtual void visit(Grid3DPtr grid, Block3DPtr block);

private:
   LBMReal nu;
};

#endif

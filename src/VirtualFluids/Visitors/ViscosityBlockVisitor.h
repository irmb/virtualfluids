#ifndef ViscosityBlockVisitor_h
#define ViscosityBlockVisitor_h

#include "Block3DVisitor.h"
#include "LBMKernel.h"


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

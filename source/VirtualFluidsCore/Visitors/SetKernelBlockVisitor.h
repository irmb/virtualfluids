#ifndef SetKernelBlockVisitor_h
#define SetKernelBlockVisitor_h

#include "Block3DVisitor.h"
#include "LBMKernel3D.h"


class SetKernelBlockVisitor : public Block3DVisitor
{
public:
   enum Action { New, Change, ChangeBC};
public:
   //SetKernelBlockVisitor(LBMKernel3DPtr kernel, LBMReal nue);
   
   //SetKernelBlockVisitor(LBMKernel3DPtr kernel, LBMReal nue, double availMem, double needMem);

   SetKernelBlockVisitor(LBMKernel3DPtr kernel, LBMReal nue, double availMem, double needMem, SetKernelBlockVisitor::Action action = SetKernelBlockVisitor::New);

   virtual ~SetKernelBlockVisitor() {}

   virtual void visit(Grid3DPtr grid, Block3DPtr block);

   void setNoDataSetFlag(bool flag);

private:
   LBMKernel3DPtr kernel;
   LBMReal nue;
   Action action;
   bool dataSetFlag;
};

#endif

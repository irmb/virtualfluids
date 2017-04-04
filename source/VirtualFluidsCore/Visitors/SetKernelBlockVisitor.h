#ifndef SetKernelBlockVisitor_h
#define SetKernelBlockVisitor_h

#include "Block3DVisitor.h"
#include "LBMKernel.h"


class SetKernelBlockVisitor : public Block3DVisitor
{
public:
   enum Action { NewKernel, ChangeKernel, ChangeKernelWithData};
public:
   //SetKernelBlockVisitor(LBMKernel3DPtr kernel, LBMReal nue);
   
   //SetKernelBlockVisitor(LBMKernel3DPtr kernel, LBMReal nue, double availMem, double needMem);

   SetKernelBlockVisitor(LBMKernelPtr kernel, LBMReal nue, double availMem, double needMem, SetKernelBlockVisitor::Action action = SetKernelBlockVisitor::NewKernel);

   virtual ~SetKernelBlockVisitor() {}

   virtual void visit(Grid3DPtr grid, Block3DPtr block);

   void setNoDataSetFlag(bool flag);

private:
   LBMKernelPtr kernel;
   LBMReal nue;
   Action action;
   bool dataSetFlag;
};

#endif

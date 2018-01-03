#ifndef SetKernelBlockVisitor_h
#define SetKernelBlockVisitor_h

#include <memory>

#include "Block3DVisitor.h"
#include "LBMSystem.h"

class Grid3D;
class Block3D;
class LBMKernel;

class SetKernelBlockVisitor : public Block3DVisitor
{
public:
   enum Action { NewKernel, ChangeKernel, ChangeKernelWithData};

   //SetKernelBlockVisitor(LBMKernel3DPtr kernel, LBMReal nue);
   
   //SetKernelBlockVisitor(LBMKernel3DPtr kernel, LBMReal nue, double availMem, double needMem);

   SetKernelBlockVisitor(std::shared_ptr<LBMKernel> kernel, LBMReal nue, double availMem, double needMem, SetKernelBlockVisitor::Action action = SetKernelBlockVisitor::NewKernel);
   virtual ~SetKernelBlockVisitor() {}

   void visit(std::shared_ptr<Grid3D> grid, std::shared_ptr<Block3D> block) override;

   void setNoDataSetFlag(bool flag);

private:
   std::shared_ptr<LBMKernel> kernel;
   LBMReal nue;
   Action action;
   bool dataSetFlag;
};

#endif

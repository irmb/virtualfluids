#ifndef SetKernelBlockVisitor_h
#define SetKernelBlockVisitor_h

#include <PointerDefinitions.h>

#include "Block3DVisitor.h"
#include "LBMSystem.h"

class Grid3D;

class Block3D;

class LBMKernel;

class SetKernelBlockVisitor : public Block3DVisitor
{
public:
    enum Action
    {
        NewKernel, ChangeKernel, ChangeKernelWithData
    };

    //SetKernelBlockVisitor(LBMKernel3DPtr kernel, LBMReal nue);

    //SetKernelBlockVisitor(LBMKernel3DPtr kernel, LBMReal nue, double availMem, double needMem);

    SetKernelBlockVisitor(SPtr<LBMKernel> kernel, LBMReal nue, double availMem, double needMem,
                          SetKernelBlockVisitor::Action action = SetKernelBlockVisitor::NewKernel);

    SetKernelBlockVisitor(SPtr<LBMKernel> kernel, LBMReal nue, int &numberOfProcesses,
                          SetKernelBlockVisitor::Action action = SetKernelBlockVisitor::NewKernel);

    virtual ~SetKernelBlockVisitor()
    {}

    void visit(SPtr<Grid3D> grid, SPtr<Block3D> block) override;

    void setNoDataSetFlag(bool flag);

private:
    SPtr<LBMKernel> kernel;
    LBMReal nue;
    Action action;
    bool dataSetFlag;

    int numberOfProcesses{1};

    double getRequiredPhysicalMemory(const SPtr<Grid3D> &grid) const;

    void throwExceptionIfNotEnoughMemory(const SPtr<Grid3D> &grid);
};

#endif

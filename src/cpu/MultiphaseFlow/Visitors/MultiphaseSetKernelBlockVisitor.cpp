#include "MultiphaseSetKernelBlockVisitor.h"
#include "D3Q27System.h"
#include "LBMSystem.h"
#include "Block3D.h"
#include "Grid3D.h"

//SetKernelBlockVisitor::SetKernelBlockVisitor(LBMKernel3DPtr kernel, LBMReal nue) : 
//                        Block3DVisitor(0, D3Q27System::MAXLEVEL), kernel(kernel), nue(nue)
//{
//
//}
//////////////////////////////////////////////////////////////////////////
//SetKernelBlockVisitor::SetKernelBlockVisitor( LBMKernel3DPtr kernel, LBMReal nue, double availMem, double needMem ) : 
//                                              Block3DVisitor(0, D3Q27System::MAXLEVEL), kernel(kernel), nue(nue)
//{
//   if (needMem > availMem)
//   {
//      throw UbException(UB_EXARGS,"SetKernelBlockVisitor: Not enough memory!!!");
//   }
//}
//////////////////////////////////////////////////////////////////////////
MultiphaseSetKernelBlockVisitor::MultiphaseSetKernelBlockVisitor(SPtr<LBMKernel> kernel, real nuL, real nuG, real availMem, real needMem, MultiphaseSetKernelBlockVisitor::Action action /*= SetKernelBlockVisitor::New*/) :
	Block3DVisitor(0, D3Q27System::MAXLEVEL), kernel(kernel), nuL(nuL), nuG(nuG), action(action), dataSetFlag(true)
{
	if (needMem > availMem)
	{
		throw UbException(UB_EXARGS, "SetKernelBlockVisitor: Not enough memory!!!");
	}
}
//////////////////////////////////////////////////////////////////////////
void MultiphaseSetKernelBlockVisitor::visit(SPtr<Grid3D> grid, SPtr<Block3D> block)
{
	if(kernel && (block->getRank() == grid->getRank()))
	{
		real collFactorL = LBMSystem::calcCollisionFactor(nuL, block->getLevel());
		real collFactorG = LBMSystem::calcCollisionFactor(nuG, block->getLevel());
		kernel->setCollisionFactorMultiphase(collFactorL, collFactorG);

		kernel->setIndex(block->getX1(), block->getX2(), block->getX3());
		kernel->setDeltaT(LBMSystem::getDeltaT(block->getLevel()));
		kernel->setBlock(block);
        UbTupleInt3 blockNX = grid->getBlockNX();
        kernel->setNX(std::array<int, 3>{ { val<1>(blockNX), val<2>(blockNX), val<3>(blockNX) } });
        SPtr<LBMKernel> newKernel = kernel->clone();

		switch (action)
		{
		case MultiphaseSetKernelBlockVisitor::NewKernel:
			block->setKernel(newKernel);
			break;
		case MultiphaseSetKernelBlockVisitor::ChangeKernel:
			{
                SPtr<DataSet3D> dataSet = block->getKernel()->getDataSet();
				if (!dataSet)
				{
					UB_THROW(UbException(UB_EXARGS, "It is not possible to change a DataSet in kernel! Old DataSet is not exist!"));
				}

				newKernel->setDataSet(dataSet);

				SPtr<BCSet> bcProc = block->getKernel()->getBCSet();
				if (!bcProc)
				{
					UB_THROW(UbException(UB_EXARGS, "It is not possible to change a BCSet in kernel! Old BCSet is not exist!"));
				}
				newKernel->setBCSet(bcProc);
				block->setKernel(newKernel);
			}
			break;

		case MultiphaseSetKernelBlockVisitor::ChangeKernelWithData:
			{
				SPtr<BCSet> bcProc = block->getKernel()->getBCSet();
				if (!bcProc)
				{
					UB_THROW(UbException(UB_EXARGS, "It is not possible to change a BCSet in kernel! Old BCSet is not exist!"));
				}
				newKernel->setBCSet(bcProc);
				block->setKernel(newKernel);
			}
			break;
		}

	}
}

void MultiphaseSetKernelBlockVisitor::setNoDataSetFlag(bool flag)
{
	dataSetFlag = flag;
}

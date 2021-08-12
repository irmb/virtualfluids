#include "MultiphaseSetKernelBlockVisitor.h"
#include "Grid3DSystem.h"
#include "LBMSystem.h"
#include "Block3D.h"
#include "Grid3D.h"

//SetKernelBlockVisitor::SetKernelBlockVisitor(LBMKernel3DPtr kernel, LBMReal nue) : 
//                        Block3DVisitor(0, Grid3DSystem::MAXLEVEL), kernel(kernel), nue(nue)
//{
//
//}
//////////////////////////////////////////////////////////////////////////
//SetKernelBlockVisitor::SetKernelBlockVisitor( LBMKernel3DPtr kernel, LBMReal nue, double availMem, double needMem ) : 
//                                              Block3DVisitor(0, Grid3DSystem::MAXLEVEL), kernel(kernel), nue(nue)
//{
//   if (needMem > availMem)
//   {
//      throw UbException(UB_EXARGS,"SetKernelBlockVisitor: Not enough memory!!!");
//   }
//}
//////////////////////////////////////////////////////////////////////////
MultiphaseSetKernelBlockVisitor::MultiphaseSetKernelBlockVisitor(SPtr<LBMKernel> kernel, LBMReal nuL, LBMReal nuG, LBMReal densityRatio, LBMReal beta, LBMReal kappa,
	LBMReal contactAngle, double availMem, double needMem, MultiphaseSetKernelBlockVisitor::Action action /*= SetKernelBlockVisitor::New*/) :
	Block3DVisitor(0, Grid3DSystem::MAXLEVEL), kernel(kernel), nuL(nuL), nuG(nuG), densityRatio(densityRatio), beta(beta), kappa(kappa), contactAngle(contactAngle), action(action), dataSetFlag(true)
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
		LBMReal collFactorL = LBMSystem::calcCollisionFactor(nuL, block->getLevel());
		LBMReal collFactorG = LBMSystem::calcCollisionFactor(nuG, block->getLevel());

		//kernel->setCollisionFactorMultiphase(collFactorL, collFactorG);
		//kernel->setDensityRatio(densityRatio);
		//kernel->setMultiphaseModelParameters(beta, kappa);
		//kernel->setContactAngle(contactAngle);

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

				SPtr<BCProcessor> bcProc = block->getKernel()->getBCProcessor();
				if (!bcProc)
				{
					UB_THROW(UbException(UB_EXARGS, "It is not possible to change a BCProcessor in kernel! Old BCProcessor is not exist!"));
				}
				newKernel->setBCProcessor(bcProc);
				block->setKernel(newKernel);
			}
			break;

		case MultiphaseSetKernelBlockVisitor::ChangeKernelWithData:
			{
				SPtr<BCProcessor> bcProc = block->getKernel()->getBCProcessor();
				if (!bcProc)
				{
					UB_THROW(UbException(UB_EXARGS, "It is not possible to change a BCProcessor in kernel! Old BCProcessor is not exist!"));
				}
				newKernel->setBCProcessor(bcProc);
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


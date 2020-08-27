#include "ThixotropyFullDirectConnector.h"
#include "LBMKernel.h"
#include "D3Q27EsoTwist3DSplittedVector.h"
#include "DataSet3D.h"
#include "EsoTwist3D.h"


using namespace std;

ThixotropyFullDirectConnector::ThixotropyFullDirectConnector(SPtr<Block3D> from, SPtr<Block3D> to, int sendDir)
	: LocalBlock3DConnector(from, to, sendDir)

{

}
//////////////////////////////////////////////////////////////////////////
void ThixotropyFullDirectConnector::init()
{
	maxX1 = (int)this->from.lock()->getKernel()->getDataSet()->getFdistributions()->getNX1() - 1;
	maxX2 = (int)this->from.lock()->getKernel()->getDataSet()->getFdistributions()->getNX2() - 1;
	maxX3 = (int)this->from.lock()->getKernel()->getDataSet()->getFdistributions()->getNX3() - 1;

	fFrom =dynamicPointerCast<EsoTwist3D>(from.lock()->getKernel()->getDataSet()->getFdistributions());
	fTo = dynamicPointerCast<EsoTwist3D>(to.lock()->getKernel()->getDataSet()->getFdistributions());
	hFrom = dynamicPointerCast<EsoTwist3D>(from.lock()->getKernel()->getDataSet()->getHdistributions());
	hTo = dynamicPointerCast<EsoTwist3D>(to.lock()->getKernel()->getDataSet()->getHdistributions());

	//celltypematrixFrom = from.lock()->getKernel()->getBCProcessor()->getCellsTypeArray();
	//fillVolumeMatrixFrom = from.lock()->getKernel()->getBCProcessor()->getFillVolumeArray();
	//normalMatrixFrom = from.lock()->getKernel()->getBCProcessor()->getNormalArray();
	//changeInConcMatrixFrom = from.lock()->getKernel()->getBCProcessor()->getChangeInConcentrationArray();
	//deltaVolumeMatrixFrom = from.lock()->getKernel()->getBCProcessor()->getChangeInVolumeArray();

	//celltypematrixTo = to.lock()->getKernel()->getBCProcessor()->getCellsTypeArray();
	//fillVolumeMatrixTo = to.lock()->getKernel()->getBCProcessor()->getFillVolumeArray();
	//normalMatrixTo = to.lock()->getKernel()->getBCProcessor()->getNormalArray();
	//changeInConcMatrixTo = to.lock()->getKernel()->getBCProcessor()->getChangeInConcentrationArray();
	//deltaVolumeMatrixTo = to.lock()->getKernel()->getBCProcessor()->getChangeInVolumeArray();

	//bcArrayFrom = from.lock()->getKernel()->getBCProcessor()->getBCArray();
	//bcArrayTo = to.lock()->getKernel()->getBCProcessor()->getBCArray();

	//bcArrayReactiveFrom = from.lock()->getKernel()->getBCProcessor()->getBCArrayReactive();
	//bcArrayReactiveTo = to.lock()->getKernel()->getBCProcessor()->getBCArrayReactive();

	// 	maxX1 = (int)this->from.lock()->getKernel()->getBCProcessor()->getCellsTypeArray()->getNX1()/* - 1*/;
	// 	maxX2 = (int)this->from.lock()->getKernel()->getBCProcessor()->getCellsTypeArray()->getNX2()/* - 1*/;
	// 	maxX3 = (int)this->from.lock()->getKernel()->getBCProcessor()->getCellsTypeArray()->getNX3()/* - 1*/;
}
//////////////////////////////////////////////////////////////////////////
void ThixotropyFullDirectConnector::sendVectors()
{
	localDistributionsFromf = dynamicPointerCast<D3Q27EsoTwist3DSplittedVector>(this->fFrom)->getLocalDistributions();
	nonLocalDistributionsFromf = dynamicPointerCast<D3Q27EsoTwist3DSplittedVector>(this->fFrom)->getNonLocalDistributions();
	zeroDistributionsFromf = dynamicPointerCast<D3Q27EsoTwist3DSplittedVector>(this->fFrom)->getZeroDistributions();

	localDistributionsTof = dynamicPointerCast<D3Q27EsoTwist3DSplittedVector>(this->fTo)->getLocalDistributions();
	nonLocalDistributionsTof = dynamicPointerCast<D3Q27EsoTwist3DSplittedVector>(this->fTo)->getNonLocalDistributions();
	zeroDistributionsTof = dynamicPointerCast<D3Q27EsoTwist3DSplittedVector>(this->fTo)->getZeroDistributions();

	localDistributionsFromh = dynamicPointerCast<D3Q27EsoTwist3DSplittedVector>(this->hFrom)->getLocalDistributions();
	nonLocalDistributionsFromh = dynamicPointerCast<D3Q27EsoTwist3DSplittedVector>(this->hFrom)->getNonLocalDistributions();
	zeroDistributionsFromh = dynamicPointerCast<D3Q27EsoTwist3DSplittedVector>(this->hFrom)->getZeroDistributions();

	localDistributionsToh = dynamicPointerCast<D3Q27EsoTwist3DSplittedVector>(this->hTo)->getLocalDistributions();
	nonLocalDistributionsToh = dynamicPointerCast<D3Q27EsoTwist3DSplittedVector>(this->hTo)->getNonLocalDistributions();
	zeroDistributionsToh = dynamicPointerCast<D3Q27EsoTwist3DSplittedVector>(this->hTo)->getZeroDistributions();

	bool con = /*(from.lock()->getGlobalID()==11) &&*/ (to.lock()->getGlobalID() == 1);

	//EAST
	if (sendDir == D3Q27System::E)
	{
		for (int x3 = 1; x3 < maxX3; x3++)
		{
			for (int x2 = 1; x2 < maxX2; x2++)
			{
				exchangeData(maxX1 - 1, x2, x3, 0, x2, x3);
				//exchangeDataforCells(maxX1 - 1, x2, x3, 0, x2, x3);
			}
		}
	}
	//WEST
	else if (sendDir == D3Q27System::W)
	{
		for (int x3 = 1; x3 < maxX3; x3++)
		{
			for (int x2 = 1; x2 < maxX2; x2++)
			{
				exchangeData(1, x2, x3, maxX1, x2, x3);
				//exchangeDataforCells(1, x2, x3, maxX1, x2, x3);
			}
		}
	}
	//NORTH
	else if (sendDir == D3Q27System::N)
	{
		for (int x3 = 1; x3 < maxX3; x3++)
		{
			for (int x1 = 1; x1 < maxX1; x1++)
			{
				exchangeData(x1, maxX2 - 1, x3, x1, 0, x3);
				//exchangeDataforCells(x1, maxX2 - 1, x3, x1, 0, x3);
			}
		}
	}
	//SOUTH
	else if (sendDir == D3Q27System::S)
	{
		for (int x3 = 1; x3 < maxX3; x3++)
		{
			for (int x1 = 1; x1 < maxX1; x1++)
			{
				exchangeData(x1, 1, x3, x1, maxX2, x3);
				//exchangeDataforCells(x1, 1, x3, x1, maxX2, x3);
			}
		}
	}

	//TOP
	else if (sendDir == D3Q27System::T)
	{
		for (int x2 = 1; x2 < maxX2; x2++)
		{
			for (int x1 = 1; x1 < maxX1; x1++)
			{
				exchangeData(x1, x2, maxX3 - 1, x1, x2, 0);
				//exchangeDataforCells(x1, x2, maxX3 - 1, x1, x2, 0);
			}
		}
	}
	//BOTTOM
	else if (sendDir == D3Q27System::B)
	{
		for (int x2 = 1; x2 < maxX2; x2++)
		{
			for (int x1 = 1; x1 < maxX1; x1++)
			{
				exchangeData(x1, x2, 1, x1, x2, maxX3);
				//exchangeDataforCells(x1, x2, 1, x1, x2, maxX3);
			}
		}
	}
	//NORTHEAST
	else if (sendDir == D3Q27System::NE)
	{
		for (int x3 = 1; x3 < maxX3; x3++)
		{
			exchangeData(maxX1 - 1, maxX2 - 1, x3, 0, 0, x3);
			//exchangeDataforCells(maxX1 - 1, maxX2 - 1, x3, 0, 0, x3);

		}
	}
	//NORTHWEST
	else if (sendDir == D3Q27System::NW)
	{
		for (int x3 = 1; x3 < maxX3; x3++)
		{
			exchangeData(1, maxX2 - 1, x3, maxX1, 0, x3);
			//exchangeDataforCells(1, maxX2 - 1, x3, maxX1, 0, x3);
		}
	}
	//SOUTHWEST
	else if (sendDir == D3Q27System::SW)
	{
		for (int x3 = 1; x3 < maxX3; x3++)
		{
			exchangeData(1, 1, x3, maxX1, maxX2, x3);
			//exchangeDataforCells(1, 1, x3, maxX1, maxX2, x3);

		}
	}
	//SOUTHEAST
	else if (sendDir == D3Q27System::SE)
	{
		for (int x3 = 1; x3 < maxX3; x3++)
		{
			exchangeData(maxX1 - 1, 1, x3, 0, maxX2, x3);
			//exchangeDataforCells(maxX1 - 1, 1, x3, 0, maxX2, x3);

		}
	}
	else if (sendDir == D3Q27System::TE)
		for (int x2 = 1; x2 < maxX2; x2++)
		{
			exchangeData(maxX1 - 1, x2, maxX3 - 1, 0, x2, 0);
			//exchangeDataforCells(maxX1 - 1, x2, maxX3 - 1, 0, x2, 0);

		}
	else if (sendDir == D3Q27System::BW)
		for (int x2 = 1; x2 < maxX2; x2++)
		{
			exchangeData(1, x2, 1, maxX1, x2, maxX3);
			//exchangeDataforCells(1, x2, 1, maxX1, x2, maxX3);

		}
	else if (sendDir == D3Q27System::BE)
		for (int x2 = 1; x2 < maxX2; x2++)
		{
			exchangeData(maxX1 - 1, x2, 1, 0, x2, maxX3);
			//exchangeDataforCells(maxX1 - 1, x2, 1, 0, x2, maxX3);

		}
	else if (sendDir == D3Q27System::TW)
		for (int x2 = 1; x2 < maxX2; x2++)
		{
			exchangeData(1, x2, maxX3 - 1, maxX1, x2, 0);
			//exchangeDataforCells(1, x2, maxX3 - 1, maxX1, x2, 0);

		}
	else if (sendDir == D3Q27System::TN)
		for (int x1 = 1; x1 < maxX1; x1++)
		{
			exchangeData(x1, maxX2 - 1, maxX3 - 1, x1, 0, 0);
			//exchangeDataforCells(x1, maxX2 - 1, maxX3 - 1, x1, 0, 0);

		}
	else if (sendDir == D3Q27System::BS)
		for (int x1 = 1; x1 < maxX1; x1++)
		{
			exchangeData(x1, 1, 1, x1, maxX2, maxX3);
			//exchangeDataforCells(x1, 1, 1, x1, maxX2, maxX3);

		}
	else if (sendDir == D3Q27System::BN)
		for (int x1 = 1; x1 < maxX1; x1++)
		{
			exchangeData(x1, maxX2 - 1, 1, x1, 0, maxX3);
			//exchangeDataforCells(x1, maxX2 - 1, 1, x1, 0, maxX3);

		}

	else if (sendDir == D3Q27System::TS)
		for (int x1 = 1; x1 < maxX1; x1++)
		{
			exchangeData(x1, 1, maxX3 - 1, x1, maxX2, 0);
			//exchangeDataforCells(x1, 1, maxX3 - 1, x1, maxX2, 0);

		}

	else if (sendDir == D3Q27System::TSW)
	{
		exchangeData(1, 1, maxX3 - 1, maxX1, maxX2, 0);
		//exchangeDataforCells(1, 1, maxX3 - 1, maxX1, maxX2, 0);

	}
	else if (sendDir == D3Q27System::TSE)
	{
		exchangeData(maxX1 - 1, 1, maxX3 - 1, 0, maxX2, 0);
		//exchangeDataforCells(maxX1 - 1, 1, maxX3 - 1, 0, maxX2, 0);

	}
	else if (sendDir == D3Q27System::TNW)
	{
		exchangeData(1, maxX2 - 1, maxX3 - 1, maxX1, 0, 0);
		//exchangeDataforCells(1, maxX2 - 1, maxX3 - 1, maxX1, 0, 0);

	}
	else if (sendDir == D3Q27System::TNE)
	{
		exchangeData(maxX1 - 1, maxX2 - 1, maxX3 - 1, 0, 0, 0);
		//exchangeDataforCells(maxX1 - 1, maxX2 - 1, maxX3 - 1, 0, 0, 0);

	}
	else if (sendDir == D3Q27System::BSW)
	{
		exchangeData(1, 1, 1, maxX1, maxX2, maxX3);
		//exchangeDataforCells(1, 1, 1, maxX1, maxX2, maxX3);

	}
	else if (sendDir == D3Q27System::BSE)
	{
		exchangeData(maxX1 - 1, 1, 1, 0, maxX2, maxX3);
		//exchangeDataforCells(maxX1 - 1, 1, 1, 0, maxX2, maxX3);

	}
	else if (sendDir == D3Q27System::BNW)
	{
		exchangeData(1, maxX2 - 1, 1, maxX1, 0, maxX3);
		//exchangeDataforCells(1, maxX2 - 1, 1, maxX1, 0, maxX3);

	}
	else if (sendDir == D3Q27System::BNE)
	{
		exchangeData(maxX1 - 1, maxX2 - 1, 1, 0, 0, maxX3);
		//exchangeDataforCells(maxX1 - 1, maxX2 - 1, 1, 0, 0, maxX3);

	}
	else UB_THROW(UbException(UB_EXARGS, "unknown dir"));
}


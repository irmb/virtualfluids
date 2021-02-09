/**
* @file ThixotropyFullDirectConnector.h
* @brief Connector send and receive full distribution in shared memory
*
* @author Hussein Alihussein
* @date 28.11.2018
*/
#ifndef ThixotropyFullDirectConnector_H
#define ThixotropyFullDirectConnector_H

#include "LocalBlock3DConnector.h"
#include "Block3D.h"
#include "D3Q27System.h"
#include "basics/container/CbArray3D.h"
#include "basics/container/CbArray4D.h"
#include "BCArray3D.h"

class EsoTwist3D;

//! \brief   Exchange data between blocks. 
//! \details Connector send and receive full distributions between two blocks in shared memory.
//! \author  Hussein Alihussein

class ThixotropyFullDirectConnector : public LocalBlock3DConnector
{
public:
	ThixotropyFullDirectConnector(SPtr<Block3D> from, SPtr<Block3D> to, int sendDir);
	void init();
	void sendVectors();
	//void sendCellVectors();

protected:
	inline void exchangeData(int x1From, int x2From, int x3From, int x1To, int x2To, int x3To);
	//inline void exchangeDataforCells(int x1From, int x2From, int x3From, int x1To, int x2To, int x3To);
	//void ExchangeDatabcArray(SPtr<BCArray3D> bcArrayFrom, SPtr<BCArray3D> bcArrayTo, int x1From, int x2From, int x3From, int x1To, int x2To, int x3To);
	//void ExchangeDatabcArrayReactive(BCArray3DReactivePtr bcArrayFrom, BCArray3DReactivePtr bcArrayTo, int x1From, int x2From, int x3From, int x1To, int x2To, int x3To);
private:
	int maxX1;
	int maxX2;
	int maxX3;

	//CbArray3D<int>::CbArray3DPtr celltypematrixFrom;
	//CbArray3D<double>::CbArray3DPtr fillVolumeMatrixFrom;
	//CbArray4D<LBMReal, IndexerX4X3X2X1>::CbArray4DPtr normalMatrixFrom;
	//CbArray3D<double>::CbArray3DPtr changeInConcMatrixFrom;
	//CbArray3D<double>::CbArray3DPtr deltaVolumeMatrixFrom;

	//CbArray3D<int>::CbArray3DPtr celltypematrixTo;
	//CbArray3D<double>::CbArray3DPtr fillVolumeMatrixTo;
	//CbArray4D<LBMReal, IndexerX4X3X2X1>::CbArray4DPtr normalMatrixTo;
	//CbArray3D<double>::CbArray3DPtr changeInConcMatrixTo;
	//CbArray3D<double>::CbArray3DPtr deltaVolumeMatrixTo;

	CbArray4D <LBMReal, IndexerX4X3X2X1>::CbArray4DPtr localDistributionsFromf;
	CbArray4D <LBMReal, IndexerX4X3X2X1>::CbArray4DPtr nonLocalDistributionsFromf;
	CbArray3D <LBMReal, IndexerX3X2X1>::CbArray3DPtr   zeroDistributionsFromf;

	CbArray4D <LBMReal, IndexerX4X3X2X1>::CbArray4DPtr localDistributionsTof;
	CbArray4D <LBMReal, IndexerX4X3X2X1>::CbArray4DPtr nonLocalDistributionsTof;
	CbArray3D <LBMReal, IndexerX3X2X1>::CbArray3DPtr   zeroDistributionsTof;

	CbArray4D <LBMReal, IndexerX4X3X2X1>::CbArray4DPtr localDistributionsFromh;
	CbArray4D <LBMReal, IndexerX4X3X2X1>::CbArray4DPtr nonLocalDistributionsFromh;
	CbArray3D <LBMReal, IndexerX3X2X1>::CbArray3DPtr   zeroDistributionsFromh;

	CbArray4D <LBMReal, IndexerX4X3X2X1>::CbArray4DPtr localDistributionsToh;
	CbArray4D <LBMReal, IndexerX4X3X2X1>::CbArray4DPtr nonLocalDistributionsToh;
	CbArray3D <LBMReal, IndexerX3X2X1>::CbArray3DPtr   zeroDistributionsToh;

	SPtr<EsoTwist3D>  fFrom, hFrom;
	SPtr<EsoTwist3D>  fTo, hTo;

	//BCArray3DPtr bcArrayFrom;
	//BCArray3DPtr bcArrayTo;

	//BCArray3DReactivePtr bcArrayReactiveFrom;
	//BCArray3DReactivePtr bcArrayReactiveTo;
};


//////////////////////////////////////////////////////////////////////////
inline void ThixotropyFullDirectConnector::exchangeData(int x1From, int x2From, int x3From, int x1To, int x2To, int x3To)
{


	(*this->localDistributionsTof)(D3Q27System::ET_E, x1To, x2To, x3To) = (*this->localDistributionsFromf)(D3Q27System::ET_E, x1From, x2From, x3From);
	(*this->localDistributionsTof)(D3Q27System::ET_N, x1To, x2To, x3To) = (*this->localDistributionsFromf)(D3Q27System::ET_N, x1From, x2From, x3From);
	(*this->localDistributionsTof)(D3Q27System::ET_T, x1To, x2To, x3To) = (*this->localDistributionsFromf)(D3Q27System::ET_T, x1From, x2From, x3From);
	(*this->localDistributionsTof)(D3Q27System::ET_NE, x1To, x2To, x3To) = (*this->localDistributionsFromf)(D3Q27System::ET_NE, x1From, x2From, x3From);
	(*this->localDistributionsTof)(D3Q27System::ET_NW, x1To + 1, x2To, x3To) = (*this->localDistributionsFromf)(D3Q27System::ET_NW, x1From + 1, x2From, x3From);
	(*this->localDistributionsTof)(D3Q27System::ET_TE, x1To, x2To, x3To) = (*this->localDistributionsFromf)(D3Q27System::ET_TE, x1From, x2From, x3From);
	(*this->localDistributionsTof)(D3Q27System::ET_TW, x1To + 1, x2To, x3To) = (*this->localDistributionsFromf)(D3Q27System::ET_TW, x1From + 1, x2From, x3From);
	(*this->localDistributionsTof)(D3Q27System::ET_TN, x1To, x2To, x3To) = (*this->localDistributionsFromf)(D3Q27System::ET_TN, x1From, x2From, x3From);
	(*this->localDistributionsTof)(D3Q27System::ET_TS, x1To, x2To + 1, x3To) = (*this->localDistributionsFromf)(D3Q27System::ET_TS, x1From, x2From + 1, x3From);
	(*this->localDistributionsTof)(D3Q27System::ET_TNE, x1To, x2To, x3To) = (*this->localDistributionsFromf)(D3Q27System::ET_TNE, x1From, x2From, x3From);
	(*this->localDistributionsTof)(D3Q27System::ET_TNW, x1To + 1, x2To, x3To) = (*this->localDistributionsFromf)(D3Q27System::ET_TNW, x1From + 1, x2From, x3From);
	(*this->localDistributionsTof)(D3Q27System::ET_TSE, x1To, x2To + 1, x3To) = (*this->localDistributionsFromf)(D3Q27System::ET_TSE, x1From, x2From + 1, x3From);
	(*this->localDistributionsTof)(D3Q27System::ET_TSW, x1To + 1, x2To + 1, x3To) = (*this->localDistributionsFromf)(D3Q27System::ET_TSW, x1From + 1, x2From + 1, x3From);

	(*this->nonLocalDistributionsTof)(D3Q27System::ET_W, x1To + 1, x2To, x3To) = (*this->nonLocalDistributionsFromf)(D3Q27System::ET_W, x1From + 1, x2From, x3From);
	(*this->nonLocalDistributionsTof)(D3Q27System::ET_S, x1To, x2To + 1, x3To) = (*this->nonLocalDistributionsFromf)(D3Q27System::ET_S, x1From, x2From + 1, x3From);
	(*this->nonLocalDistributionsTof)(D3Q27System::ET_B, x1To, x2To, x3To + 1) = (*this->nonLocalDistributionsFromf)(D3Q27System::ET_B, x1From, x2From, x3From + 1);
	(*this->nonLocalDistributionsTof)(D3Q27System::ET_SW, x1To + 1, x2To + 1, x3To) = (*this->nonLocalDistributionsFromf)(D3Q27System::ET_SW, x1From + 1, x2From + 1, x3From);
	(*this->nonLocalDistributionsTof)(D3Q27System::ET_SE, x1To, x2To + 1, x3To) = (*this->nonLocalDistributionsFromf)(D3Q27System::ET_SE, x1From, x2From + 1, x3From);
	(*this->nonLocalDistributionsTof)(D3Q27System::ET_BW, x1To + 1, x2To, x3To + 1) = (*this->nonLocalDistributionsFromf)(D3Q27System::ET_BW, x1From + 1, x2From, x3From + 1);
	(*this->nonLocalDistributionsTof)(D3Q27System::ET_BE, x1To, x2To, x3To + 1) = (*this->nonLocalDistributionsFromf)(D3Q27System::ET_BE, x1From, x2From, x3From + 1);
	(*this->nonLocalDistributionsTof)(D3Q27System::ET_BS, x1To, x2To + 1, x3To + 1) = (*this->nonLocalDistributionsFromf)(D3Q27System::ET_BS, x1From, x2From + 1, x3From + 1);
	(*this->nonLocalDistributionsTof)(D3Q27System::ET_BN, x1To, x2To, x3To + 1) = (*this->nonLocalDistributionsFromf)(D3Q27System::ET_BN, x1From, x2From, x3From + 1);
	(*this->nonLocalDistributionsTof)(D3Q27System::ET_BSW, x1To + 1, x2To + 1, x3To + 1) = (*this->nonLocalDistributionsFromf)(D3Q27System::ET_BSW, x1From + 1, x2From + 1, x3From + 1);
	(*this->nonLocalDistributionsTof)(D3Q27System::ET_BSE, x1To, x2To + 1, x3To + 1) = (*this->nonLocalDistributionsFromf)(D3Q27System::ET_BSE, x1From, x2From + 1, x3From + 1);
	(*this->nonLocalDistributionsTof)(D3Q27System::ET_BNW, x1To + 1, x2To, x3To + 1) = (*this->nonLocalDistributionsFromf)(D3Q27System::ET_BNW, x1From + 1, x2From, x3From + 1);
	(*this->nonLocalDistributionsTof)(D3Q27System::ET_BNE, x1To, x2To, x3To + 1) = (*this->nonLocalDistributionsFromf)(D3Q27System::ET_BNE, x1From, x2From, x3From + 1);

	(*this->zeroDistributionsTof)(x1To, x2To, x3To) = (*this->zeroDistributionsFromf)(x1From, x2From, x3From);


	(*this->localDistributionsToh)(D3Q27System::ET_E, x1To, x2To, x3To) = (*this->localDistributionsFromh)(D3Q27System::ET_E, x1From, x2From, x3From);
	(*this->localDistributionsToh)(D3Q27System::ET_N, x1To, x2To, x3To) = (*this->localDistributionsFromh)(D3Q27System::ET_N, x1From, x2From, x3From);
	(*this->localDistributionsToh)(D3Q27System::ET_T, x1To, x2To, x3To) = (*this->localDistributionsFromh)(D3Q27System::ET_T, x1From, x2From, x3From);
	(*this->localDistributionsToh)(D3Q27System::ET_NE, x1To, x2To, x3To) = (*this->localDistributionsFromh)(D3Q27System::ET_NE, x1From, x2From, x3From);
	(*this->localDistributionsToh)(D3Q27System::ET_NW, x1To + 1, x2To, x3To) = (*this->localDistributionsFromh)(D3Q27System::ET_NW, x1From + 1, x2From, x3From);
	(*this->localDistributionsToh)(D3Q27System::ET_TE, x1To, x2To, x3To) = (*this->localDistributionsFromh)(D3Q27System::ET_TE, x1From, x2From, x3From);
	(*this->localDistributionsToh)(D3Q27System::ET_TW, x1To + 1, x2To, x3To) = (*this->localDistributionsFromh)(D3Q27System::ET_TW, x1From + 1, x2From, x3From);
	(*this->localDistributionsToh)(D3Q27System::ET_TN, x1To, x2To, x3To) = (*this->localDistributionsFromh)(D3Q27System::ET_TN, x1From, x2From, x3From);
	(*this->localDistributionsToh)(D3Q27System::ET_TS, x1To, x2To + 1, x3To) = (*this->localDistributionsFromh)(D3Q27System::ET_TS, x1From, x2From + 1, x3From);
	(*this->localDistributionsToh)(D3Q27System::ET_TNE, x1To, x2To, x3To) = (*this->localDistributionsFromh)(D3Q27System::ET_TNE, x1From, x2From, x3From);
	(*this->localDistributionsToh)(D3Q27System::ET_TNW, x1To + 1, x2To, x3To) = (*this->localDistributionsFromh)(D3Q27System::ET_TNW, x1From + 1, x2From, x3From);
	(*this->localDistributionsToh)(D3Q27System::ET_TSE, x1To, x2To + 1, x3To) = (*this->localDistributionsFromh)(D3Q27System::ET_TSE, x1From, x2From + 1, x3From);
	(*this->localDistributionsToh)(D3Q27System::ET_TSW, x1To + 1, x2To + 1, x3To) = (*this->localDistributionsFromh)(D3Q27System::ET_TSW, x1From + 1, x2From + 1, x3From);

	(*this->nonLocalDistributionsToh)(D3Q27System::ET_W, x1To + 1, x2To, x3To) = (*this->nonLocalDistributionsFromh)(D3Q27System::ET_W, x1From + 1, x2From, x3From);
	(*this->nonLocalDistributionsToh)(D3Q27System::ET_S, x1To, x2To + 1, x3To) = (*this->nonLocalDistributionsFromh)(D3Q27System::ET_S, x1From, x2From + 1, x3From);
	(*this->nonLocalDistributionsToh)(D3Q27System::ET_B, x1To, x2To, x3To + 1) = (*this->nonLocalDistributionsFromh)(D3Q27System::ET_B, x1From, x2From, x3From + 1);
	(*this->nonLocalDistributionsToh)(D3Q27System::ET_SW, x1To + 1, x2To + 1, x3To) = (*this->nonLocalDistributionsFromh)(D3Q27System::ET_SW, x1From + 1, x2From + 1, x3From);
	(*this->nonLocalDistributionsToh)(D3Q27System::ET_SE, x1To, x2To + 1, x3To) = (*this->nonLocalDistributionsFromh)(D3Q27System::ET_SE, x1From, x2From + 1, x3From);
	(*this->nonLocalDistributionsToh)(D3Q27System::ET_BW, x1To + 1, x2To, x3To + 1) = (*this->nonLocalDistributionsFromh)(D3Q27System::ET_BW, x1From + 1, x2From, x3From + 1);
	(*this->nonLocalDistributionsToh)(D3Q27System::ET_BE, x1To, x2To, x3To + 1) = (*this->nonLocalDistributionsFromh)(D3Q27System::ET_BE, x1From, x2From, x3From + 1);
	(*this->nonLocalDistributionsToh)(D3Q27System::ET_BS, x1To, x2To + 1, x3To + 1) = (*this->nonLocalDistributionsFromh)(D3Q27System::ET_BS, x1From, x2From + 1, x3From + 1);
	(*this->nonLocalDistributionsToh)(D3Q27System::ET_BN, x1To, x2To, x3To + 1) = (*this->nonLocalDistributionsFromh)(D3Q27System::ET_BN, x1From, x2From, x3From + 1);
	(*this->nonLocalDistributionsToh)(D3Q27System::ET_BSW, x1To + 1, x2To + 1, x3To + 1) = (*this->nonLocalDistributionsFromh)(D3Q27System::ET_BSW, x1From + 1, x2From + 1, x3From + 1);
	(*this->nonLocalDistributionsToh)(D3Q27System::ET_BSE, x1To, x2To + 1, x3To + 1) = (*this->nonLocalDistributionsFromh)(D3Q27System::ET_BSE, x1From, x2From + 1, x3From + 1);
	(*this->nonLocalDistributionsToh)(D3Q27System::ET_BNW, x1To + 1, x2To, x3To + 1) = (*this->nonLocalDistributionsFromh)(D3Q27System::ET_BNW, x1From + 1, x2From, x3From + 1);
	(*this->nonLocalDistributionsToh)(D3Q27System::ET_BNE, x1To, x2To, x3To + 1) = (*this->nonLocalDistributionsFromh)(D3Q27System::ET_BNE, x1From, x2From, x3From + 1);

	(*this->zeroDistributionsToh)(x1To, x2To, x3To) = (*this->zeroDistributionsFromh)(x1From, x2From, x3From);
}
//////////////////////////////////////////////////////////////////////////
//inline void ThixotropyFullDirectConnector::exchangeDataforCells(int x1From, int x2From, int x3From, int x1To, int x2To, int x3To)
//{
//	//(*this->celltypematrixTo)( x1To, x2To, x3To) = (*this->celltypematrixFrom)( x1From, x2From, x3From);
//
//	(*this->celltypematrixTo)(x1To, x2To, x3To) = (*this->celltypematrixFrom)(x1From, x2From, x3From);
//	(*this->fillVolumeMatrixTo)(x1To, x2To, x3To) = (*this->fillVolumeMatrixFrom)(x1From, x2From, x3From);
//	(*this->changeInConcMatrixTo)(x1To, x2To, x3To) = (*this->changeInConcMatrixFrom)(x1From, x2From, x3From);
//	(*this->deltaVolumeMatrixTo)(x1To, x2To, x3To) = (*this->deltaVolumeMatrixFrom)(x1From, x2From, x3From);
//
//	(*this->normalMatrixTo)(0, x1To, x2To, x3To) = (*this->normalMatrixFrom)(0, x1From, x2From, x3From);
//	(*this->normalMatrixTo)(1, x1To, x2To, x3To) = (*this->normalMatrixFrom)(1, x1From, x2From, x3From);
//	(*this->normalMatrixTo)(2, x1To, x2To, x3To) = (*this->normalMatrixFrom)(2, x1From, x2From, x3From);
//	//bcArrayTo->bcindexmatrix(x1To, x2To, x3To) = bcArrayFrom->bcindexmatrix(x1From, x2From, x3From);
//	ExchangeDatabcArray(bcArrayFrom, bcArrayTo, x1From, x2From, x3From, x1To, x2To, x3To);
//	ExchangeDatabcArrayReactive(bcArrayReactiveFrom, bcArrayReactiveTo, x1From, x2From, x3From, x1To, x2To, x3To);
//}
////////////////////////////////////////////////////////////////////////////
//inline void ThixotropyFullDirectConnector::ExchangeDatabcArray(BCArray3DPtr bcArrayFrom, BCArray3DPtr bcArrayTo, int x1From, int x2From, int x3From, int x1To, int x2To, int x3To)
//{
//	if (bcArrayFrom->isFluid(x1From, x2From, x3From))
//	{
//		if (bcArrayFrom->hasBC(x1From, x2From, x3From))
//		{
//			BoundaryConditionsPtr bc = bcArrayFrom->getBC(x1From, x2From, x3From);
//
//			bcArrayTo->setBC(x1To, x2To, x3To, bc);
//
//		}
//		else if (bcArrayFrom->isFluidWithoutBC(x1From, x2From, x3From))
//		{
//			bcArrayTo->setFluid(x1To, x2To, x3To);
//		}
//	}
//	else if (bcArrayFrom->isSolid(x1From, x2From, x3From))
//	{
//		bcArrayTo->setSolid(x1To, x2To, x3To);
//	}
//	else
//	{
//		bcArrayTo->setUndefined(x1To, x2To, x3To);
//	}
//
//}
////////////////////////////////////////////////////////////////////////////
//inline void ThixotropyFullDirectConnector::ExchangeDatabcArrayReactive(BCArray3DReactivePtr bcArrayFrom, BCArray3DReactivePtr bcArrayTo, int x1From, int x2From, int x3From, int x1To, int x2To, int x3To)
//{
//	if (bcArrayFrom->isFluid(x1From, x2From, x3From))
//	{
//		if (bcArrayFrom->hasBC(x1From, x2From, x3From))
//		{
//			BoundaryConditionsPtr bc = bcArrayFrom->getBC(x1From, x2From, x3From);
//
//			bcArrayTo->setBC(x1To, x2To, x3To, bc);
//
//		}
//		else if (bcArrayFrom->isFluidWithoutBC(x1From, x2From, x3From))
//		{
//			bcArrayTo->setFluid(x1To, x2To, x3To);
//		}
//	}
//	else if (bcArrayFrom->isSolid(x1From, x2From, x3From))
//	{
//		bcArrayTo->setSolid(x1To, x2To, x3To);
//	}
//	else if (bcArrayFrom->isReactiveSolid(x1From, x2From, x3From))
//	{
//		bcArrayTo->setReactiveSolidType1(x1To, x2To, x3To);
//	}
//	else
//	{
//		bcArrayTo->setUndefined(x1To, x2To, x3To);
//	}
//
//}
#endif
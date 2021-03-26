#ifndef OneDistributionFullVectorConnector_H
#define OneDistributionFullVectorConnector_H

#include <vector>

#include "Block3D.h"
#include "D3Q27System.h"
#include "EsoTwist3D.h"
//#include "EsoTwistD3Q27System.h"
#include "LBMKernel.h"
#include "FullVectorConnector.h"
#include "D3Q27EsoTwist3DSplittedVector.h"
#include "basics/container/CbArray3D.h"
#include "basics/container/CbArray4D.h"

// daten werden in einen vector (dieser befindet sich im transmitter) kopiert
// der vector wird via transmitter uebertragen
// transmitter kann ein lokal, MPI, RCG, CTL oder was auch immer fuer ein
// transmitter sein, der von Transmitter abgeleitet ist ;-)
class OneDistributionFullVectorConnector : public FullVectorConnector
{
public:
    OneDistributionFullVectorConnector(SPtr<Block3D> block, VectorTransmitterPtr sender, VectorTransmitterPtr receiver,
                               int sendDir);

    void init() override;

protected:
    inline void updatePointers() override;
    inline void fillData(vector_type &sdata, int &index, int x1, int x2, int x3) override;
    inline void distributeData(vector_type &rdata, int &index, int x1, int x2, int x3) override;

private:
    CbArray4D<LBMReal, IndexerX4X3X2X1>::CbArray4DPtr localDistributions;
    CbArray4D<LBMReal, IndexerX4X3X2X1>::CbArray4DPtr nonLocalDistributions;
    CbArray3D<LBMReal, IndexerX3X2X1>::CbArray3DPtr zeroDistributions;

    SPtr<EsoTwist3D> fDis;
};
//////////////////////////////////////////////////////////////////////////
inline void OneDistributionFullVectorConnector::updatePointers()
{
    localDistributions    = dynamicPointerCast<D3Q27EsoTwist3DSplittedVector>(this->fDis)->getLocalDistributions();
    nonLocalDistributions = dynamicPointerCast<D3Q27EsoTwist3DSplittedVector>(this->fDis)->getNonLocalDistributions();
    zeroDistributions     = dynamicPointerCast<D3Q27EsoTwist3DSplittedVector>(this->fDis)->getZeroDistributions();
}
//////////////////////////////////////////////////////////////////////////
inline void OneDistributionFullVectorConnector::fillData(vector_type &sdata, int &index, int x1, int x2, int x3)
{
    sdata[index++] = (*this->localDistributions)(D3Q27System::ET_E, x1, x2, x3);
    sdata[index++] = (*this->localDistributions)(D3Q27System::ET_N, x1, x2, x3);
    sdata[index++] = (*this->localDistributions)(D3Q27System::ET_T, x1, x2, x3);
    sdata[index++] = (*this->localDistributions)(D3Q27System::ET_NE, x1, x2, x3);
    sdata[index++] = (*this->localDistributions)(D3Q27System::ET_NW, x1 + 1, x2, x3);
    sdata[index++] = (*this->localDistributions)(D3Q27System::ET_TE, x1, x2, x3);
    sdata[index++] = (*this->localDistributions)(D3Q27System::ET_TW, x1 + 1, x2, x3);
    sdata[index++] = (*this->localDistributions)(D3Q27System::ET_TN, x1, x2, x3);
    sdata[index++] = (*this->localDistributions)(D3Q27System::ET_TS, x1, x2 + 1, x3);
    sdata[index++] = (*this->localDistributions)(D3Q27System::ET_TNE, x1, x2, x3);
    sdata[index++] = (*this->localDistributions)(D3Q27System::ET_TNW, x1 + 1, x2, x3);
    sdata[index++] = (*this->localDistributions)(D3Q27System::ET_TSE, x1, x2 + 1, x3);
    sdata[index++] = (*this->localDistributions)(D3Q27System::ET_TSW, x1 + 1, x2 + 1, x3);

    sdata[index++] = (*this->nonLocalDistributions)(D3Q27System::ET_W, x1 + 1, x2, x3);
    sdata[index++] = (*this->nonLocalDistributions)(D3Q27System::ET_S, x1, x2 + 1, x3);
    sdata[index++] = (*this->nonLocalDistributions)(D3Q27System::ET_B, x1, x2, x3 + 1);
    sdata[index++] = (*this->nonLocalDistributions)(D3Q27System::ET_SW, x1 + 1, x2 + 1, x3);
    sdata[index++] = (*this->nonLocalDistributions)(D3Q27System::ET_SE, x1, x2 + 1, x3);
    sdata[index++] = (*this->nonLocalDistributions)(D3Q27System::ET_BW, x1 + 1, x2, x3 + 1);
    sdata[index++] = (*this->nonLocalDistributions)(D3Q27System::ET_BE, x1, x2, x3 + 1);
    sdata[index++] = (*this->nonLocalDistributions)(D3Q27System::ET_BS, x1, x2 + 1, x3 + 1);
    sdata[index++] = (*this->nonLocalDistributions)(D3Q27System::ET_BN, x1, x2, x3 + 1);
    sdata[index++] = (*this->nonLocalDistributions)(D3Q27System::ET_BSW, x1 + 1, x2 + 1, x3 + 1);
    sdata[index++] = (*this->nonLocalDistributions)(D3Q27System::ET_BSE, x1, x2 + 1, x3 + 1);
    sdata[index++] = (*this->nonLocalDistributions)(D3Q27System::ET_BNW, x1 + 1, x2, x3 + 1);
    sdata[index++] = (*this->nonLocalDistributions)(D3Q27System::ET_BNE, x1, x2, x3 + 1);

    sdata[index++] = (*this->zeroDistributions)(x1, x2, x3);
}
//////////////////////////////////////////////////////////////////////////
inline void OneDistributionFullVectorConnector::distributeData(vector_type &rdata, int &index, int x1, int x2, int x3)
{
    (*this->localDistributions)(D3Q27System::ET_E, x1, x2, x3)           = rdata[index++];
    (*this->localDistributions)(D3Q27System::ET_N, x1, x2, x3)           = rdata[index++];
    (*this->localDistributions)(D3Q27System::ET_T, x1, x2, x3)           = rdata[index++];
    (*this->localDistributions)(D3Q27System::ET_NE, x1, x2, x3)          = rdata[index++];
    (*this->localDistributions)(D3Q27System::ET_NW, x1 + 1, x2, x3)      = rdata[index++];
    (*this->localDistributions)(D3Q27System::ET_TE, x1, x2, x3)          = rdata[index++];
    (*this->localDistributions)(D3Q27System::ET_TW, x1 + 1, x2, x3)      = rdata[index++];
    (*this->localDistributions)(D3Q27System::ET_TN, x1, x2, x3)          = rdata[index++];
    (*this->localDistributions)(D3Q27System::ET_TS, x1, x2 + 1, x3)      = rdata[index++];
    (*this->localDistributions)(D3Q27System::ET_TNE, x1, x2, x3)         = rdata[index++];
    (*this->localDistributions)(D3Q27System::ET_TNW, x1 + 1, x2, x3)     = rdata[index++];
    (*this->localDistributions)(D3Q27System::ET_TSE, x1, x2 + 1, x3)     = rdata[index++];
    (*this->localDistributions)(D3Q27System::ET_TSW, x1 + 1, x2 + 1, x3) = rdata[index++];

    (*this->nonLocalDistributions)(D3Q27System::ET_W, x1 + 1, x2, x3)           = rdata[index++];
    (*this->nonLocalDistributions)(D3Q27System::ET_S, x1, x2 + 1, x3)           = rdata[index++];
    (*this->nonLocalDistributions)(D3Q27System::ET_B, x1, x2, x3 + 1)           = rdata[index++];
    (*this->nonLocalDistributions)(D3Q27System::ET_SW, x1 + 1, x2 + 1, x3)      = rdata[index++];
    (*this->nonLocalDistributions)(D3Q27System::ET_SE, x1, x2 + 1, x3)          = rdata[index++];
    (*this->nonLocalDistributions)(D3Q27System::ET_BW, x1 + 1, x2, x3 + 1)      = rdata[index++];
    (*this->nonLocalDistributions)(D3Q27System::ET_BE, x1, x2, x3 + 1)          = rdata[index++];
    (*this->nonLocalDistributions)(D3Q27System::ET_BS, x1, x2 + 1, x3 + 1)      = rdata[index++];
    (*this->nonLocalDistributions)(D3Q27System::ET_BN, x1, x2, x3 + 1)          = rdata[index++];
    (*this->nonLocalDistributions)(D3Q27System::ET_BSW, x1 + 1, x2 + 1, x3 + 1) = rdata[index++];
    (*this->nonLocalDistributions)(D3Q27System::ET_BSE, x1, x2 + 1, x3 + 1)     = rdata[index++];
    (*this->nonLocalDistributions)(D3Q27System::ET_BNW, x1 + 1, x2, x3 + 1)     = rdata[index++];
    (*this->nonLocalDistributions)(D3Q27System::ET_BNE, x1, x2, x3 + 1)         = rdata[index++];

    (*this->zeroDistributions)(x1, x2, x3) = rdata[index++];
}

#endif // OneDistributionFullVectorConnector_H

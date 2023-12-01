#ifndef OneDistributionFullVectorConnector_H
#define OneDistributionFullVectorConnector_H

#include <vector>

#include "Block3D.h"
#include "D3Q27System.h"
#include "EsoTwist3D.h"
//#include "EsoTwistD3Q27System.h"
#include "LBMKernel.h"
#include "FullVectorConnector.h"
#include "EsoSplit.h"
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

    void fillData() override
    {
        FullVectorConnector::fillData();
    }

    void distributeData() override
    {
        FullVectorConnector::distributeData();
    }

protected:
    inline void updatePointers() override;
    inline void fillData(vector_type &sdata, int &index, int x1, int x2, int x3) override;
    inline void distributeData(vector_type &rdata, int &index, int x1, int x2, int x3) override;

private:
    CbArray4D<real, IndexerX4X3X2X1>::CbArray4DPtr localDistributions;
    CbArray4D<real, IndexerX4X3X2X1>::CbArray4DPtr nonLocalDistributions;
    CbArray3D<real, IndexerX3X2X1>::CbArray3DPtr zeroDistributions;

    SPtr<EsoTwist3D> fDis;
};
//////////////////////////////////////////////////////////////////////////
inline void OneDistributionFullVectorConnector::updatePointers()
{
    localDistributions    = dynamicPointerCast<EsoSplit>(this->fDis)->getLocalDistributions();
    nonLocalDistributions = dynamicPointerCast<EsoSplit>(this->fDis)->getNonLocalDistributions();
    zeroDistributions     = dynamicPointerCast<EsoSplit>(this->fDis)->getZeroDistributions();
}
//////////////////////////////////////////////////////////////////////////
inline void OneDistributionFullVectorConnector::fillData(vector_type &sdata, int &index, int x1, int x2, int x3)
{
    using namespace vf::lbm::dir;

    sdata[index++] = (*this->localDistributions)(eP00, x1, x2, x3);
    sdata[index++] = (*this->localDistributions)(e0P0, x1, x2, x3);
    sdata[index++] = (*this->localDistributions)(e00P, x1, x2, x3);
    sdata[index++] = (*this->localDistributions)(ePP0, x1, x2, x3);
    sdata[index++] = (*this->localDistributions)(eMP0, x1 + 1, x2, x3);
    sdata[index++] = (*this->localDistributions)(eP0P, x1, x2, x3);
    sdata[index++] = (*this->localDistributions)(eM0P, x1 + 1, x2, x3);
    sdata[index++] = (*this->localDistributions)(e0PP, x1, x2, x3);
    sdata[index++] = (*this->localDistributions)(e0MP, x1, x2 + 1, x3);
    sdata[index++] = (*this->localDistributions)(ePPP, x1, x2, x3);
    sdata[index++] = (*this->localDistributions)(eMPP, x1 + 1, x2, x3);
    sdata[index++] = (*this->localDistributions)(ePMP, x1, x2 + 1, x3);
    sdata[index++] = (*this->localDistributions)(eMMP, x1 + 1, x2 + 1, x3);

    sdata[index++] = (*this->nonLocalDistributions)(eM00, x1 + 1, x2, x3);
    sdata[index++] = (*this->nonLocalDistributions)(e0M0, x1, x2 + 1, x3);
    sdata[index++] = (*this->nonLocalDistributions)(e00M, x1, x2, x3 + 1);
    sdata[index++] = (*this->nonLocalDistributions)(eMM0, x1 + 1, x2 + 1, x3);
    sdata[index++] = (*this->nonLocalDistributions)(ePM0, x1, x2 + 1, x3);
    sdata[index++] = (*this->nonLocalDistributions)(eM0M, x1 + 1, x2, x3 + 1);
    sdata[index++] = (*this->nonLocalDistributions)(eP0M, x1, x2, x3 + 1);
    sdata[index++] = (*this->nonLocalDistributions)(e0MM, x1, x2 + 1, x3 + 1);
    sdata[index++] = (*this->nonLocalDistributions)(e0PM, x1, x2, x3 + 1);
    sdata[index++] = (*this->nonLocalDistributions)(eMMM, x1 + 1, x2 + 1, x3 + 1);
    sdata[index++] = (*this->nonLocalDistributions)(ePMM, x1, x2 + 1, x3 + 1);
    sdata[index++] = (*this->nonLocalDistributions)(eMPM, x1 + 1, x2, x3 + 1);
    sdata[index++] = (*this->nonLocalDistributions)(ePPM, x1, x2, x3 + 1);

    sdata[index++] = (*this->zeroDistributions)(x1, x2, x3);
}
//////////////////////////////////////////////////////////////////////////
inline void OneDistributionFullVectorConnector::distributeData(vector_type &rdata, int &index, int x1, int x2, int x3)
{
    using namespace vf::lbm::dir;

    (*this->localDistributions)(eP00, x1, x2, x3)           = rdata[index++];
    (*this->localDistributions)(e0P0, x1, x2, x3)           = rdata[index++];
    (*this->localDistributions)(e00P, x1, x2, x3)           = rdata[index++];
    (*this->localDistributions)(ePP0, x1, x2, x3)          = rdata[index++];
    (*this->localDistributions)(eMP0, x1 + 1, x2, x3)      = rdata[index++];
    (*this->localDistributions)(eP0P, x1, x2, x3)          = rdata[index++];
    (*this->localDistributions)(eM0P, x1 + 1, x2, x3)      = rdata[index++];
    (*this->localDistributions)(e0PP, x1, x2, x3)          = rdata[index++];
    (*this->localDistributions)(e0MP, x1, x2 + 1, x3)      = rdata[index++];
    (*this->localDistributions)(ePPP, x1, x2, x3)         = rdata[index++];
    (*this->localDistributions)(eMPP, x1 + 1, x2, x3)     = rdata[index++];
    (*this->localDistributions)(ePMP, x1, x2 + 1, x3)     = rdata[index++];
    (*this->localDistributions)(eMMP, x1 + 1, x2 + 1, x3) = rdata[index++];

    (*this->nonLocalDistributions)(eM00, x1 + 1, x2, x3)           = rdata[index++];
    (*this->nonLocalDistributions)(e0M0, x1, x2 + 1, x3)           = rdata[index++];
    (*this->nonLocalDistributions)(e00M, x1, x2, x3 + 1)           = rdata[index++];
    (*this->nonLocalDistributions)(eMM0, x1 + 1, x2 + 1, x3)      = rdata[index++];
    (*this->nonLocalDistributions)(ePM0, x1, x2 + 1, x3)          = rdata[index++];
    (*this->nonLocalDistributions)(eM0M, x1 + 1, x2, x3 + 1)      = rdata[index++];
    (*this->nonLocalDistributions)(eP0M, x1, x2, x3 + 1)          = rdata[index++];
    (*this->nonLocalDistributions)(e0MM, x1, x2 + 1, x3 + 1)      = rdata[index++];
    (*this->nonLocalDistributions)(e0PM, x1, x2, x3 + 1)          = rdata[index++];
    (*this->nonLocalDistributions)(eMMM, x1 + 1, x2 + 1, x3 + 1) = rdata[index++];
    (*this->nonLocalDistributions)(ePMM, x1, x2 + 1, x3 + 1)     = rdata[index++];
    (*this->nonLocalDistributions)(eMPM, x1 + 1, x2, x3 + 1)     = rdata[index++];
    (*this->nonLocalDistributions)(ePPM, x1, x2, x3 + 1)         = rdata[index++];

    (*this->zeroDistributions)(x1, x2, x3) = rdata[index++];
}

#endif // OneDistributionFullVectorConnector_H

#include "OneDistributionFullVectorConnector.h"
#include "DataSet3D.h"
#include "LBMKernel.h"
//////////////////////////////////////////////////////////////////////////
OneDistributionFullVectorConnector::OneDistributionFullVectorConnector(SPtr<Block3D> block, VectorTransmitterPtr sender,
                                                       VectorTransmitterPtr receiver, int sendDir)
    : FullVectorConnector(block, sender, receiver, sendDir)
{
    if (!block || !sender || !receiver)
        UB_THROW(UbException(UB_EXARGS, "sender or receiver == NULL!!"));
}
//////////////////////////////////////////////////////////////////////////
void OneDistributionFullVectorConnector::init()
{
    FullVectorConnector::init();
    
    fDis = dynamicPointerCast<EsoTwist3D>(block.lock()->getKernel()->getDataSet()->getFdistributions());

    int anz = 27;
    switch (sendDir) {
        case D3Q27System::DIR_000:
            UB_THROW(UbException(UB_EXARGS, "ZERO not allowed"));
            break;
        case D3Q27System::DIR_P00:
        case D3Q27System::DIR_M00:
            sender->getData().resize(maxX2 * maxX3 * anz, 0.0);
            break;
        case D3Q27System::DIR_0P0:
        case D3Q27System::DIR_0M0:
            sender->getData().resize(maxX1 * maxX3 * anz, 0.0);
            break;
        case D3Q27System::DIR_00P:
        case D3Q27System::DIR_00M:
            sender->getData().resize(maxX1 * maxX2 * anz, 0.0);
            break;

        case D3Q27System::DIR_PP0:
        case D3Q27System::DIR_MM0:
        case D3Q27System::DIR_PM0:
        case D3Q27System::DIR_MP0:
            sender->getData().resize(maxX3 * anz, 0.0);
            break;

        case D3Q27System::DIR_P0P:
        case D3Q27System::DIR_M0M:
        case D3Q27System::DIR_P0M:
        case D3Q27System::DIR_M0P:
            sender->getData().resize(maxX2 * anz, 0.0);
            break;

        case D3Q27System::DIR_0PP:
        case D3Q27System::DIR_0MM:
        case D3Q27System::DIR_0PM:
        case D3Q27System::DIR_0MP:
            sender->getData().resize(maxX1 * anz, 0.0);
            break;

        case D3Q27System::DIR_PPP:
        case D3Q27System::DIR_MMM:
        case D3Q27System::DIR_PPM:
        case D3Q27System::DIR_MMP:
        case D3Q27System::DIR_PMP:
        case D3Q27System::DIR_MPM:
        case D3Q27System::DIR_PMM:
        case D3Q27System::DIR_MPP:
            sender->getData().resize(anz, 0.0);
            break;

        default:
            UB_THROW(UbException(UB_EXARGS, "unknown sendDir"));
    }
}
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

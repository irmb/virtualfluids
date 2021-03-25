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
    fDis = dynamicPointerCast<EsoTwist3D>(block.lock()->getKernel()->getDataSet()->getFdistributions());

    int anz = 27;
    switch (sendDir) {
        case D3Q27System::REST:
            UB_THROW(UbException(UB_EXARGS, "ZERO not allowed"));
            break;
        case D3Q27System::E:
        case D3Q27System::W:
            sender->getData().resize(maxX2 * maxX3 * anz, 0.0);
            break;
        case D3Q27System::N:
        case D3Q27System::S:
            sender->getData().resize(maxX1 * maxX3 * anz, 0.0);
            break;
        case D3Q27System::T:
        case D3Q27System::B:
            sender->getData().resize(maxX1 * maxX2 * anz, 0.0);
            break;

        case D3Q27System::NE:
        case D3Q27System::SW:
        case D3Q27System::SE:
        case D3Q27System::NW:
            sender->getData().resize(maxX3 * anz, 0.0);
            break;

        case D3Q27System::TE:
        case D3Q27System::BW:
        case D3Q27System::BE:
        case D3Q27System::TW:
            sender->getData().resize(maxX2 * anz, 0.0);
            break;

        case D3Q27System::TN:
        case D3Q27System::BS:
        case D3Q27System::BN:
        case D3Q27System::TS:
            sender->getData().resize(maxX1 * anz, 0.0);
            break;

        case D3Q27System::TNE:
        case D3Q27System::BSW:
        case D3Q27System::BNE:
        case D3Q27System::TSW:
        case D3Q27System::TSE:
        case D3Q27System::BNW:
        case D3Q27System::BSE:
        case D3Q27System::TNW:
            sender->getData().resize(anz, 0.0);
            break;

        default:
            UB_THROW(UbException(UB_EXARGS, "unknown sendDir"));
    }
}
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

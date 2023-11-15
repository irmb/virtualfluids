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
    using namespace vf::lbm::dir;
    using namespace vf::basics::constant;

    FullVectorConnector::init();
    
    fDis = dynamicPointerCast<EsoTwist3D>(block.lock()->getKernel()->getDataSet()->getFdistributions());

    int anz = 27;
    switch (sendDir) {
        case d000:
            UB_THROW(UbException(UB_EXARGS, "ZERO not allowed"));
            break;
        case dP00:
        case dM00:
            sender->getData().resize(maxX2 * maxX3 * anz, c0o1);
            break;
        case d0P0:
        case d0M0:
            sender->getData().resize(maxX1 * maxX3 * anz, c0o1);
            break;
        case d00P:
        case d00M:
            sender->getData().resize(maxX1 * maxX2 * anz, c0o1);
            break;

        case dPP0:
        case dMM0:
        case dPM0:
        case dMP0:
            sender->getData().resize(maxX3 * anz, c0o1);
            break;

        case dP0P:
        case dM0M:
        case dP0M:
        case dM0P:
            sender->getData().resize(maxX2 * anz, c0o1);
            break;

        case d0PP:
        case d0MM:
        case d0PM:
        case d0MP:
            sender->getData().resize(maxX1 * anz, c0o1);
            break;

        case dPPP:
        case dMMM:
        case dPPM:
        case dMMP:
        case dPMP:
        case dMPM:
        case dPMM:
        case dMPP:
            sender->getData().resize(anz, c0o1);
            break;

        default:
            UB_THROW(UbException(UB_EXARGS, "unknown sendDir"));
    }
}
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

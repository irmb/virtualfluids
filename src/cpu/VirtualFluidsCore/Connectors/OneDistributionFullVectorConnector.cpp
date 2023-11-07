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
        case DIR_0P0:
        case DIR_0M0:
            sender->getData().resize(maxX1 * maxX3 * anz, c0o1);
            break;
        case DIR_00P:
        case DIR_00M:
            sender->getData().resize(maxX1 * maxX2 * anz, c0o1);
            break;

        case DIR_PP0:
        case DIR_MM0:
        case DIR_PM0:
        case DIR_MP0:
            sender->getData().resize(maxX3 * anz, c0o1);
            break;

        case DIR_P0P:
        case DIR_M0M:
        case DIR_P0M:
        case DIR_M0P:
            sender->getData().resize(maxX2 * anz, c0o1);
            break;

        case DIR_0PP:
        case DIR_0MM:
        case DIR_0PM:
        case DIR_0MP:
            sender->getData().resize(maxX1 * anz, c0o1);
            break;

        case DIR_PPP:
        case DIR_MMM:
        case DIR_PPM:
        case DIR_MMP:
        case DIR_PMP:
        case DIR_MPM:
        case DIR_PMM:
        case DIR_MPP:
            sender->getData().resize(anz, c0o1);
            break;

        default:
            UB_THROW(UbException(UB_EXARGS, "unknown sendDir"));
    }
}
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

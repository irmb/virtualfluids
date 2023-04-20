#include "VoidLBMKernel.h"
#include "BCSet.h"
#include "DataSet3D.h"
#include "VoidData3D.h"
#include "D3Q27System.h"

VoidLBMKernel::VoidLBMKernel() = default;
//////////////////////////////////////////////////////////////////////////
VoidLBMKernel::~VoidLBMKernel() = default;
//////////////////////////////////////////////////////////////////////////
void VoidLBMKernel::initDataSet()
{
    SPtr<DistributionArray3D> d(new VoidData3D(nx[0] + 2, nx[1] + 2, nx[2] + 2, -999.9));
    dataSet->setFdistributions(d);
}
//////////////////////////////////////////////////////////////////////////
SPtr<LBMKernel> VoidLBMKernel::clone()
{
    SPtr<LBMKernel> kernel(new VoidLBMKernel());
    kernel->setNX(nx);
    dynamicPointerCast<VoidLBMKernel>(kernel)->initDataSet();
    kernel->setCollisionFactor(this->collFactor);
    kernel->setBCSet(bcSet->clone(kernel));
    kernel->setWithForcing(withForcing);
    kernel->setForcingX1(muForcingX1);
    kernel->setForcingX2(muForcingX2);
    kernel->setForcingX3(muForcingX3);
    kernel->setIndex(ix1, ix2, ix3);
    kernel->setDeltaT(deltaT);
    return kernel;
}
//////////////////////////////////////////////////////////////////////////
void VoidLBMKernel::calculate(int step) {}
//////////////////////////////////////////////////////////////////////////
real VoidLBMKernel::getCalculationTime() { return vf::lbm::constant::c0o1; }

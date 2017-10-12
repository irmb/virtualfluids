#include "VoidLBMKernel.h"
#include "VoidData3D.h"

VoidLBMKernel::VoidLBMKernel()
{

}

VoidLBMKernel::VoidLBMKernel(int nx1, int nx2, int nx3) : nx1(nx1), nx2(nx2), nx3(nx3)
{
   DistributionArray3DPtr d(new VoidData3D(nx1+2, nx2+2, nx3+2, -999.0));
   dataSet->setFdistributions(d);
}

VoidLBMKernel::~VoidLBMKernel()
{

}

LBMKernelPtr VoidLBMKernel::clone()
{
   LBMKernelPtr kernel(new VoidLBMKernel(nx1, nx2, nx3));
   kernel->setCollisionFactor(this->collFactor);
   kernel->setBCProcessor(bcProcessor->clone(kernel));
   kernel->setWithForcing(withForcing);
   kernel->setForcingX1(muForcingX1);
   kernel->setForcingX2(muForcingX2);
   kernel->setForcingX3(muForcingX3);
   kernel->setIndex(ix1, ix2, ix3);
   kernel->setDeltaT(deltaT);
   return kernel;
}

void VoidLBMKernel::calculate()
{

}

void VoidLBMKernel::swapDistributions()
{

}

double VoidLBMKernel::getCallculationTime()
{
   return 0.0;
}

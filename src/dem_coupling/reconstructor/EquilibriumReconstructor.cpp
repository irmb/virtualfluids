#include "EquilibriumReconstructor.h"

#include "ILBMKernel.h"
#include "D3Q27System.h"
#include "DataSet3D.h"
#include "BCProcessor.h"
#include "BCArray3D.h"

#include <dem_coupling/physicsEngineAdapter/PhysicsEngineGeometryAdapter.h>

void EquilibriumReconstructor::reconstructNode(const int& x1, const int& x2, const int& x3,
                                               const Vector3D& worldCoordinates, std::shared_ptr<PhysicsEngineGeometryAdapter> physicsEngineGeometry, std::shared_ptr<ILBMKernel> kernel) const
{
    const double averageDensity = this->getLocalAverageDensity(x1, x2, x3, kernel);
    LBMReal feq[27];
    const Vector3D boundaryVelocity = physicsEngineGeometry->getVelocityAtPosition(worldCoordinates);

    if (kernel->getCompressible())
        D3Q27System::calcCompFeq(feq, averageDensity, boundaryVelocity[0], boundaryVelocity[1], boundaryVelocity[2]);
    else
        D3Q27System::calcIncompFeq(feq, averageDensity, boundaryVelocity[0], boundaryVelocity[1], boundaryVelocity[2]);


    DistributionArray3DPtr distributions = kernel->getDataSet()->getFdistributions();
    distributions->setDistribution(feq, x1, x2, x3);
    distributions->setDistributionInv(feq, x1, x2, x3);
}


double EquilibriumReconstructor::getLocalAverageDensity(const int &x1, const int &x2, const int &x3, std::shared_ptr<ILBMKernel> kernel) const
{
    int nAverage = 0;
    double averageDensity = 0.0;

    BCArray3DPtr bcArray = kernel->getBCProcessor()->getBCArray();

    LBMReal f[D3Q27System::ENDF + 1];
    DistributionArray3DPtr distributions = kernel->getDataSet()->getFdistributions();

    int neighborX1, neighborX2, neighborX3;
    for (int fDir = D3Q27System::FSTARTDIR; fDir <= D3Q27System::FENDDIR; fDir++)
    {
        neighborX1 = x1 + D3Q27System::DX1[fDir];
        neighborX2 = x2 + D3Q27System::DX2[fDir];
        neighborX3 = x3 + D3Q27System::DX3[fDir];

        if (bcArray->isFluid(neighborX1, neighborX2, neighborX3))
        {
            distributions->getDistribution(f, neighborX1, neighborX2, neighborX3);
            averageDensity += D3Q27System::getDensity(f);
            ++nAverage;
        }
    }
    return (nAverage > 0) ? averageDensity / nAverage : 0.0;
}


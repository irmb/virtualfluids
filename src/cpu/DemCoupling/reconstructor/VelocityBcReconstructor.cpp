#include "VelocityBcReconstructor.h"

#include <exception>

#include "BCArray3D.h"
#include "BCProcessor.h"
#include "D3Q27System.h"
#include "DataSet3D.h"
#include "EsoTwist3D.h"
#include "ILBMKernel.h"

#include "PhysicsEngineGeometryAdapter.h"

void VelocityBcReconstructor::reconstructNode(const int &x1, const int &x2, const int &x3,
                                              const Vector3D &worldCoordinates,
                                              std::shared_ptr<PhysicsEngineGeometryAdapter> physicsEngineGeometry,
                                              std::shared_ptr<ILBMKernel> kernel) const
{
    if (kernel->getCompressible())
        throw std::runtime_error("not implemented yet!");

    const Vector3D boundaryVelocity = physicsEngineGeometry->getVelocityAtPosition(worldCoordinates);
    // TODO: move to D3Q27 system
    LBMReal wijk[D3Q27System::ENDF + 1];
    D3Q27System::calcIncompFeq(wijk, 1, 0, 0, 0);

    SPtr<BCArray3D> bcArray = kernel->getBCProcessor()->getBCArray();

    SPtr<BoundaryConditions> bc = SPtr<BoundaryConditions>(new BoundaryConditions());
    bc->setBoundaryVelocityX1((float)boundaryVelocity[0]);
    bc->setBoundaryVelocityX2((float)boundaryVelocity[1]);
    bc->setBoundaryVelocityX3((float)boundaryVelocity[2]);

    LBMReal feqNullRho[D3Q27System::ENDF + 1];
    D3Q27System::calcIncompFeq(feqNullRho, 0, boundaryVelocity[0], boundaryVelocity[1], boundaryVelocity[2]);

    LBMReal fpre[D3Q27System::ENDF + 1];
    LBMReal fpost[D3Q27System::ENDF + 1];
    SPtr<DistributionArray3D> distributions = kernel->getDataSet()->getFdistributions();

    distributions->swap();
    distributions->getDistributionInv(fpost, x1, x2, x3);
    distributions->swap();
    distributions->getDistribution(fpre, x1, x2, x3);

    int neighborX1, neighborX2, neighborX3;
    int neighborX1Inv, neighborX2Inv, neighborX3Inv;

    double sumRho = 0, sumWijk = 0;
    double collFactor = kernel->getCollisionFactor();

    for (int fDir = D3Q27System::FSTARTDIR; fDir <= D3Q27System::FENDDIR; fDir++) {
        neighborX1 = x1 + D3Q27System::DX1[fDir];
        neighborX2 = x2 + D3Q27System::DX2[fDir];
        neighborX3 = x3 + D3Q27System::DX3[fDir];

        if (bcArray->isFluid(neighborX1, neighborX2, neighborX3)) {
            int invDir = D3Q27System::INVDIR[fDir];

            neighborX1Inv = x1 + D3Q27System::DX1[invDir];
            neighborX2Inv = x2 + D3Q27System::DX2[invDir];
            neighborX3Inv = x3 + D3Q27System::DX3[invDir];
            if (!bcArray->isFluid(neighborX1Inv, neighborX2Inv, neighborX3Inv)) {

                double velocity = bc->getBoundaryVelocity(invDir);

                fpre[fDir]   = fpre[invDir] - velocity;
                double Omega = fpost[fDir] - fpre[fDir];

                sumRho += Omega / collFactor + fpre[fDir] - feqNullRho[fDir];
                sumWijk += wijk[fDir];
            }
        }
    }

    double rho = 0.0;
    if (sumWijk > 0.0)
        rho = sumRho / sumWijk;

    for (int fDir = D3Q27System::FSTARTDIR; fDir <= D3Q27System::FENDDIR; fDir++) {
        neighborX1 = x1 + D3Q27System::DX1[fDir];
        neighborX2 = x2 + D3Q27System::DX2[fDir];
        neighborX3 = x3 + D3Q27System::DX3[fDir];

        if (!bcArray->isFluid(neighborX1, neighborX2, neighborX3)) {
            int invDir    = D3Q27System::INVDIR[fDir];
            neighborX1Inv = x1 + D3Q27System::DX1[invDir];
            neighborX2Inv = x2 + D3Q27System::DX2[invDir];
            neighborX3Inv = x3 + D3Q27System::DX3[invDir];
            if (!bcArray->isFluid(neighborX1Inv, neighborX2Inv, neighborX3Inv)) {
                fpre[fDir] = D3Q27System::getIncompFeqForDirection(
                    fDir, rho, bc->getBoundaryVelocityX1(), bc->getBoundaryVelocityX2(), bc->getBoundaryVelocityX3());
            }
        }
    }
}

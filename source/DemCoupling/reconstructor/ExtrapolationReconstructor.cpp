#include "ExtrapolationReconstructor.h"

#include "ILBMKernel.h"
#include "D3Q27System.h"
#include "DataSet3D.h"
#include "BCProcessor.h"
#include "BCArray3D.h"

#include "PhysicsEngineGeometryAdapter.h"
#include "DistributionArray3D.h"

void ExtrapolationReconstructor::setAlternativeReconstructor(std::shared_ptr<Reconstructor> alternativeReconstructor)
{
    this->alternativeReconstructor = alternativeReconstructor;
}


ExtrapolationReconstructor::ExtrapolationReconstructor(std::shared_ptr<Reconstructor> alternativeReconstructor) : alternativeReconstructor(alternativeReconstructor)
{
}

void ExtrapolationReconstructor::reconstructNode(const int& x1, const int& x2, const int& x3,
                                                 const Vector3D& worldCoordinates, std::shared_ptr<PhysicsEngineGeometryAdapter> physicsEngineGeometry, std::shared_ptr<ILBMKernel> kernel) const
{
    const UbTupleInt3 extrapolationDirection = getSphereDirection(worldCoordinates, physicsEngineGeometry);
    const int numberOfCellsForExtrapolation = getNumberOfExtrapolationCells(x1, x2, x3, extrapolationDirection, kernel);

    if (numberOfCellsForExtrapolation < 2)
        alternativeReconstructor->reconstructNode(x1, x2, x3, worldCoordinates, physicsEngineGeometry, kernel);
    else
    {
        //UBLOG(logINFO, "point (x,y,z) " << val<1>(worldCoordinates) << ", " << val<2>(worldCoordinates) << ", " << val<3>(worldCoordinates));
        //UBLOG(logINFO, "extradir (x,y,z) " << val<1>(extrapolationDirection) << ", " << val<2>(extrapolationDirection) << ", " << val<3>(extrapolationDirection));
        //UBLOG(logINFO, "numberOfCellsForExtrapolation: " << numberOfCellsForExtrapolation );

        this->extrapolatePdFs(x1, x2, x3, extrapolationDirection, numberOfCellsForExtrapolation, kernel);
    }
        
}

UbTupleInt3 ExtrapolationReconstructor::getSphereDirection(const Vector3D& worldCoordinates, std::shared_ptr<PhysicsEngineGeometryAdapter> physicsEngineGeometry) const
{
    const Vector3D spherePosition = physicsEngineGeometry->getPosition();
    const Vector3D bodyNormal = worldCoordinates - spherePosition;
    return this->getCorrespondingLatticeDirection(bodyNormal);
}

UbTupleInt3 ExtrapolationReconstructor::getCorrespondingLatticeDirection(const Vector3D& direction) const
{
    int correspondingDirection = 0;
    double innerProduct = 0.0;
    for (int fDir = D3Q27System::FSTARTDIR; fDir <= D3Q27System::FENDDIR; fDir++)
    {
        // compute inner product <dir,c_i>
        const double temporaryInnerProduct = direction[0] * D3Q27System::cNorm[0][fDir] + direction[1] * D3Q27System::cNorm[1][fDir] + direction[2] * D3Q27System::cNorm[2][fDir];
        if (temporaryInnerProduct > innerProduct)
        {
            innerProduct = temporaryInnerProduct;
            correspondingDirection = fDir;
        }
    }

    return UbTupleInt3(D3Q27System::DX1[correspondingDirection], D3Q27System::DX2[correspondingDirection], D3Q27System::DX3[correspondingDirection]);
}

int ExtrapolationReconstructor::getNumberOfExtrapolationCells(const int x1, const int x2, const int x3, const UbTupleInt3& extrapolationDirection, std::shared_ptr<ILBMKernel> kernel) const
{
    if (extrapolationDirection == UbTupleInt3(0, 0, 0))
        return 0;

    const int desiredCellsInExtrapolationDirection = 3;
   
    for (int numCells = 1; numCells <= desiredCellsInExtrapolationDirection; ++numCells)
    {
        UbTupleInt3 neighbor(x1 + numCells * val<1>(extrapolationDirection), x2 + numCells * val<2>(extrapolationDirection), x3 + numCells * val<3>(extrapolationDirection));

        if(!kernel->isInsideOfDomain(val<1>(neighbor), val<2>(neighbor), val<3>(neighbor)))
            return numCells - 1;


        if (!kernel->getBCProcessor()->getBCArray()->isFluid(val<1>(neighbor), val<2>(neighbor), val<3>(neighbor)))
            return numCells - 1;
    }
    return desiredCellsInExtrapolationDirection;
}


void ExtrapolationReconstructor::extrapolatePdFs(const int x1, const int x2, const int x3,
    const UbTupleInt3& extrapolationDirection, int numberOfCellsForExtrapolation, std::shared_ptr<ILBMKernel> kernel) const
{
    SPtr<DistributionArray3D> distributions = kernel->getDataSet()->getFdistributions();

    const int nx1 = val<1>(extrapolationDirection);
    const int nx2 = val<2>(extrapolationDirection);
    const int nx3 = val<3>(extrapolationDirection);

    LBMReal pdf[D3Q27System::ENDF + 1];
    LBMReal pdfNeighbor1[D3Q27System::ENDF + 1];
    LBMReal pdfNeighbor2[D3Q27System::ENDF + 1];

    distributions->getDistribution(pdf, x1, x2, x3);
    distributions->getDistribution(pdfNeighbor1, x1 + nx1, x2 + nx2, x3 + nx3);
    distributions->getDistribution(pdfNeighbor2, x1 + 2 * nx1, x2 + 2 * nx2, x3 + 2 * nx3);

    if (numberOfCellsForExtrapolation == 3) // quadratic normal extrapolation
    {
        LBMReal pdfNeighbor3[D3Q27System::ENDF + 1];
        distributions->getDistribution(pdfNeighbor3, x1 + 3 * nx1, x2 + 3 * nx2, x3 + 3 * nx3);

        for (int fDir = D3Q27System::FSTARTDIR; fDir <= D3Q27System::FENDDIR; fDir++)
            pdf[fDir] = 3 * pdfNeighbor1[fDir] - 3 * pdfNeighbor2[fDir] + pdfNeighbor3[fDir];
    }
    else  // numberOfCellsForExtrapolation == 2 // linear normal extrapolation
    { 
        for (int fDir = D3Q27System::FSTARTDIR; fDir <= D3Q27System::FENDDIR; fDir++)
            pdf[fDir] = 2 * pdfNeighbor1[fDir] - pdfNeighbor2[fDir];
    }

    distributions->setDistribution(pdf, x1, x2, x3);
}

#ifndef INTERPOLATIONPROCESSOR_H
#define INTERPOLATIONPROCESSOR_H

#include <memory>

#include "BCArray3D.h"
#include "BoundaryConditions.h"
#include "DistributionArray3D.h"
#include "LBMSystem.h"

struct D3Q27ICell {
    real TSW[27];
    real TNW[27];
    real TNE[27];
    real TSE[27];
    real BSW[27];
    real BNW[27];
    real BNE[27];
    real BSE[27];
};

class InterpolationProcessor;
using InterpolationProcessorPtr = std::shared_ptr<InterpolationProcessor>;


class InterpolationProcessor
{
public:
    InterpolationProcessor();
    virtual ~InterpolationProcessor();
    virtual InterpolationProcessorPtr clone()                                    = 0;
    virtual void setOmegas(real omegaC, real omegaF)                       = 0;
    virtual void interpolateCoarseToFine(D3Q27ICell &icellC, D3Q27ICell &icellF) = 0;
    virtual void interpolateCoarseToFine(D3Q27ICell &icellC, D3Q27ICell &icellF, real xoff, real yoff,
                                         real zoff)                           = 0;
    virtual void interpolateFineToCoarse(D3Q27ICell &icellF, real *icellC)    = 0;
    virtual void interpolateFineToCoarse(D3Q27ICell &icellF, real *icellC, real xoff, real yoff,
                                         real zoff)                           = 0;

    static void readICell(SPtr<DistributionArray3D> f, D3Q27ICell &icell, int x1, int x2, int x3);
    static void writeICell(SPtr<DistributionArray3D> f, const D3Q27ICell &icell, int x1, int x2, int x3);
    static void writeICellInv(SPtr<DistributionArray3D> f, const D3Q27ICell &icell, int x1, int x2, int x3);
    static void writeINode(SPtr<DistributionArray3D> f, const real *const inode, int x1, int x2, int x3);
    static void writeINodeInv(SPtr<DistributionArray3D> f, const real *const inode, int x1, int x2, int x3);
    static bool iCellHasSolid(const SPtr<BCArray3D> bcArray, int x1, int x2, int x3);
    static int iCellHowManySolids(const SPtr<BCArray3D> bcArray, int x1, int x2, int x3);

    bool findNeighborICell(const SPtr<BCArray3D> bcArray, SPtr<DistributionArray3D> f, D3Q27ICell &icell, int maxX1,
                           int maxX2, int maxX3, int x1, int x2, int x3, real &xoff, real &yoff, real &zoff);

protected:
    virtual void calcInterpolatedCoefficiets(const D3Q27ICell &icell, real omega, real eps_new) {}
    virtual void calcInterpolatedNodeFC(real *f, real omega) {}
    virtual void calcInterpolatedVelocity(real x, real y, real z, real &vx1, real &vx2, real &vx3) {}
    virtual void calcInterpolatedShearStress(real x, real y, real z, real &tauxx, real &tauyy,
                                             real &tauzz, real &tauxy, real &tauxz, real &tauyz)
    {
    }
    virtual void setOffsets(real xoff, real yoff, real zoff) {}
    friend class InterpolationHelper;

private:
    bool inRange(int x1, int x2, int x3);
    int m_maxX1, m_maxX2, m_maxX3;
};

//////////////////////////////////////////////////////////////////////////
inline bool InterpolationProcessor::inRange(int x1, int x2, int x3)
{
    return x1 >= 0 && x1 < m_maxX1 && x2 >= 0 && x2 < m_maxX2 && x3 >= 0 && x3 < m_maxX3;
}

#endif

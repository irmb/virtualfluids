#ifndef CompressibleOffsetMomentsInterpolationProcessor_H_
#define CompressibleOffsetMomentsInterpolationProcessor_H_

#include "Interpolator.h"
#include "D3Q27System.h"

//////////////////////////////////////////////////////////////////////////
//it works only for cascaded LBM
//super compact interpolation method by Martin Geier
//////////////////////////////////////////////////////////////////////////

class CompressibleOffsetMomentsInterpolator : public Interpolator
{
public:
    CompressibleOffsetMomentsInterpolator() = default;
    CompressibleOffsetMomentsInterpolator(real omegaC, real omegaF);
    InterpolationProcessorPtr clone() override;

    void setOmegas(real omegaC, real omegaF) override;

    void interpolateCoarseToFine(D3Q27ICell &icellC, D3Q27ICell &icellF) override;
    void interpolateCoarseToFine(D3Q27ICell &icellC, D3Q27ICell &icellF, real xoff, real yoff, real zoff) override;
    void interpolateFineToCoarse(D3Q27ICell &icellF, real *icellC) override;
    void interpolateFineToCoarse(D3Q27ICell &icellF, real *icellC, real xoff, real yoff, real zoff) override;

private:
    real omegaC{ 0.0 }, omegaF{ 0.0 };
};

//////////////////////////////////////////////////////////////////////////
inline void CompressibleOffsetMomentsInterpolator::interpolateCoarseToFine(D3Q27ICell &icellC, D3Q27ICell &icellF)
{
    this->interpolateCoarseToFine(icellC, icellF, 0.0, 0.0, 0.0);
}
//////////////////////////////////////////////////////////////////////////
inline void CompressibleOffsetMomentsInterpolator::interpolateFineToCoarse(D3Q27ICell &icellF, real *icellC)
{
    this->interpolateFineToCoarse(icellF, icellC, 0.0, 0.0, 0.0);
}

#endif

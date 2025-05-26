#ifndef LBM_ADVECTION_DIFFUSION_BOUNDARY_CONDITIONS_H
#define LBM_ADVECTION_DIFFUSION_BOUNDARY_CONDITIONS_H
#include <basics/DataTypes.h>
#include <basics/constants/NumericConstants.h>
#include <lbm/advectionDiffusion/Equilibrium.h>
#include <lbm/constants/D3Q27.h>

namespace vf::lbm::advection_diffusion
{

template <size_t direction>
constexpr real computeInterpolatedDistribution(const real* populations, const real concentration, const real velocityX,
                                               const real velocityY, const real velocityZ, const real subgridDistance,
                                               const real relaxationFrequency)
{
    using namespace vf::basics::constant;
    using namespace vf::lbm::dir;

    const size_t inverseDirection = inverseDir<direction>();
    const real population = populations[direction];
    const real populationInverseDirection = populations[inverseDirection];
    const real equilibrium = computeEquilibrium<direction>(concentration, velocityX, velocityY, velocityZ);
    return ((c1o1 - subgridDistance) * ((population - equilibrium * relaxationFrequency) / (c1o1 - relaxationFrequency)) +
            subgridDistance * (population + populationInverseDirection)) /
           (subgridDistance + c1o1);
}

template <size_t direction>
constexpr real computePopulationSimpleBounceBackWithFlux(const real* populations, const real fluxX, const real fluxY,
                                                         const real fluxZ)
{
    using namespace vf::basics::constant;
    using namespace vf::lbm::dir;

    return populations[direction] - c6o1 * getWeight<direction>() * getVelocity<direction>(fluxX, fluxY, fluxZ);
}

template <size_t direction>
constexpr real computePopulationInterpolatedBounceBackWithFlux(const real subgridDistance, const real* populations,
                                                               const real vx1, const real vx2, const real vx3,
                                                               const real relaxationFrequency, const real concentration,
                                                               const real fluxX, const real fluxY, const real fluxZ)
{
    using namespace vf::basics::constant;
    using namespace vf::lbm::dir;

    const real flux = getVelocity<direction>(fluxX, fluxY, fluxZ);
    const real interpolated = computeInterpolatedDistribution<direction>(populations, concentration, vx1, vx2, vx3,
                                                                         subgridDistance, relaxationFrequency);
    return interpolated - c6o1 * getWeight<direction>() * flux / (subgridDistance + c1o1);
}

template <size_t direction>
constexpr real computePopulationSimpleAntiBounceBack(const real* populations, const real concentrationWall,
                                                     const real velocityWallX, const real velocityWallY,
                                                     const real velocityWallZ)
{
    using namespace vf::basics::constant;

    const real equilibriumWall =
        computeEquilibrium<direction>(concentrationWall, velocityWallX, velocityWallY, velocityWallZ);

    return -populations[direction] + c2o1 * equilibriumWall;
}

template <size_t direction>
constexpr real computePopulationInterpolatedAntiBounceBack(const real subgridDistance, const real* populations,
                                                           const real concentrationNode, const real concentrationWall,
                                                           const real vx1, real vx2, const real vx3,
                                                           const real velocityWallX, const real velocityWallY,
                                                           const real velocityWallZ, const real relaxationFrequency)
{
    using namespace vf::basics::constant;

    const real equilibriumWall =
        computeEquilibrium<direction>(concentrationWall, velocityWallX, velocityWallY, velocityWallZ);
    const real interpolated = computeInterpolatedDistribution<direction>(populations, concentrationNode, vx1, vx2, vx3,
                                                                         subgridDistance, relaxationFrequency);
    return -interpolated + c2o1 * equilibriumWall / (subgridDistance + c1o1);
}

} // namespace vf::lbm::advection_diffusion

#endif
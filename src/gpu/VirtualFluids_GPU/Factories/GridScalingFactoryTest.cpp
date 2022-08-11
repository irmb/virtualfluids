#include <gmock/gmock.h>
#include <typeindex>

#include "GPU/GPU_Interface.h"
#include "GridScalingFactory.h"

using gridScalingFCtarget = void (*)(LBMSimulationParameter *, LBMSimulationParameter *, ICellFC *,
                                     CUstream_st *stream);
using gridScalingCFtarget = void (*)(LBMSimulationParameter *, LBMSimulationParameter *, ICellCF *, OffCF,
                                     CUstream_st *stream);

// tests for default scaling functions
TEST(GridScalingFactoryTest, defaultGridScalingFC)
{
    auto scalingFactory = GridScalingFactory();
    auto scalingFunction = scalingFactory.getGridScalingFC();
    EXPECT_THAT(scalingFunction, testing::Eq(nullptr));
    EXPECT_THROW(scalingFunction(nullptr, nullptr, nullptr, nullptr), std::bad_function_call);
}

// TEST(GridScalingFactoryTest, defaultGridScalingCF)
// {
//     auto scalingFactory = GridScalingFactory();
//     auto scalingFunction = scalingFactory.getGridScalingCF();
//     EXPECT_THAT(scalingFunction, testing::Eq(nullptr));
//     EXPECT_THROW(scalingFunction(nullptr, nullptr, nullptr, OffCF(), nullptr), std::bad_function_call);
// }

// tests for boundary conditions which are set by the user (tests both set and get functions)
gridScalingFCtarget getScalingFCTarget(GridScalingFactory &scalingFactory)
{
    auto scaling = scalingFactory.getGridScalingFC();
    void (*scalingTarget)(LBMSimulationParameter *, LBMSimulationParameter *, ICellFC *, CUstream_st * stream) =
        (*scaling
              .target<void (*)(LBMSimulationParameter *, LBMSimulationParameter *, ICellFC *, CUstream_st * stream)>());
    return scalingTarget;
}

TEST(GridScalingFactoryTest, gridScalingFC)
{
    auto scalingFactory = GridScalingFactory();

    scalingFactory.setScalingFactory(GridScalingFactory::GridScaling::ScaleK17);
    EXPECT_TRUE(*(getScalingFCTarget(scalingFactory)) == ScaleFC_K17_redesigned)
        << "The returned boundary condition is not the expected function ScaleFC_K17_redesigned.";

    scalingFactory.setScalingFactory(GridScalingFactory::GridScaling::ScaleRhoSq);
    EXPECT_TRUE(*(getScalingFCTarget(scalingFactory)) == ScaleFC_RhoSq_comp_27)
        << "The returned boundary condition is not the expected function ScaleFC_RhoSq_comp_27.";
}

// gridScalingCFtarget getScalingCFTarget(GridScalingFactory &scalingFactory)
// {
//     auto scaling = scalingFactory.getGridScalingCF();
//     void (*scalingTarget)(LBMSimulationParameter *, LBMSimulationParameter *, ICellCF *, OffCF, CUstream_st *stream)
//     =
//         (*scaling.target<void (*)(LBMSimulationParameter *, LBMSimulationParameter *, ICellCF *, OffCF, CUstream_st
//         *stream)>());
//     return scalingTarget;
// }

// TEST(GridScalingFactoryTest, gridScalingCF)
// {
//     auto scalingFactory = GridScalingFactory();

//     scalingFactory.setScalingFactory(GridScalingFactory::GridScaling::ScaleK17);
//     EXPECT_TRUE( *(getScalingCFTarget(scalingFactory)) == ScaleCF_K17_redesigned)
//         << "The returned boundary condition is not the expected function ScaleFC_K17_redesigned.";

//     scalingFactory.setScalingFactory(GridScalingFactory::GridScaling::ScaleRhoSq);
//     EXPECT_TRUE( *(getScalingCFTarget(scalingFactory)) == ScaleCF_RhoSq_comp_27)
//         << "The returned boundary condition is not the expected function ScaleFC_RhoSq_comp_27.";
// }

#include <gmock/gmock.h>
#include <typeindex>

#include "BoundaryConditionFactory.h"
#include "GPU/GPU_Interface.h"
#include "gpu/GridGenerator/grid/BoundaryConditions/BoundaryCondition.h"

using bcFunction = void (*)(LBMSimulationParameter *, QforBoundaryConditions *);
using bcFunctionParamter = void (*)(Parameter *, QforBoundaryConditions *, const int level);

// tests for default boundary conditions
TEST(BoundaryConditionFactoryTest, defaultVelocityBC)
{
    auto bcFactory = BoundaryConditionFactory();
    auto bc = bcFactory.getVelocityBoundaryConditionPost();
    EXPECT_THAT(bc, testing::Eq(nullptr));
    EXPECT_THROW(bc(nullptr, nullptr), std::bad_function_call);
}

TEST(BoundaryConditionFactoryTest, defaultNoSlipBC)
{
    auto bcFactory = BoundaryConditionFactory();
    auto bc = bcFactory.getNoSlipBoundaryConditionPost();
    EXPECT_NO_THROW(bc(nullptr, nullptr)); // empty lambda function should not throw
}

TEST(BoundaryConditionFactoryTest, defaultSlipBC)
{
    auto bcFactory = BoundaryConditionFactory();
    auto bc = bcFactory.getSlipBoundaryConditionPost();
    EXPECT_THAT(bc, testing::Eq(nullptr));
    EXPECT_THROW(bc(nullptr, nullptr), std::bad_function_call);
}

TEST(BoundaryConditionFactoryTest, defaultPressureBC)
{
    auto bcFactory = BoundaryConditionFactory();
    auto bc = bcFactory.getPressureBoundaryConditionPre();
    EXPECT_THAT(bc, testing::Eq(nullptr));
    EXPECT_THROW(bc(nullptr, nullptr), std::bad_function_call);
}

TEST(BoundaryConditionFactoryTest, defaultGeometryBC)
{
    auto bcFactory = BoundaryConditionFactory();
    auto bc = bcFactory.getGeometryBoundaryConditionPost();
    EXPECT_NO_THROW(bc(nullptr, nullptr)); // empty lambda function should not throw
}

TEST(BoundaryConditionFactoryTest, defaultStressBC)
{
    auto bcFactory = BoundaryConditionFactory();
    auto bc = bcFactory.getStressBoundaryConditionPost();
    EXPECT_THAT(bc, testing::Eq(nullptr));
    EXPECT_THROW(bc(nullptr, nullptr, 0), std::bad_function_call);
}

bcFunction getVelocityBcTarget(BoundaryConditionFactory &bcFactory)
{
    auto bc = bcFactory.getVelocityBoundaryConditionPost();
    void (*bcTarget)(LBMSimulationParameter *, QforBoundaryConditions *) =
        (*bc.target<void (*)(LBMSimulationParameter *, QforBoundaryConditions *)>());
    return bcTarget;
}

// tests for boundary conditions whcih are set by the user (tests both set and get functions)

TEST(BoundaryConditionFactoryTest, velocityBC)
{
    auto bcFactory = BoundaryConditionFactory();

    bcFactory.setVelocityBoundaryCondition(BoundaryConditionFactory::VelocityBC::VelocitySimpleBounceBackCompressible);
    EXPECT_THAT(*(getVelocityBcTarget(bcFactory)), testing::Eq(QVelDevicePlainBB27))
        << "The returned boundary condition is not the expected function QVelDevicePlainBB27.";

    bcFactory.setVelocityBoundaryCondition(BoundaryConditionFactory::VelocityBC::VelocityIncompressible);
    EXPECT_THAT(*(getVelocityBcTarget(bcFactory)), testing::Eq(QVelDev27))
        << "The returned boundary condition is not the expected function QVelDev27.";

    bcFactory.setVelocityBoundaryCondition(BoundaryConditionFactory::VelocityBC::VelocityCompressible);
    EXPECT_THAT(*(getVelocityBcTarget(bcFactory)), testing::Eq(QVelDevComp27))
        << "The returned boundary condition is not the expected function QVelDevComp27.";

    bcFactory.setVelocityBoundaryCondition(BoundaryConditionFactory::VelocityBC::VelocityAndPressureCompressible);
    EXPECT_THAT(*(getVelocityBcTarget(bcFactory)), testing::Eq(QVelDevCompZeroPress27))
        << "The returned boundary condition is not the expected function QVelDevCompZeroPress27.";
}

bcFunction getNoSlipBcTarget(BoundaryConditionFactory &bcFactory)
{
    auto bc = bcFactory.getNoSlipBoundaryConditionPost();
    void (*bcTarget)(LBMSimulationParameter *, QforBoundaryConditions *) =
        (*bc.target<void (*)(LBMSimulationParameter *, QforBoundaryConditions *)>());
    return bcTarget;
}

TEST(BoundaryConditionFactoryTest, noSlipBC)
{
    auto bcFactory = BoundaryConditionFactory();

    bcFactory.setNoSlipBoundaryCondition(BoundaryConditionFactory::NoSlipBC::NoSlipImplicitBounceBack);
    auto bc = bcFactory.getNoSlipBoundaryConditionPost();
    EXPECT_NO_THROW(bc(nullptr, nullptr)); // empty lambda function should not throw

    bcFactory.setNoSlipBoundaryCondition(BoundaryConditionFactory::NoSlipBC::NoSlipBounceBack);
    EXPECT_THAT(*(getNoSlipBcTarget(bcFactory)), testing::Eq(BBDev27))
        << "The returned boundary condition is not the expected function BBDev27.";

    bcFactory.setNoSlipBoundaryCondition(BoundaryConditionFactory::NoSlipBC::NoSlipIncompressible);
    EXPECT_THAT(*(getNoSlipBcTarget(bcFactory)), testing::Eq(QDev27))
        << "The returned boundary condition is not the expected function QDev27.";

    bcFactory.setNoSlipBoundaryCondition(BoundaryConditionFactory::NoSlipBC::NoSlipCompressible);
    EXPECT_THAT(*(getNoSlipBcTarget(bcFactory)), testing::Eq(QDevComp27))
        << "The returned boundary condition is not the expected function QDevComp27.";

    bcFactory.setNoSlipBoundaryCondition(BoundaryConditionFactory::NoSlipBC::NoSlip3rdMomentsCompressible);
    EXPECT_THAT(*(getNoSlipBcTarget(bcFactory)), testing::Eq(QDev3rdMomentsComp27))
        << "The returned boundary condition is not the expected function BBDev27.";
}

bcFunction getSlipBcTarget(BoundaryConditionFactory &bcFactory)
{
    auto bc = bcFactory.getSlipBoundaryConditionPost();
    void (*bcTarget)(LBMSimulationParameter *, QforBoundaryConditions *) =
        (*bc.target<void (*)(LBMSimulationParameter *, QforBoundaryConditions *)>());
    return bcTarget;
}

TEST(BoundaryConditionFactoryTest, slipBC)
{
    auto bcFactory = BoundaryConditionFactory();

    bcFactory.setSlipBoundaryCondition(BoundaryConditionFactory::SlipBC::SlipIncompressible);
    EXPECT_THAT(*(getSlipBcTarget(bcFactory)), testing::Eq(QSlipDev27))
        << "The returned boundary condition is not the expected function QSlipDev27.";

    bcFactory.setSlipBoundaryCondition(BoundaryConditionFactory::SlipBC::SlipCompressible);
    EXPECT_THAT(*(getSlipBcTarget(bcFactory)), testing::Eq(QSlipDevComp27))
        << "The returned boundary condition is not the expected function QSlipDevComp27.";

    bcFactory.setSlipBoundaryCondition(BoundaryConditionFactory::SlipBC::SlipCompressibleTurbulentViscosity);
    EXPECT_THAT(*(getSlipBcTarget(bcFactory)), testing::Eq(QSlipDevCompTurbulentViscosity27))
        << "The returned boundary condition is not the expected function QSlipDevCompTurbulentViscosity27.";
}

bcFunction getPressureBcTarget(BoundaryConditionFactory &bcFactory)
{
    auto bc = bcFactory.getPressureBoundaryConditionPre();
    void (*bcTarget)(LBMSimulationParameter *, QforBoundaryConditions *) =
        (*bc.target<void (*)(LBMSimulationParameter *, QforBoundaryConditions *)>());
    return bcTarget;
}

TEST(BoundaryConditionFactoryTest, pressureBC)
{
    auto bcFactory = BoundaryConditionFactory();

    bcFactory.setPressureBoundaryCondition(BoundaryConditionFactory::PressureBC::PressureEquilibrium);
    EXPECT_THAT(*(getPressureBcTarget(bcFactory)), testing::Eq(QPressDev27))
        << "The returned boundary condition is not the expected function QPressDev27.";

    bcFactory.setPressureBoundaryCondition(BoundaryConditionFactory::PressureBC::PressureEquilibrium2);
    EXPECT_THAT(*(getPressureBcTarget(bcFactory)), testing::Eq(QPressDevEQZ27))
        << "The returned boundary condition is not the expected function QPressDevEQZ27.";

    bcFactory.setPressureBoundaryCondition(BoundaryConditionFactory::PressureBC::PressureNonEquilibriumIncompressible);
    EXPECT_THAT(*(getPressureBcTarget(bcFactory)), testing::Eq(QPressDevIncompNEQ27))
        << "The returned boundary condition is not the expected function QPressDevIncompNEQ27.";

    bcFactory.setPressureBoundaryCondition(BoundaryConditionFactory::PressureBC::PressureNonEquilibriumCompressible);
    EXPECT_THAT(*(getPressureBcTarget(bcFactory)), testing::Eq(QPressDevNEQ27))
        << "The returned boundary condition is not the expected function QPressDevNEQ27.";

    bcFactory.setPressureBoundaryCondition(BoundaryConditionFactory::PressureBC::OutflowNonReflective);
    EXPECT_THAT(*(getPressureBcTarget(bcFactory)), testing::Eq(QPressNoRhoDev27))
        << "The returned boundary condition is not the expected function QPressNoRhoDev27.";
}

bcFunction getGeometryBcTarget(BoundaryConditionFactory &bcFactory)
{
    auto bc = bcFactory.getGeometryBoundaryConditionPost();
    void (*bcTarget)(LBMSimulationParameter *, QforBoundaryConditions *) =
        (*bc.target<void (*)(LBMSimulationParameter *, QforBoundaryConditions *)>());
    return bcTarget;
}

TEST(BoundaryConditionFactoryTest, geometryBC)
{
    auto bcFactory = BoundaryConditionFactory();

    // velocity
    bcFactory.setGeometryBoundaryCondition(BoundaryConditionFactory::VelocityBC::VelocityIncompressible);
    EXPECT_THAT(*(getGeometryBcTarget(bcFactory)), testing::Eq(QVelDev27))
        << "The returned boundary condition is not the expected function QVelDev27.";

    bcFactory.setGeometryBoundaryCondition(BoundaryConditionFactory::VelocityBC::VelocityCompressible);
    EXPECT_THAT(*(getGeometryBcTarget(bcFactory)), testing::Eq(QVelDevComp27))
        << "The returned boundary condition is not the expected function QVelDevComp27.";

    bcFactory.setGeometryBoundaryCondition(BoundaryConditionFactory::VelocityBC::VelocityAndPressureCompressible);
    EXPECT_THAT(*(getGeometryBcTarget(bcFactory)), testing::Eq(QVelDevCompZeroPress27))
        << "The returned boundary condition is not the expected function QVelDevCompZeroPress27.";

    // no slip
    bcFactory.setGeometryBoundaryCondition(BoundaryConditionFactory::NoSlipBC::NoSlipImplicitBounceBack);
    auto bc = bcFactory.getGeometryBoundaryConditionPost();
    EXPECT_NO_THROW(bc(nullptr, nullptr)); // empty lambda function should not throw

    bcFactory.setGeometryBoundaryCondition(BoundaryConditionFactory::NoSlipBC::NoSlipIncompressible);
    EXPECT_THAT(*(getGeometryBcTarget(bcFactory)), testing::Eq(QDev27))
        << "The returned boundary condition is not the expected function QDev27.";

    bcFactory.setGeometryBoundaryCondition(BoundaryConditionFactory::NoSlipBC::NoSlipCompressible);
    EXPECT_THAT(*(getGeometryBcTarget(bcFactory)), testing::Eq(QDevComp27))
        << "The returned boundary condition is not the expected function QDevComp27.";

    bcFactory.setGeometryBoundaryCondition(BoundaryConditionFactory::NoSlipBC::NoSlip3rdMomentsCompressible);
    EXPECT_THAT(*(getGeometryBcTarget(bcFactory)), testing::Eq(QDev3rdMomentsComp27))
        << "The returned boundary condition is not the expected function QDev3rdMomentsComp27.";
}

TEST(BoundaryConditionFactoryTest, stressBoundaryConditions)
{
    auto bcFactory = BoundaryConditionFactory();

    bcFactory.setStressBoundaryCondition(BoundaryConditionFactory::StressBC::StressBounceBack);
    auto bc = bcFactory.getStressBoundaryConditionPost();
    auto bcTarget = *bc.target<bcFunctionParamter>();
    EXPECT_THAT(*bcTarget, testing::Eq(BBStressDev27))
        << "The returned boundary condition is not the expected function BBStressDev27.";

    bcFactory.setStressBoundaryCondition(BoundaryConditionFactory::StressBC::StressCompressible);
    bc = bcFactory.getStressBoundaryConditionPost();
    bcTarget = *bc.target<bcFunctionParamter>();
    EXPECT_THAT(*bcTarget, testing::Eq(QStressDevComp27))
        << "The returned boundary condition is not the expected function QStressDevComp27.";
}
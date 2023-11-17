#include <gmock/gmock.h>

#include "ActuatorFarmInlines.h"

TEST(ActuatorFarmInlinesTest, calcNodeIndexInBladeArrays)
{
    const uint numberOfNodesPerBlade = 4;
    const uint numberOfBlades = 3;

    // first node on first blade
    uint bladeNode = 0;
    uint blade = 0;
    uint turbine = 0;
    auto nodeIndexInBladeArrays = calcNodeIndexInBladeArrays(bladeNode, numberOfNodesPerBlade, blade, numberOfBlades, turbine);
    EXPECT_THAT(nodeIndexInBladeArrays, testing::Eq(0));
    nodeIndexInBladeArrays = calcNodeIndexInBladeArrays({turbine, blade, bladeNode}, numberOfNodesPerBlade, numberOfBlades);
    EXPECT_THAT(nodeIndexInBladeArrays, testing::Eq(0));

    // last node on first blade
    bladeNode = 3;
    blade = 0;
    turbine = 0;
    nodeIndexInBladeArrays = calcNodeIndexInBladeArrays(bladeNode, numberOfNodesPerBlade, blade, numberOfBlades, turbine);
    EXPECT_THAT(nodeIndexInBladeArrays, testing::Eq(3));
    nodeIndexInBladeArrays = calcNodeIndexInBladeArrays({turbine, blade, bladeNode}, numberOfNodesPerBlade, numberOfBlades);
    EXPECT_THAT(nodeIndexInBladeArrays, testing::Eq(3));

    // first node on third blade
    bladeNode = 0;
    blade = 2;
    turbine = 0;
    nodeIndexInBladeArrays = calcNodeIndexInBladeArrays(bladeNode, numberOfNodesPerBlade, blade, numberOfBlades, turbine);
    EXPECT_THAT(nodeIndexInBladeArrays, testing::Eq(8));
    nodeIndexInBladeArrays = calcNodeIndexInBladeArrays({turbine, blade, bladeNode}, numberOfNodesPerBlade, numberOfBlades);
    EXPECT_THAT(nodeIndexInBladeArrays, testing::Eq(8));

    // last node on third blade, also last node on first turbine
    bladeNode = 3;
    blade = 2;
    turbine = 0;
    nodeIndexInBladeArrays = calcNodeIndexInBladeArrays(bladeNode, numberOfNodesPerBlade, blade, numberOfBlades, turbine);
    EXPECT_THAT(nodeIndexInBladeArrays, testing::Eq(11));
    nodeIndexInBladeArrays = calcNodeIndexInBladeArrays({turbine, blade, bladeNode}, numberOfNodesPerBlade, numberOfBlades);
    EXPECT_THAT(nodeIndexInBladeArrays, testing::Eq(11));

    // first node on second turbine
    bladeNode = 0;
    blade = 0;
    turbine = 1;
    nodeIndexInBladeArrays = calcNodeIndexInBladeArrays(bladeNode, numberOfNodesPerBlade, blade, numberOfBlades, turbine);
    EXPECT_THAT(nodeIndexInBladeArrays, testing::Eq(12));
    nodeIndexInBladeArrays = calcNodeIndexInBladeArrays({turbine, blade, bladeNode}, numberOfNodesPerBlade, numberOfBlades);
    EXPECT_THAT(nodeIndexInBladeArrays, testing::Eq(12));

    // last node on second turbine
    bladeNode = 3;
    blade = 2;
    turbine = 1;
    nodeIndexInBladeArrays = calcNodeIndexInBladeArrays(bladeNode, numberOfNodesPerBlade, blade, numberOfBlades, turbine);
    EXPECT_THAT(nodeIndexInBladeArrays, testing::Eq(23));
    nodeIndexInBladeArrays = calcNodeIndexInBladeArrays({turbine, blade, bladeNode}, numberOfNodesPerBlade, numberOfBlades);
    EXPECT_THAT(nodeIndexInBladeArrays, testing::Eq(23));
}

TEST(ActuatorFarmInlinesTest, calcTurbineBladeAndBladeNode)
{
    const uint numberOfNodesPerBlade = 4;
    const uint numberOfBlades = 3;

    uint bladeNode;
    uint blade;
    uint turbine;

    TurbineNodeIndex result;

    uint node = 0; // first node on first blade
    calcTurbineBladeAndBladeNode(node, bladeNode, numberOfNodesPerBlade, blade, numberOfBlades, turbine);
    EXPECT_THAT(bladeNode, testing::Eq(0));
    EXPECT_THAT(blade, testing::Eq(0));
    EXPECT_THAT(turbine, testing::Eq(0));
    result = calcTurbineBladeAndBladeNode(node, numberOfNodesPerBlade, numberOfBlades);
    EXPECT_THAT(result.bladeNode, testing::Eq(0));
    EXPECT_THAT(result.blade, testing::Eq(0));
    EXPECT_THAT(result.turbine, testing::Eq(0));

    node = 3; // last node on first blade
    calcTurbineBladeAndBladeNode(node, bladeNode, numberOfNodesPerBlade, blade, numberOfBlades, turbine);
    EXPECT_THAT(bladeNode, testing::Eq(3));
    EXPECT_THAT(blade, testing::Eq(0));
    EXPECT_THAT(turbine, testing::Eq(0));
    result = calcTurbineBladeAndBladeNode(node, numberOfNodesPerBlade, numberOfBlades);
    EXPECT_THAT(result.bladeNode, testing::Eq(3));
    EXPECT_THAT(result.blade, testing::Eq(0));
    EXPECT_THAT(result.turbine, testing::Eq(0));

    node = 8; // first node on third blade
    calcTurbineBladeAndBladeNode(node, bladeNode, numberOfNodesPerBlade, blade, numberOfBlades, turbine);
    EXPECT_THAT(bladeNode, testing::Eq(0));
    EXPECT_THAT(blade, testing::Eq(2));
    EXPECT_THAT(turbine, testing::Eq(0));
    result = calcTurbineBladeAndBladeNode(node, numberOfNodesPerBlade, numberOfBlades);
    EXPECT_THAT(result.bladeNode, testing::Eq(0));
    EXPECT_THAT(result.blade, testing::Eq(2));
    EXPECT_THAT(result.turbine, testing::Eq(0));

    node = 11; // last node on third blade, also last node on first turbine
    calcTurbineBladeAndBladeNode(node, bladeNode, numberOfNodesPerBlade, blade, numberOfBlades, turbine);
    EXPECT_THAT(bladeNode, testing::Eq(3));
    EXPECT_THAT(blade, testing::Eq(2));
    EXPECT_THAT(turbine, testing::Eq(0));
    result = calcTurbineBladeAndBladeNode(node, numberOfNodesPerBlade, numberOfBlades);
    EXPECT_THAT(result.bladeNode, testing::Eq(3));
    EXPECT_THAT(result.blade, testing::Eq(2));
    EXPECT_THAT(result.turbine, testing::Eq(0));

    node = 12; // first node on second turbine
    calcTurbineBladeAndBladeNode(node, bladeNode, numberOfNodesPerBlade, blade, numberOfBlades, turbine);
    EXPECT_THAT(bladeNode, testing::Eq(0));
    EXPECT_THAT(blade, testing::Eq(0));
    EXPECT_THAT(turbine, testing::Eq(1));
    result = calcTurbineBladeAndBladeNode(node, numberOfNodesPerBlade, numberOfBlades);
    EXPECT_THAT(result.bladeNode, testing::Eq(0));
    EXPECT_THAT(result.blade, testing::Eq(0));
    EXPECT_THAT(result.turbine, testing::Eq(1));

    node = 23; // last node on second turbine
    calcTurbineBladeAndBladeNode(node, bladeNode, numberOfNodesPerBlade, blade, numberOfBlades, turbine);
    EXPECT_THAT(bladeNode, testing::Eq(3));
    EXPECT_THAT(blade, testing::Eq(2));
    EXPECT_THAT(turbine, testing::Eq(1));
    result = calcTurbineBladeAndBladeNode(node, numberOfNodesPerBlade, numberOfBlades);
    EXPECT_THAT(result.bladeNode, testing::Eq(3));
    EXPECT_THAT(result.blade, testing::Eq(2));
    EXPECT_THAT(result.turbine, testing::Eq(1));
}

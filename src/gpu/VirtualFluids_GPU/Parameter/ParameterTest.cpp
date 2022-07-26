#include <gmock/gmock.h>

#include <filesystem>
#include <iostream>
#include <string>

#include "Parameter.h"
#include "basics/config/ConfigurationFile.h"

auto RealEq = [](auto value) {
#ifdef VF_DOUBLE_ACCURACY
    return testing::DoubleEq(value);
#else
    return testing::FloatEq(value);
#endif
};

TEST(ParameterTest, passingEmptyFileWithoutPath_ShouldNotThrow)
{
    // assuming that the config files is stored parallel to this file.
    std::filesystem::path filePath = __FILE__;
    filePath.replace_filename("parameterTest_emptyfile.cfg");

    vf::basics::ConfigurationFile config;
    config.load(filePath.string());

    EXPECT_NO_THROW(Parameter para(config, 1, 0));
}

// TODO: test setPossNeighborFilesX
// TODO: test default values

TEST(ParameterTest, check_all_Parameter_CanBePassedToConstructor)
{
    // assuming that the config files is stored parallel to this file.
    std::filesystem::path filePath = __FILE__;
    filePath.replace_filename("parameterTest.cfg");

    vf::basics::ConfigurationFile config;
    config.load(filePath.string());

    Parameter para(config);

    // test optional parameter
    EXPECT_THAT(para.getOutputPath(), testing::Eq("/output/path/"));
    EXPECT_THAT(para.getGridPath(), testing::Eq("/path/to/grid/")); // ... all grid files (e.g. multi-gpu/ multi-level) could be tested as well
    EXPECT_THAT(para.getgeoVec(), testing::Eq("/path/to/grid/geoVec.dat"));
    EXPECT_THAT(para.getMaxDev(), testing::Eq(2));
    EXPECT_THAT(para.getDevices(), testing::ElementsAreArray({ 2, 3 }));
    EXPECT_THAT(para.getOutputPrefix(), testing::Eq("MyPrefix"));
    EXPECT_THAT(para.getPrintFiles(), testing::Eq(true));
    EXPECT_THAT(para.getIsGeometryValues(), testing::Eq(true));
    EXPECT_THAT(para.getCalc2ndOrderMoments(), testing::Eq(true));
    EXPECT_THAT(para.getCalc3rdOrderMoments(), testing::Eq(true));
    EXPECT_THAT(para.getCalcHighOrderMoments(), testing::Eq(true));
    EXPECT_THAT(para.getCalcMedian(), testing::Eq(true));
    EXPECT_THAT(para.getCalcCp(), testing::Eq(true));
    EXPECT_THAT(para.getCalcDragLift(), testing::Eq(true));
    EXPECT_THAT(para.getWriteVeloASCIIfiles(), testing::Eq(true));
    EXPECT_THAT(para.getCalcPlaneConc(), testing::Eq(true));
    EXPECT_THAT(para.getConcFile(), testing::Eq(true));
    EXPECT_THAT(para.isStreetVelocityFile(), testing::Eq(true));
    EXPECT_THAT(para.getUseMeasurePoints(), testing::Eq(true));
    EXPECT_THAT(para.getUseWale(), testing::Eq(true));
    EXPECT_THAT(para.getUseInitNeq(), testing::Eq(true));
    EXPECT_THAT(para.getSimulatePorousMedia(), testing::Eq(true));

    EXPECT_THAT(para.getD3Qxx(), testing::Eq(99));
    EXPECT_THAT(para.getTimestepEnd(), testing::Eq(33));
    EXPECT_THAT(para.getTimestepOut(), testing::Eq(22));
    EXPECT_THAT(para.getTimestepStartOut(), testing::Eq(11));
    EXPECT_THAT(para.getTimeCalcMedStart(), testing::Eq(22));
    EXPECT_THAT(para.getTimeCalcMedEnd(), testing::Eq(44));
    EXPECT_THAT(para.getPressInID(), testing::Eq(25));
    EXPECT_THAT(para.getPressOutID(), testing::Eq(26));
    EXPECT_THAT(para.getPressInZ(), testing::Eq(27));
    EXPECT_THAT(para.getPressOutZ(), testing::Eq(28));

    EXPECT_THAT(para.getDiffOn(), testing::Eq(true));
    EXPECT_THAT(para.getDiffMod(), testing::Eq(99));
    EXPECT_THAT(para.getDiffusivity(), RealEq(1.11));
    EXPECT_THAT(para.getTemperatureInit(), RealEq(2.22));
    EXPECT_THAT(para.getTemperatureBC(), RealEq(3.33));

    EXPECT_THAT(para.getViscosity(), RealEq(4.44));
    EXPECT_THAT(para.getVelocity(), RealEq(5.55));
    EXPECT_THAT(para.getViscosityRatio(), RealEq(6.66));
    EXPECT_THAT(para.getVelocityRatio(), RealEq(7.77));
    EXPECT_THAT(para.getDensityRatio(), RealEq(8.88));
    EXPECT_THAT(para.getPressureRatio(), RealEq(9.99));

    EXPECT_THAT(para.getRealX(), RealEq(0.1));
    EXPECT_THAT(para.getRealY(), RealEq(0.2));
    EXPECT_THAT(para.getFactorPressBC(), RealEq(0.3));

    EXPECT_THAT(para.getReadGeo(), testing::Eq(true));
    EXPECT_THAT(para.getGeometryFileC(), testing::Eq("/pass/to/c"));
    EXPECT_THAT(para.getGeometryFileM(), testing::Eq("/pass/to/m"));
    EXPECT_THAT(para.getGeometryFileF(), testing::Eq("/pass/to/f"));

    EXPECT_THAT(para.getclockCycleForMP(), RealEq(0.4));
    EXPECT_THAT(para.getTimestepForMP(), testing::Eq(4));

    std::vector<real> forces{ 2.0, 2.1, 2.2 };
    double *forces_actual = para.getForcesDouble();
    for (size_t i = 0; i < forces.size(); ++i) {
        EXPECT_THAT((real)forces_actual[i], RealEq(forces[i]));
    }

    std::vector<real> limiters{ 3.0, 3.1, 3.2 };
    double *limiters_actual = para.getQuadricLimitersDouble();
    for (size_t i = 0; i < limiters.size(); ++i) {
        EXPECT_THAT((real)limiters_actual[i], RealEq(limiters[i]));
    }

    EXPECT_THAT(para.getCalcParticles(), testing::Eq(true));
    EXPECT_THAT(para.getParticleBasicLevel(), testing::Eq(1));
    EXPECT_THAT(para.getParticleInitLevel(), testing::Eq(2));
    EXPECT_THAT(para.getNumberOfParticles(), testing::Eq(1111));
    EXPECT_THAT(para.getStartXHotWall(), RealEq(4.1));
    EXPECT_THAT(para.getEndXHotWall(), RealEq(4.2));

    EXPECT_THAT(para.getTimeDoCheckPoint(), testing::Eq(33));
    EXPECT_THAT(para.getTimeDoRestart(), testing::Eq(44));
    EXPECT_THAT(para.getDoCheckPoint(), testing::Eq(true));
    EXPECT_THAT(para.getDoRestart(), testing::Eq(true));
    EXPECT_THAT(para.getMaxLevel(), testing::Eq(1)); // NOGL - 1

    EXPECT_THAT(para.getGridX(), testing::ElementsAreArray({ 100, 101 }));
    EXPECT_THAT(para.getGridY(), testing::ElementsAreArray({ 200, 201 }));
    EXPECT_THAT(para.getGridZ(), testing::ElementsAreArray({ 300, 301 }));
    EXPECT_THAT(para.getDistX(), testing::ElementsAreArray({ 400, 401 }));
    EXPECT_THAT(para.getDistY(), testing::ElementsAreArray({ 500, 501 }));
    EXPECT_THAT(para.getDistZ(), testing::ElementsAreArray({ 600, 601 }));

    EXPECT_THAT(para.getMainKernel(), testing::Eq("KernelName"));
    EXPECT_THAT(para.getMultiKernelOn(), testing::Eq(true));
    EXPECT_THAT(para.getMultiKernelLevel(), testing::ElementsAreArray({ 3, 2, 1 }));

    std::vector<std::string> kernel{ "Kernel1", "Kernel2", "Kernel3" };
    auto kernel_actual = para.getMultiKernel();
    for (size_t i = 0; i < kernel.size(); ++i) {
        EXPECT_THAT(kernel_actual[i], testing::Eq(kernel[i]));
    }

    EXPECT_THAT(para.getCoarse(), testing::Eq(0));
    EXPECT_THAT(para.getFine(), testing::Eq(1)); // NOGL - 1
    EXPECT_THAT(para.parH.size(), testing::Eq(2));
    EXPECT_THAT(para.parD.size(), testing::Eq(2));
}

TEST(ParameterTest, defaultGridPath)
{
    Parameter para;
    EXPECT_THAT(para.getGridPath(), testing::Eq("grid/"));
    EXPECT_THAT(para.getConcentration(), testing::Eq("grid/conc.dat"));
}

TEST(ParameterTest, defaultGridPathMultiGPU)
{
    Parameter para(2, 1);

    EXPECT_THAT(para.getGridPath(), testing::Eq("grid/1/"));
    EXPECT_THAT(para.getConcentration(), testing::Eq("grid/1/conc.dat"));
}

TEST(ParameterTest, setGridPathOverridesDefaultGridPath)
{
    Parameter para(2, 1);
    para.setGridPath("gridPathTest");

    EXPECT_THAT( para.getGridPath(), testing::Eq("gridPathTest/1/"));
    EXPECT_THAT(para.getConcentration(), testing::Eq("gridPathTest/1/conc.dat"));
}

TEST(ParameterTest, setGridPathOverridesConfigFile)
{
    // assuming that the config files is stored parallel to this file.
    std::filesystem::path filePath = __FILE__;
    filePath.replace_filename("parameterTest.cfg");
    vf::basics::ConfigurationFile config;
    config.load(filePath.string());
    auto para = Parameter(config, 2, 0);
    para.setGridPath("gridPathTest");

    EXPECT_THAT( para.getGridPath(), testing::Eq("gridPathTest/0/"));
    EXPECT_THAT(para.getConcentration(), testing::Eq("gridPathTest/0/conc.dat"));

}

TEST(ParameterTest, userMissedSlash)
{
    Parameter para;
    para.setGridPath("gridPathTest");

    EXPECT_THAT(para.getGridPath(), testing::Eq("gridPathTest/"));
    EXPECT_THAT(para.getConcentration(), testing::Eq("gridPathTest/conc.dat"));

}

TEST(ParameterTest, userMissedSlashMultiGPU)
{
    Parameter para(2, 0);
    para.setGridPath("gridPathTest");

    EXPECT_THAT(para.getGridPath(), testing::Eq("gridPathTest/0/"));
    EXPECT_THAT(para.getConcentration(), testing::Eq("gridPathTest/0/conc.dat"));
}
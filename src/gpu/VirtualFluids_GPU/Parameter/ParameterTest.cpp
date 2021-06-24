#include <gmock/gmock.h>

#include <string>

#include "Parameter.h"
#include <basics/config/ConfigurationFile.h>


auto RealEq = [](auto value) { 
#ifdef VF_DOUBLE_ACCURACY
    return testing::DoubleEq(value); 
#else 
    return testing::FloatEq(value);
#endif
};


TEST(ParameterTest, passingEmptyFileWithoutPath_ShouldThrow)
{
    vf::basics::ConfigurationFile config;
    std::string targetPath = __FILE__;
    targetPath = targetPath.substr(0, targetPath.find_last_of('/') + 1);
    targetPath += "parameterTest_emptyfile.cfg";

    config.load(targetPath);

    EXPECT_THROW(Parameter para(config, 1, 0), std::runtime_error);
}

TEST(ParameterTest, check_outputPath)
{
    vf::basics::ConfigurationFile config;
    std::string targetPath = __FILE__;
    targetPath = targetPath.substr(0, targetPath.find_last_of('/') + 1);
    targetPath += "parameterTest.cfg";

    config.load(targetPath);

    Parameter para(config, 1, 0);

    // this two parameters need to be defined in each config file
    EXPECT_THAT(para.getOutputPath(), testing::Eq("/output/path"));
    EXPECT_THAT(para.getgeoVec(), testing::Eq("/path/to/grid/geoVec.dat"));
    // ... all grid files could be tested as well

    // test optional parameter
    EXPECT_THAT(para.getMaxDev(), testing::Eq(2));
    EXPECT_THAT(para.getDevices(), testing::ElementsAreArray({2,3}));
    EXPECT_THAT(para.getOutputPrefix(), testing::Eq("MyPrefix"));
    EXPECT_THAT(para.getPrintFiles(), testing::Eq(true));
    EXPECT_THAT(para.getIsGeometryValues(), testing::Eq(true));
    EXPECT_THAT(para.getCalc2ndOrderMoments(), testing::Eq(true));
    EXPECT_THAT(para.getCalc3rdOrderMoments(), testing::Eq(true));
    EXPECT_THAT(para.getCalcHighOrderMoments(), testing::Eq(true));
    EXPECT_THAT(para.getCalcMedian(), testing::Eq(true));
    EXPECT_THAT(para.getCalcCp(), testing::Eq(true));
    EXPECT_THAT(para.getCalcDragLift(), testing::Eq(true));
}




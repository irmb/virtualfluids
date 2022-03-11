#include <gmock/gmock.h>

#include <iostream>
#include <algorithm>
#include <filesystem>

#include <basics/config/ConfigurationFile.h>
#include <Parameter/Parameter.h>

#include <DataStructureInitializer/GridReaderGenerator/IndexRearrangementForStreams.h>
#include <gpu/GridGenerator/grid/GridBuilder/LevelGridBuilder.h>
#include <gpu/GridGenerator/grid/GridImp.h>


auto RealEq = [](auto value) { 
#ifdef VF_DOUBLE_ACCURACY
    return testing::DoubleEq(value); 
#else 
    return testing::FloatEq(value);
#endif
};


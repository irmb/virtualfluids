//=======================================================================================
// ____          ____    __    ______     __________   __      __       __        __
// \    \       |    |  |  |  |   _   \  |___    ___| |  |    |  |     /  \      |  |
//  \    \      |    |  |  |  |  |_)   |     |  |     |  |    |  |    /    \     |  |
//   \    \     |    |  |  |  |   _   /      |  |     |  |    |  |   /  /\  \    |  |
//    \    \    |    |  |  |  |  | \  \      |  |     |   \__/   |  /  ____  \   |  |____
//     \    \   |    |  |__|  |__|  \__\     |__|      \________/  /__/    \__\  |_______|
//      \    \  |    |   ________________________________________________________________
//       \    \ |    |  |  ______________________________________________________________|
//        \    \|    |  |  |         __          __     __     __     ______      _______
//         \         |  |  |_____   |  |        |  |   |  |   |  |   |   _  \    /  _____)
//          \        |  |   _____|  |  |        |  |   |  |   |  |   |  | \  \   \_______
//           \       |  |  |        |  |_____   |   \_/   |   |  |   |  |_/  /    _____  |
//            \ _____|  |__|        |________|   \_______/    |__|   |______/    (_______/
//
//  This file is part of VirtualFluids. VirtualFluids is free software: you can
//  redistribute it and/or modify it under the terms of the GNU General Public
//  License as published by the Free Software Foundation, either version 3 of
//  the License, or (at your option) any later version.
//
//  VirtualFluids is distributed in the hope that it will be useful, but WITHOUT
//  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//  FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
//  for more details.
//
//  SPDX-License-Identifier: GPL-3.0-or-later
//  SPDX-FileCopyrightText: Copyright © VirtualFluids Project contributors, see AUTHORS.md in root folder
//
//! \addtogroup gpu_kernel_tests Kernel
//! \ingroup gpu_core_tests core
//! \{

#include "Kernel/KernelImp.h"
#include "Kernel/Compressible/NavierStokes/K17/K17CompressibleNavierStokes.h"
#include "Kernel/KernelFactory/KernelFactoryImp.h"
#include "Parameter/Parameter.h"
#include <basics/tests/LogRedirector.h>
#include <basics/tests/testUtilities.h>
#include <logger/Logger.h>

#include "../testUtilitiesGPU.h"

using namespace vf::collisionKernel;

class KernelImpTest : public testing::TestWithParam<std::string>
{
protected:
    void SetUp() override
    {
        kernelName = GetParam();
        SPtr<Parameter> para = testing::vf::createParameterForLevel(level);
        para->getParD(level)->numberofthreads = 2;
        para->getParD(level)->numberOfNodes = 2;
        KernelFactoryImp kernelFactory;
        kernel = kernelFactory.makeKernel(para, kernelName, level);
    }

    static bool containsCriticalTooHigh(std::string& stringForTest)
    {
        bool containsCritical = stringForTest.find("critical") != std::string::npos;
        bool containsTooHigh = stringForTest.find("too high") != std::string::npos;
        return containsCritical && containsTooHigh;
    }

    static bool containsWarningForRecommendedValue(std::string& stringForTest)
    {
        bool containsWarning = stringForTest.find("warning") != std::string::npos;
        bool containsRecommended = stringForTest.find("larger than the recommended value") != std::string::npos;
        return containsWarning && containsRecommended;
    }

    static bool containsInfoOnRecommendedValue(std::string& stringForTest)
    {
        return stringForTest.find("is recommended") != std::string::npos;
    }

    static bool containsInfoOnHardLimit(std::string& stringForTest)
    {
        return stringForTest.find("The hard limit for") != std::string::npos;
    }

    static bool containsNoViscosityMaximumDefined(std::string& stringForTest)
    {
        return stringForTest.find("no viscosity maximum defined") != std::string::npos;
    }

    const uint level = 0;
    const real viscosityLBOnFinestLevelValid = 0.0009;
    const real velocityLBValid = 0.049;
    std::string kernelName;
    SPtr<Kernel> kernel;
    testing::vf::LogRedirector logRedirector;
};

TEST_P(KernelImpTest, VelocityAndViscosityAreOk_expectNoWarning)
{
    kernel->checkKernelParameters(level, velocityLBValid, viscosityLBOnFinestLevelValid);

    std::string log = logRedirector.getLoggerOutput();
    EXPECT_FALSE(containsCriticalTooHigh(log)) << "recorded log:\n" << log;
    EXPECT_FALSE(containsWarningForRecommendedValue(log)) << "recorded log:\n" << log;
    EXPECT_FALSE(containsInfoOnHardLimit(log)) << "recorded log:\n" << log;
    EXPECT_FALSE(containsInfoOnRecommendedValue(log)) << "recorded log:\n" << log;

    if (kernelName == compressible::BGK || kernelName == compressible::BGKPlus || kernelName == incompressible::BGK ||
        kernelName == incompressible::BGKPlus) {
        EXPECT_TRUE(containsNoViscosityMaximumDefined(log)) << "recorded log:\n" << log;
    } else {
        EXPECT_FALSE(containsNoViscosityMaximumDefined(log)) << "recorded log:\n" << log;
    }
}

INSTANTIATE_TEST_SUITE_P(Compressible_VelocityAndViscosityAreOk, KernelImpTest,
                         testing::ValuesIn(compressible::listOfKernels));

INSTANTIATE_TEST_SUITE_P(Incompressible_VelocityAndViscosityAreOk, KernelImpTest,
                         testing::ValuesIn(incompressible::listOfKernels));

TEST_P(KernelImpTest, VelocityIsTooHigh_expectCriticalOutput)
{
    real velocityLBTest = 0.423;

    kernel->checkKernelParameters(level, velocityLBTest, viscosityLBOnFinestLevelValid);

    std::string log = logRedirector.getLoggerOutput();

    EXPECT_TRUE(containsCriticalTooHigh(log)) << "recorded log:\n" << log;
    EXPECT_TRUE(containsInfoOnRecommendedValue(log)) << "recorded log:\n" << log;
    EXPECT_FALSE(containsWarningForRecommendedValue(log)) << "recorded log:\n" << log;
    EXPECT_FALSE(containsInfoOnHardLimit(log)) << "recorded log:\n" << log;
}

INSTANTIATE_TEST_SUITE_P(Compressible_VelocityIsTooHigh_expectCriticalOutput, KernelImpTest,
                         testing::ValuesIn(compressible::listOfKernels));

INSTANTIATE_TEST_SUITE_P(Incompressible_VelocityIsTooHigh_expectCriticalOutput, KernelImpTest,
                         testing::ValuesIn(incompressible::listOfKernels));

TEST_P(KernelImpTest, VelocityIsAboveRecommended_expectWarningWithRecommendation)
{
    real velocityLBTest = 0.10001;

    kernel->checkKernelParameters(level, velocityLBTest, viscosityLBOnFinestLevelValid);

    std::string log = logRedirector.getLoggerOutput();

    EXPECT_TRUE(containsWarningForRecommendedValue(log)) << "recorded log:\n" << log;
    EXPECT_TRUE(containsInfoOnHardLimit(log)) << "recorded log:\n" << log;
    EXPECT_FALSE(containsCriticalTooHigh(log)) << "recorded log:\n" << log;
    EXPECT_FALSE(containsInfoOnRecommendedValue(log)) << "recorded log:\n" << log;
}

INSTANTIATE_TEST_SUITE_P(Compressible_VelocityIsAboveRecommended_expectWarningWithRecommendation, KernelImpTest,
                         testing::ValuesIn(compressible::listOfKernels));

INSTANTIATE_TEST_SUITE_P(Incompressible_VelocityIsAboveRecommended_expectWarningWithRecommendation, KernelImpTest,
                         testing::ValuesIn(incompressible::listOfKernels));

TEST_P(KernelImpTest, ViscosityIsTooHigh_expectCriticalOutput)
{
    real viscosityLBOnFinestLevelTest = 1.1 / 6;

    if (kernelName == compressible::K17CompressibleNavierStokes) {
        viscosityLBOnFinestLevelTest = 1.1 / 42;
    }

    kernel->checkKernelParameters(level, velocityLBValid, viscosityLBOnFinestLevelTest);

    std::string log = logRedirector.getLoggerOutput();

    if (kernelName == compressible::K17CompressibleNavierStokes) {
        EXPECT_TRUE(containsCriticalTooHigh(log)) << "recorded log:\n" << log;
        EXPECT_FALSE(containsWarningForRecommendedValue(log)) << "recorded log:\n" << log;
        EXPECT_TRUE(containsInfoOnRecommendedValue(log)) << "recorded log:\n" << log;
        EXPECT_FALSE(containsNoViscosityMaximumDefined(log)) << "recorded log:\n" << log;
    } else if (kernelName == compressible::K15CompressibleNavierStokes || kernelName == incompressible::CumulantK15) {
        EXPECT_FALSE(containsCriticalTooHigh(log)) << "recorded log:\n" << log;
        EXPECT_TRUE(containsWarningForRecommendedValue(log)) << "recorded log:\n" << log;
        EXPECT_FALSE(containsInfoOnRecommendedValue(log)) << "recorded log:\n" << log;
        EXPECT_FALSE(containsNoViscosityMaximumDefined(log)) << "recorded log:\n" << log;
    } else if (kernelName == compressible::BGK || kernelName == compressible::BGKPlus || kernelName == incompressible::BGK ||
               kernelName == incompressible::BGKPlus) {
        EXPECT_FALSE(containsCriticalTooHigh(log)) << "recorded log:\n" << log;
        EXPECT_FALSE(containsWarningForRecommendedValue(log)) << "recorded log:\n" << log;
        EXPECT_FALSE(containsInfoOnRecommendedValue(log)) << "recorded log:\n" << log;
        EXPECT_TRUE(containsNoViscosityMaximumDefined(log)) << "recorded log:\n" << log;
    }
    EXPECT_FALSE(containsInfoOnHardLimit(log)) << "recorded log:\n" << log;
}

INSTANTIATE_TEST_SUITE_P(Compressible_ViscosityIsTooHigh_expectCriticalOutput, KernelImpTest,
                         testing::ValuesIn(compressible::listOfKernels));

INSTANTIATE_TEST_SUITE_P(Incompressible_ViscosityIsTooHigh_expectCriticalOutput, KernelImpTest,
                         testing::ValuesIn(incompressible::listOfKernels));

TEST_P(KernelImpTest, ViscosityIsAboveRecommended_expectWarningWithRecommendation)
{
    real viscosityLBOnFinestLevelTest = 1.1/6;

    if (kernelName == compressible::K17CompressibleNavierStokes) {
        viscosityLBOnFinestLevelTest = 0.002;
    }

    kernel->checkKernelParameters(level, velocityLBValid, viscosityLBOnFinestLevelTest);

    std::string log = logRedirector.getLoggerOutput();

    if (kernelName == compressible::K17CompressibleNavierStokes) {
        EXPECT_TRUE(containsWarningForRecommendedValue(log)) << "recorded log:\n" << log;
        EXPECT_TRUE(containsInfoOnHardLimit(log)) << "recorded log:\n" << log;
        EXPECT_FALSE(containsNoViscosityMaximumDefined(log)) << "recorded log:\n" << log;
    } else if (kernelName == compressible::K15CompressibleNavierStokes || kernelName == incompressible::CumulantK15) {
        EXPECT_TRUE(containsWarningForRecommendedValue(log)) << "recorded log:\n" << log;
        EXPECT_FALSE(containsInfoOnHardLimit(log)) << "recorded log:\n" << log;
        EXPECT_FALSE(containsNoViscosityMaximumDefined(log)) << "recorded log:\n" << log;
    } else if (kernelName == compressible::BGK || kernelName == compressible::BGKPlus || kernelName == incompressible::BGK ||
               kernelName == incompressible::BGKPlus) {
        EXPECT_TRUE(containsNoViscosityMaximumDefined(log)) << "recorded log:\n" << log;
    }
    EXPECT_FALSE(containsInfoOnRecommendedValue(log)) << "recorded log:\n" << log;
    EXPECT_FALSE(containsCriticalTooHigh(log)) << "recorded log:\n" << log;
}

INSTANTIATE_TEST_SUITE_P(Compressible_ViscosityIsAboveRecommended_expectWarningWithRecommendation, KernelImpTest,
                         testing::ValuesIn(compressible::listOfKernels));

INSTANTIATE_TEST_SUITE_P(Incompressible_ViscosityIsAboveRecommended_expectWarningWithRecommendation, KernelImpTest,
                         testing::ValuesIn(incompressible::listOfKernels));

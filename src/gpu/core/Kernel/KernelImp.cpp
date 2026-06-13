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
//! \addtogroup gpu_Kernel Kernel
//! \ingroup gpu_core core
//! \{
#include "KernelImp.h"

#include <spdlog/fmt/fmt.h>

#include "Calculation/Calculation.h"
#include "logger/Logger.h"

namespace vf::gpu {

void KernelImp::runOnIndices(const unsigned int *indices, unsigned int sizeIndices, CollisionTemplate collisionTemplate, CudaStreamIndex streamIndex)
{
    printf("Method not implemented for this Kernel \n");
}

std::vector<PreProcessorType> KernelImp::getPreProcessorTypes()
{
    return myPreProcessorTypes;
}

bool KernelImp::getKernelUsesFluidNodeIndices() const
{
    return this->kernelUsesFluidNodeIndices;
}

std::string KernelImp::realToString(real kernelParameter)
{
    // use fmt library to control the formatting of the number
    return fmt::format("{:1.4g}", kernelParameter);
}

std::string KernelImp::composeWarningForMaximumKernelParameterHardLimit(const std::string& kernelParameterName,
                                                                        real kernelParameter, real maximumHard)
{
    return "The " + kernelParameterName + " (in LB units) is too high. It was set to " + realToString(kernelParameter) +
           " but should be smaller than " + realToString(maximumHard) + " (hard limit).";
}

std::string KernelImp::composeWarningForMaximumKernelParameterRecommendation(const std::string& kernelParameterName,
                                                                             real kernelParameter, real maximumRecommended)
{
    return "The " + kernelParameterName + " is " + realToString(kernelParameter) +
           ", which is larger than the recommended value of " + realToString(maximumRecommended) + ".";
}

std::string KernelImp::composeRecommendationForMaximumKernelParameter(const std::string& kernelParameterName,
                                                                      real maximumRecommended)
{
    return "A " + kernelParameterName + " smaller than " + realToString(maximumRecommended) + " is recommended.";
}

std::string KernelImp::composeInfoOnHardLimitForMaximumKernelParameter(const std::string& kernelParameterName,
                                                                       real maximumHard)
{
    return "The hard limit for the " + kernelParameterName + " is " + realToString(maximumHard) + ".";
}

void KernelImp::checkViscosity(real viscosityLBOnFinestLevel, uint maxLevel) const
{
    const std::string kernelParameterName = "viscosity";
    std::string message = "At level " + std::to_string(maxLevel) + ": ";

    if (!viscosityMaximumHard.has_value() && !viscosityMaximumRecommended.has_value()) {
        VF_LOG_INFO("There is no viscosity maximum defined for this kernel.");
        return;
    }

    if (viscosityMaximumHard.has_value() && viscosityLBOnFinestLevel > viscosityMaximumHard) {
        message += composeWarningForMaximumKernelParameterHardLimit(kernelParameterName, viscosityLBOnFinestLevel,
                                                                    viscosityMaximumHard.value());
        if (viscosityMaximumRecommended.has_value()) {
            message += "\n" + composeRecommendationForMaximumKernelParameter(kernelParameterName,
                                                                             viscosityMaximumRecommended.value());
        }
        VF_LOG_CRITICAL(message);
        return;
    }

    if (viscosityMaximumRecommended.has_value() && viscosityLBOnFinestLevel > viscosityMaximumRecommended.value()) {
        message += composeWarningForMaximumKernelParameterRecommendation(kernelParameterName, viscosityLBOnFinestLevel,
                                                                         viscosityMaximumRecommended.value());
        if (viscosityMaximumHard.has_value()) {
            message +=
                "\n" + composeInfoOnHardLimitForMaximumKernelParameter(kernelParameterName, viscosityMaximumHard.value());
        }
        VF_LOG_WARNING(message);
    }
}

void KernelImp::checkVelocity(real velocityLB) const
{
    const std::string kernelParameterName = "velocity";
    std::string message;

    if (velocityLB > velocityMaximumHard) {
        VF_LOG_CRITICAL(
            composeWarningForMaximumKernelParameterHardLimit(kernelParameterName, velocityLB, velocityMaximumHard) + "\n" +
            composeRecommendationForMaximumKernelParameter(kernelParameterName, velocityMaximumRecommended));
        return;
    }

    if (velocityLB > velocityMaximumRecommended) {
        VF_LOG_WARNING(composeWarningForMaximumKernelParameterRecommendation(kernelParameterName, velocityLB,
                                                                             velocityMaximumRecommended) +
                       "\n" + composeInfoOnHardLimitForMaximumKernelParameter(kernelParameterName, velocityMaximumHard));
        return;
    }
}

void KernelImp::checkKernelParameters(uint maxLevel, real velocityLB, real viscosityLBOnFinestLevel) const
{
    checkVelocity(velocityLB);
    checkViscosity(viscosityLBOnFinestLevel, maxLevel);
}

KernelImp::KernelImp(std::shared_ptr<Parameter> para, int level) : para(std::move(para)), level(level)
{
}

KernelImp::KernelImp(std::shared_ptr<Parameter> para, int level, real viscosityMaximumRecommended, real viscosityMaximumHard)
    : para(std::move(para)), level(level), viscosityMaximumRecommended(viscosityMaximumRecommended),
      viscosityMaximumHard(viscosityMaximumHard)
{
}

KernelImp::KernelImp(std::shared_ptr<Parameter> para, int level, real viscosityMaximumRecommended)
    : para(std::move(para)), level(level), viscosityMaximumRecommended(viscosityMaximumRecommended)
{
}

}

//! \}

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
#ifndef KERNEL_IMP_H
#define KERNEL_IMP_H

#include <memory>
#include <optional>

#include <constants/NumericConstants.h>

#include "Calculation/Calculation.h"
#include "Kernel.h"
#include "cuda_helper/CudaGrid.h"

class Parameter;
class CudaStreamManager;

class KernelImp : public Kernel
{
public:
    void runOnIndices(const unsigned int* indices, unsigned int sizeIndices, CollisionTemplate collisionTemplate,
                      CudaStreamIndex streamIndex = CudaStreamIndex::Legacy) override;

    void checkKernelParameters(uint maxLevel, real velocityLB, real viscosityLBOnFinestLevel) const override;

    std::vector<PreProcessorType> getPreProcessorTypes() override;

    bool getKernelUsesFluidNodeIndices() const;

protected:
    KernelImp(std::shared_ptr<Parameter> para, int level);
    //! \param viscosityMaximumRecommended in LB units
    KernelImp(std::shared_ptr<Parameter> para, int level, real viscosityMaximumRecommended);
    //! \param viscosityMaximumRecommended and \param viscosityMaximumHard are in LB units
    KernelImp(std::shared_ptr<Parameter> para, int level, real viscosityMaximumRecommended, real viscosityMaximumHard);
    KernelImp() = default;

    // the limit of the velocity in LB units is the same for all kernels
    const real velocityMaximumHard = 1 - std::sqrt(vf::basics::constant::c1o3);
    const real velocityMaximumRecommended = 0.1;

    // the limit of the viscosity in LB units depends on the kernel, and may not be defined for some kernels
    const std::optional<real> viscosityMaximumHard = std::nullopt;
    const std::optional<real> viscosityMaximumRecommended = std::nullopt;

    std::shared_ptr<Parameter> para;
    int level;
    std::vector<PreProcessorType> myPreProcessorTypes;
    vf::cuda::CudaGrid cudaGrid;

    bool kernelUsesFluidNodeIndices = false;

private:
    void checkViscosity(real viscosityLBOnFinestLevel, uint maxLevel) const;
    void checkVelocity(real velocityLB) const;

    static std::string realToString(real kernelParameter);

    static std::string composeWarningForMaximumKernelParameterHardLimit(const std::string& kernelParameterName,
                                                                        real kernelParameter, real maximumHard);
    static std::string composeWarningForMaximumKernelParameterRecommendation(const std::string& kernelParameterName,
                                                                             real kernelParameter, real maximumRecommended);

    static std::string composeRecommendationForMaximumKernelParameter(const std::string& kernelParameterName,
                                                                      real maximumRecommended);
    static std::string composeInfoOnHardLimitForMaximumKernelParameter(const std::string& kernelParameterName,
                                                                       real maximumHard);
};

#endif

//! \}

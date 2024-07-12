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
//! \addtogroup gpu_Samplers Sampplers
//! \ingroup gpu_core core
//! \{
//! \author Henry Korb
//! \brief Base class for all samplers
//=======================================================================================

#ifndef SAMPLER_H
#define SAMPLER_H

#include <functional>

#include <basics/DataTypes.h>
#include <basics/PointerDefinitions.h>

class GridProvider;

inline std::string fixOutputPath(const std::string path)
{
    if (path.back() == '/')
        return path;
    return path + "/";
}


//! \brief Base class for all samplers
class Sampler
{
public:
    Sampler(const std::string outputPath,
            const std::string probeName)
        : outputPath(fixOutputPath(outputPath)), probeName(probeName)
    {
    }
    virtual ~Sampler() = default;

    virtual void init() = 0;
    virtual void sample(int level, uint t) = 0;
    virtual void getTaggedFluidNodes(GridProvider* gridProvider) = 0;

protected:
    std::string outputPath;
    std::string probeName;
};

#endif

//! \}
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
//! \addtogroup gpu_PreProcessor PreProcessor
//! \ingroup gpu_core core
//! \{
//! \author Martin Schoenherr
//=======================================================================================
#include "PreProcessorFactoryImp.h"

#include "PreProcessor/PreProcessorImp.h"

#include "PreProcessor/PreProcessorStrategy/InitAdvectionDiffusionCompressible/InitAdvectionDiffusionCompressible.h"
#include "PreProcessor/PreProcessorStrategy/InitNavierStokesCompressible/InitNavierStokesCompressible.h"
#include "PreProcessor/PreProcessorStrategy/InitAdvectionDiffusionIncompressible/InitAdvectionDiffusionIncompressible.h"
#include "PreProcessor/PreProcessorStrategy/InitNavierStokesIncompressible/InitNavierStokesIncompressible.h"


std::shared_ptr<PreProcessor> PreProcessorFactoryImp::makePreProcessor(std::vector<PreProcessorType> preProcessorTypes, std::shared_ptr<Parameter> para)
{
    std::shared_ptr<PreProcessorImp> prePro = PreProcessorImp::getNewInstance();

    for (std::size_t i = 0; i < preProcessorTypes.size(); i++)
        prePro->addStrategy(makePreProcessorStrategy(preProcessorTypes.at(i), para));

    return prePro;
}

std::shared_ptr<PreProcessorStrategy> PreProcessorFactoryImp::makePreProcessorStrategy(PreProcessorType preProcessorType, std::shared_ptr<Parameter> para)
{
    switch (preProcessorType)
    {
    case InitNavierStokesIncompressible:
        return InitNavierStokesIncompressible::getNewInstance(para);
        break;
    case InitNavierStokesCompressible:
        return InitNavierStokesCompressible::getNewInstance(para);
        break;
    case InitAdvectionDiffusionIncompressible:
        return InitAdvectionDiffusionIncompressible::getNewInstance(para);
        break;
    case InitAdvectionDiffusionCompressible:
        return InitAdvectionDiffusionCompressible::getNewInstance(para);
        break;
    default:
        break;
    }
    throw  std::runtime_error("PreProcessorFactory does not know the PreProcessorType.");
}

//! \}

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
//            \ _____|  |__|        |________|   \_______/    |__|   |______/    (_______/)
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
//! \author Hussein Alihussein
//! \brief Shared default parameters for fast-winding solid classification and Q-related helpers.
//! \note Generated with the assistance of OpenAI Codex (GPT-5). Reviewed and adapted by author.
//=======================================================================================
#pragma once

namespace vf::grid_winding
{
//! \brief Default fast-winding accuracy scale used in `computeSolidAngle(...)`.
//!
//! Defaults in this header follow the libigl/HDK fast-winding API used in this project.
//! \ref <a href="https://doi.org/10.1145/3197517.3201337"><b>[ G. Barill et al. (2018), DOI:10.1145/3197517.3201337 ]</b></a>
//! Higher values allow more approximation in the tree traversal, which is usually faster
//! but can smooth small geometry details. Lower values are usually slower but more exact.
inline constexpr float FastWindingDefaultAccuracyScale = 2.0f;

//! \brief Default normalized threshold for solid classification.
//!
//! The winding value is normalized by `1 / (4*pi)`, so a closed, well-oriented surface is
//! usually near `1` inside and near `0` outside. `0.5` is a balanced inside/outside split.
//! Lower values mark more nodes as solid. Higher values mark fewer nodes as solid.
inline constexpr float FastWindingDefaultThreshold = 0.5f;

//! \brief Default comparison tolerance near the threshold.
//!
//! This helps avoid unstable classification for points very close to the surface because of
//! floating-point noise. Larger values are more tolerant (can mark a few more near-surface
//! nodes as solid). Smaller values are stricter but can be noisier.
inline constexpr float FastWindingDefaultTolerance = 1.0e-4f;

//! \brief Default multipole expansion order for libigl/HDK `UT_SolidAngle::init(...)`.
//!
//! This mirrors the library default (`order = 2`). Lower order is usually faster but less
//! accurate. Higher order can improve accuracy but increases setup cost.
inline constexpr int FastWindingDefaultExpansionOrder = 2;

} // namespace vf::grid_winding

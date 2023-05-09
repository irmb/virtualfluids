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
//  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
//  for more details.
//
//  You should have received a copy of the GNU General Public License along
//  with VirtualFluids (see COPYING.txt). If not, see <http://www.gnu.org/licenses/>.
//
//! \file VirtualFluids.cpp
//! \ingroup src
//! \author Henry Korb
//=======================================================================================
#include <pybind11/pybind11.h>
#include "basics/basics.cpp"
#include "lbm/lbm.cpp"
#include "logger/logger.cpp"

#ifdef VF_GPU_PYTHONBINDINGS
#include "gpu/gpu.cpp"
#endif
#ifdef VF_CPU_PYTHONBINDINGS
#include "cpu/cpu.cpp"
#endif


namespace py_bindings
{
    namespace py = pybind11;

    PYBIND11_MODULE(bindings, m)
    {
        // because we do not use the old logger (src/basics/logger) anymore and cout is not passed anymore to the old logger, we probably do not need this anymore
        // pybind11::add_ostream_redirect(m, "ostream_redirect");
        basics::makeModule(m);
        lbm::makeModule(m);
        logging::makeModule(m);
#ifdef VF_GPU_PYTHONBINDINGS
        gpu::makeModule(m);
#endif
#ifdef VF_CPU_PYTHONBINDINGS
        cpu::makeModule(m);
#endif
    }
}

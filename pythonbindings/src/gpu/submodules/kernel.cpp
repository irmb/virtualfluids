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
//! \file kernel.cpp
//! \ingroup submodules
//! \author Henry Korb
//=======================================================================================
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <gpu/core/Kernel/Utilities/KernelTypes.h>

namespace kernel
{
    namespace py = pybind11;

    void makeModule(py::module_ &parentModule)
    {
        auto kernel_module = parentModule.def_submodule("Kernel", "Kernel types");
        auto compressible = kernel_module.def_submodule("compressible", "Compressible Kernel types");
        auto incompressible = kernel_module.def_submodule("incompressible", "Incompressible Kernel types");
        auto porous_media = kernel_module.def_submodule("porous_media", "Porous Media Kernel types");
        auto wale = kernel_module.def_submodule("wale", "WALE Kernel types");

        compressible.attr("BGK") = vf::collisionKernel::compressible::BGK;
        compressible.attr("BGKUnified") = vf::collisionKernel::compressible::BGKUnified;
        compressible.attr("BGKPlus") = vf::collisionKernel::compressible::BGKPlus;
        compressible.attr("MRT") = vf::collisionKernel::compressible::MRT;
        compressible.attr("Cascade") = vf::collisionKernel::compressible::Cascade;
        compressible.attr("CumulantClassic") = vf::collisionKernel::compressible::CumulantClassic;
        compressible.attr("CumulantK15Unified") = vf::collisionKernel::compressible::CumulantK15Unified;
        compressible.attr("K17CompressibleNavierStokesUnified") = vf::collisionKernel::compressible::K17CompressibleNavierStokesUnified;
        compressible.attr("K17CompressibleNavierStokes") = vf::collisionKernel::compressible::K17CompressibleNavierStokes;
        compressible.attr("K17CompressibleNavierStokesBulkViscosity") = vf::collisionKernel::compressible::K17CompressibleNavierStokesBulkViscosity;
        compressible.attr("K17CompressibleNavierStokesChimeraLegacy") = vf::collisionKernel::compressible::K17CompressibleNavierStokesChimeraLegacy;
        compressible.attr("CumulantAll4SP27") = vf::collisionKernel::compressible::CumulantAll4SP27;
        compressible.attr("CumulantK18") = vf::collisionKernel::compressible::CumulantK18;
        compressible.attr("CumulantK20") = vf::collisionKernel::compressible::CumulantK20;
        compressible.attr("K15CompressibleNavierStokes") = vf::collisionKernel::compressible::K15CompressibleNavierStokes;
        compressible.attr("K15CompressibleNavierStokesBulk") = vf::collisionKernel::compressible::K15CompressibleNavierStokesBulk;
        compressible.attr("K15CompressibleNavierStokesSponge") = vf::collisionKernel::compressible::K15CompressibleNavierStokesSponge;

        incompressible.attr("BGK") = vf::collisionKernel::incompressible::BGK;
        incompressible.attr("BGKPlus") = vf::collisionKernel::incompressible::BGKPlus;
        incompressible.attr("MRT") = vf::collisionKernel::incompressible::MRT;
        incompressible.attr("Cascade") = vf::collisionKernel::incompressible::Cascade;
        incompressible.attr("Cumulant1h") = vf::collisionKernel::incompressible::Cumulant1h;
        incompressible.attr("CumulantIsometric") = vf::collisionKernel::incompressible::CumulantIsometric;
        incompressible.attr("CumulantK15") = vf::collisionKernel::incompressible::CumulantK15;

        porous_media.attr("CumulantOne") = vf::collisionKernel::porousMedia::CumulantOne;
        
        wale.attr("CumulantK17") = vf::collisionKernel::wale::CumulantK17;
        wale.attr("CumulantK17Debug") = vf::collisionKernel::wale::CumulantK17Debug;
        wale.attr("CumulantK15") = vf::collisionKernel::wale::CumulantK15;
        wale.attr("CumulantK15SoniMalav") = vf::collisionKernel::wale::CumulantK15SoniMalav;

    }
}

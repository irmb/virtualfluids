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
//! \author Sven Marcus, Henry Korb
//=======================================================================================
#include "BCStrategy.h"
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <BoundaryConditions/PressureNonEquilibrium.h>
#include <BoundaryConditions/OutflowNonReflecting.h>
#include <BoundaryConditions/BC.h>
#include <BoundaryConditions/NoSlipBC.h>
#include <BoundaryConditions/VelocityBC.h>
#include <BoundaryConditions/PressureBC.h>
#include <BoundaryConditions/NoSlipInterpolated.h>
#include <BoundaryConditions/VelocityInterpolated.h>

namespace boundaryconditions
{
    namespace py = pybind11;
    using namespace py::literals;

    template<class adapter, class algorithm,
            class = std::enable_if_t<std::is_base_of<BC, adapter>::value>,
            class = std::enable_if_t<std::is_base_of<BCStrategy, algorithm>::value>>
    class PyBoundaryCondition : public adapter
    {
    public:
        template<typename ...Args>
        PyBoundaryCondition(Args &&... args) : adapter(std::forward<Args>(args)...)
        {
            this->setBCStrategy(std::make_shared<algorithm>());
        }
    };

    template<class adapter, class algorithm>
    using bc_class = py::class_<PyBoundaryCondition<adapter, algorithm>, BC,
            std::shared_ptr<PyBoundaryCondition<adapter, algorithm>>>;

    void makeModule(py::module_ &parentModule)
    {
        py::module_ bcModule = parentModule.def_submodule("boundaryconditions");

        auto _ = py::class_<BC, std::shared_ptr<BC>>(bcModule, "BC");

        bc_class<NoSlipBC, NoSlipInterpolated>(bcModule, "NoSlipBoundaryCondition")
                .def(py::init());

        bc_class<VelocityBC, VelocityInterpolated>(bcModule, "VelocityBoundaryCondition")
                .def(py::init())
                .def(py::init<bool &, bool &, bool &, mu::Parser &, real &, real &>(),
                     "vx1"_a, "vx2"_a, "vx3"_a,
                     "function"_a, "start_time"_a, "end_time"_a)
                .def(py::init<bool &, bool &, bool &, mu::Parser &, mu::Parser &, mu::Parser &, real &, real &>(),
                     "vx1"_a, "vx2"_a, "vx3"_a,
                     "function_vx1"_a, "function_vx2"_a, "function_vx2"_a,
                     "start_time"_a, "end_time"_a)
                .def(py::init<real &, real &, real &, real &, real &, real &, real &, real &, real &>(),
                     "vx1"_a, "vx1_start_time"_a, "vx1_end_time"_a,
                     "vx2"_a, "vx2_start_time"_a, "vx2_end_time"_a,
                     "vx3"_a, "vx3_start_time"_a, "vx3_end_time"_a);

        bc_class<PressureBC, OutflowNonReflecting>(bcModule, "NonReflectingOutflow")
                .def(py::init());
    }

}
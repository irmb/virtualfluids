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
//! \author Henry Korb
//=======================================================================================
#include <pybind11/cast.h>
#include <pybind11/pybind11.h>

#include <parallel/Communicator.h>
#include <parallel/MPICommunicator.h>

namespace parallel
{
namespace py = pybind11;

PYBIND11_MODULE(parallel, m)
{
py::class_<vf::parallel::Communicator, std::shared_ptr<vf::parallel::Communicator>>(m, "Communicator")
        .def_static("get_instance", &vf::parallel::Communicator::getInstance)
        .def("get_process_id", py::overload_cast<>(&vf::parallel::Communicator::getProcessID, py::const_))
        .def("get_number_of_processes", &vf::parallel::Communicator::getNumberOfProcesses);

    py::class_<vf::parallel::MPICommunicator, vf::parallel::Communicator, std::shared_ptr<vf::parallel::MPICommunicator>>(m, "MPICommunicator")
        .def_static("get_instance", &vf::parallel::MPICommunicator::getInstance);
}
} // namespace parallel

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
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <gpu/core/Samplers/Probe.h>
#include <gpu/core/Samplers/WallModelProbe.h>
#include <gpu/core/Samplers/PlanarAverageProbe.h>
#include <gpu/core/Samplers/Sampler.h>

namespace probes
{
    namespace py = pybind11;

    py::module makeModule(py::module_ &parentModule)
    {
        py::module probeModule = parentModule.def_submodule("probes");

        py::module probeProbeModule = probeModule.def_submodule("Probe");

        py::enum_<Probe::Statistic>(probeProbeModule, "Statistic")
        .value("Instantaneous", Probe::Statistic::Instantaneous)
        .value("Means", Probe::Statistic::Means)
        .value("Variances", Probe::Statistic::Variances);

        py::class_<Probe, Sampler, std::shared_ptr<Probe>>(probeProbeModule, "Probe")
        .def(py::init<  SPtr<Parameter>,
                        SPtr<CudaMemoryManager>,
                        const std::string,
                        const std::string,
                        uint,
                        uint, 
                        uint,
                        uint,
                        bool,
                        bool>(), 
                        py::arg("para"),
                        py::arg("cuda_memory_manager"),
                        py::arg("probe_name"),
                        py::arg("output_path"),
                        py::arg("t_start_avg"),
                        py::arg("t_avg"),
                        py::arg("t_start_out"),
                        py::arg("t_out"),
                        py::arg("output_timeseries"),
                        py::arg("average_every_timestep"))
        .def("add_statistic", &Probe::addStatistic, py::arg("variable"))
        .def("set_file_name_to_n_out", &Probe::setFileNameToNOut)
        .def("add_all_available_statistics", &Probe::addAllAvailableStatistics)
        .def("add_probe_point", &Probe::addProbePoint, py::arg("point_coord_x"), py::arg("point_coord_y"), py::arg("point_coord_z"))
        .def("add_probe_points_from_list", &Probe::addProbePointsFromList, py::arg("point_coords_x"), py::arg("point_coords_y"), py::arg("point_coords_z"))
        .def("set_probe_plane", &Probe::addProbePlane, py::arg("pos_x"), py::arg("pos_y"), py::arg("pos_z"), py::arg("delta_x"), py::arg("delta_y"), py::arg("delta_z"));

        py::module planarAverageProbeModule = probeModule.def_submodule("PlanarAverageProbe");

        py::enum_<PlanarAverageProbe::PlaneNormal>(planarAverageProbeModule, "PlaneNormal")
        .value("x", PlanarAverageProbe::PlaneNormal::x)
        .value("y", PlanarAverageProbe::PlaneNormal::y)
        .value("z", PlanarAverageProbe::PlaneNormal::z);

        py::enum_<PlanarAverageProbe::Statistic>(planarAverageProbeModule, "Statistic")
        .value("Means", PlanarAverageProbe::Statistic::Means)
        .value("Covariances", PlanarAverageProbe::Statistic::Covariances)
        .value("Skewness", PlanarAverageProbe::Statistic::Skewness)
        .value("Flatness", PlanarAverageProbe::Statistic::Flatness);

        py::class_<PlanarAverageProbe, Sampler, std::shared_ptr<PlanarAverageProbe>>(planarAverageProbeModule, "PlanarAverageProbe")
        .def(py::init<  SPtr<Parameter>,
                        SPtr<CudaMemoryManager>,
                        const std::string,
                        const std::string,
                        uint,
                        uint,
                        uint,
                        uint,
                        uint,
                        PlanarAverageProbe::PlaneNormal,
                        bool>(),
                        py::arg("para"),
                        py::arg("cuda_memory_manager"),
                        py::arg("probe_name"),
                        py::arg("output_path"),
                        py::arg("t_start_avg"),
                        py::arg("t_start_tmp_avg"),
                        py::arg("t_avg"),
                        py::arg("t_start_out"),
                        py::arg("t_out"),
                        py::arg("plane_normal"),
                        py::arg("compute_time_averages"));


        py::class_<WallModelProbe, Sampler, std::shared_ptr<WallModelProbe>>(probeModule, "WallModelProbe")
        .def(py::init<  SPtr<Parameter>,
                        SPtr<CudaMemoryManager>,
                        const std::string,
                        const std::string,
                        uint,
                        uint, 
                        uint,
                        uint,
                        uint,
                        bool,
                        bool,
                        bool,
                        bool>(),
                        py::arg("para"),
                        py::arg("cuda_memory_manager"),
                        py::arg("probe_name"),
                        py::arg("output_path"),
                        py::arg("t_start_avg"),
                        py::arg("t_start_tmp_avg"),
                        py::arg("t_avg"),
                        py::arg("t_start_out"),
                        py::arg("t_out"),
                        py::arg("average_every_timestep"),
                        py::arg("compute_temporal_averages"),
                        py::arg("output_stress"),
                        py::arg("evaluate_pressure_gradient"));

        return probeModule;
    }
}
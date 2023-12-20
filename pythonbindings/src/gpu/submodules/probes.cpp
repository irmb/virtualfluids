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
#include <gpu/core/PreCollisionInteractor/Probes/Probe.h>
#include <gpu/core/PreCollisionInteractor/Probes/PointProbe.h>
#include <gpu/core/PreCollisionInteractor/Probes/PlaneProbe.h>
#include <gpu/core/PreCollisionInteractor/Probes/WallModelProbe.h>
#include <gpu/core/PreCollisionInteractor/Probes/PlanarAverageProbe.h>
#include <gpu/core/PreCollisionInteractor/PreCollisionInteractor.h>

namespace probes
{
    namespace py = pybind11;

    py::module makeModule(py::module_ &parentModule)
    {
        py::module probeModule = parentModule.def_submodule("probes");

        py::enum_<Statistic>(probeModule, "Statistic")
        .value("Instantaneous", Statistic::Instantaneous)
        .value("Means", Statistic::Means)
        .value("Variances", Statistic::Variances)
        .value("SpatialMeans", Statistic::SpatialMeans)
        .value("SpatioTemporalMeans", Statistic::SpatioTemporalMeans)
        .value("SpatialCovariances", Statistic::SpatialCovariances)
        .value("SpatioTemporalCovariances", Statistic::SpatioTemporalCovariances)
        .value("SpatialSkewness", Statistic::SpatialSkewness)
        .value("SpatioTemporalSkewness", Statistic::SpatioTemporalSkewness)
        .value("SpatialFlatness", Statistic::SpatialFlatness)
        .value("SpatioTemporalFlatness", Statistic::SpatioTemporalFlatness);

        py::class_<Probe, PreCollisionInteractor, std::shared_ptr<Probe>>(probeModule, "Probe")
        .def("add_statistic", &Probe::addStatistic, py::arg("variable"))
        .def("set_file_name_to_n_out", &Probe::setFileNameToNOut)
        .def("add_all_available_statistics", &Probe::addAllAvailableStatistics);

        py::class_<PointProbe, Probe, std::shared_ptr<PointProbe>>(probeModule, "PointProbe")
        .def(py::init<
                        const std::string,
                        const std::string,
                        uint,
                        uint, 
                        uint,
                        uint,
                        bool>(), 
                        py::arg("probe_name"),
                        py::arg("output_path"),
                        py::arg("t_start_avg"),
                        py::arg("t_avg"),
                        py::arg("t_start_out"),
                        py::arg("t_out"),
                        py::arg("output_timeseries"))
        .def("add_probe_point", &PointProbe::addProbePoint, py::arg("point_coord_x"), py::arg("point_coord_y"), py::arg("point_coord_z"))
        .def("add_probe_points_from_list", &PointProbe::addProbePointsFromList, py::arg("point_coords_x"), py::arg("point_coords_y"), py::arg("point_coords_z"));

        py::class_<PlaneProbe, Probe, std::shared_ptr<PlaneProbe>>(probeModule, "PlaneProbe")
        .def(py::init<
                        const std::string,
                        const std::string,
                        uint,
                        uint, 
                        uint,
                        uint>(), 
                        py::arg("probe_name"),
                        py::arg("output_path"),
                        py::arg("t_start_avg"),
                        py::arg("t_avg"),
                        py::arg("t_start_out"),
                        py::arg("t_out"))
        .def("set_probe_plane", &PlaneProbe::setProbePlane, py::arg("pos_x"), py::arg("pos_y"), py::arg("pos_z"), py::arg("delta_x"), py::arg("delta_y"), py::arg("delta_z"));

        py::class_<PlanarAverageProbe, Probe, std::shared_ptr<PlanarAverageProbe>>(probeModule, "PlanarAverageProbe")
        .def(py::init<
                        const std::string,
                        const std::string,
                        uint,
                        uint,
                        uint,
                        uint,
                        uint,
                        char>(),
                        py::arg("probe_name"),
                        py::arg("output_path"),
                        py::arg("t_start_avg"),
                        py::arg("t_start_tmp_avg"),
                        py::arg("t_avg"),
                        py::arg("t_start_out"),
                        py::arg("t_out"),
                        py::arg("plane_normal"));


        py::class_<WallModelProbe, Probe, std::shared_ptr<WallModelProbe>>(probeModule, "WallModelProbe")
        .def(py::init<
                        const std::string,
                        const std::string,
                        uint,
                        uint, 
                        uint,
                        uint,
                        uint>(), 
                        py::arg("probe_name"),
                        py::arg("output_path"),
                        py::arg("t_start_avg"),
                        py::arg("t_start_tmp_avg"),
                        py::arg("t_avg"),
                        py::arg("t_start_out"),
                        py::arg("t_out"))
        .def("set_force_output_to_stress", &WallModelProbe::setForceOutputToStress, py::arg("output_stress"))
        .def("set_evaluate_pressure_gradient", &WallModelProbe::setEvaluatePressureGradient, py::arg("eval_press_grad"));

        return probeModule;
    }
}
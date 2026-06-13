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
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <gpu/core/PreCollisionInteractor/Actuator/ActuatorFarm.h>
#include <gpu/core/PreCollisionInteractor/Actuator/ActuatorFarmStandalone.h>
#include <gpu/core/PreCollisionInteractor/Actuator/ActuatorFarmStandaloneVAWT.h>
#include <gpu/core/PreCollisionInteractor/PreCollisionInteractor.h>


class PyActuatorFarm : public vf::gpu::ActuatorFarm 
{
public:
    using ActuatorFarm::ActuatorFarm; // Inherit constructors
    void updateForcesAndCoordinates(real time, real deltaT) override 
    { 
        PYBIND11_OVERRIDE_PURE_NAME(void, ActuatorFarm, "update_forces_and_coordinates", updateForcesAndCoordinates, time, deltaT); 
    }
};

namespace gpu_bindings::actuator_farm
{
    namespace py = pybind11;
    using namespace vf::gpu;

    template<class dtype>
    dtype* np_to_arr(py::array_t<dtype, py::array::c_style> array){ return static_cast<dtype *>(array.request().ptr); };

    template<class dtype>
    intptr_t arr_to_cp(dtype* array){ return reinterpret_cast<intptr_t>(array); };

    inline void makeModule(py::module_ &parentModule)
    {
        using arr = py::array_t<real, py::array::c_style>;

        using HubConfig = ActuatorFarm::HubConfig;
        using TowerConfig = ActuatorFarm::TowerConfig;
        using VAWTConfig = ActuatorFarm::VAWTConfig;

        py::class_<HubConfig>(parentModule, "HubConfig")
            .def(py::init<real, real, real, uint>(),
                py::arg("length"),
                py::arg("radius"),
                py::arg("position_offset"),
                py::arg("number_of_points_per_turbine"))
            .def_readwrite("length", &HubConfig::length)
            .def_readwrite("radius", &HubConfig::radius)
            .def_readwrite("position_offset", &HubConfig::positionOffset)
            .def_readwrite("number_of_points_per_turbine", &HubConfig::numberOfPointsPerTurbine);

        py::class_<TowerConfig>(parentModule, "TowerConfig")
            .def(py::init<real, real, std::vector<real>, uint>(),
                py::arg("radius"),
                py::arg("offset"),
                py::arg("tower_heights"),
                py::arg("number_of_points_per_turbine"))
            .def_readwrite("radius", &TowerConfig::radius)
            .def_readwrite("offset", &TowerConfig::offset)
            .def_readwrite("tower_heights", &TowerConfig::towerHeights)
            .def_readwrite("number_of_points_per_turbine", &TowerConfig::numberOfPointsPerTurbine);

        py::class_<VAWTConfig>(parentModule, "VAWTConfig")
            .def(py::init([](real rotorHeight, bool flagLocalSmearingWidth) {
                return VAWTConfig{rotorHeight, flagLocalSmearingWidth};
            }),
                py::arg("rotor_height"),
                py::arg("flag_local_smearing_width"))
            .def_readwrite("rotor_height", &VAWTConfig::rotorHeight)
            .def_readwrite("flag_local_smearing_width", &VAWTConfig::flagLocalSmearingWidth);

        py::class_<ActuatorFarm, PreCollisionInteractor, PyActuatorFarm, std::shared_ptr<ActuatorFarm>>(parentModule, "ActuatorFarm", py::dynamic_attr())
        .def(py::init<  SPtr<Parameter>,
                        SPtr<CudaMemoryManager>,
                        const real,
                        const std::vector<real>&,
                        const std::vector<real>&,
                        const std::vector<real>&,
                        const std::vector<real>&,
                        const real,
                        const int,
                        const bool,
                        const uint,
                        const std::optional<ActuatorFarm::HubConfig>&,
                        const std::optional<ActuatorFarm::TowerConfig>&,
                        const std::optional<ActuatorFarm::VAWTConfig>&>(),
                        py::arg("para"),
                        py::arg("cuda_memory_manager"),
                        py::arg("diameter"),
                        py::arg("blade_radii"),
                        py::arg("turbine_positions_x"),
                        py::arg("turbine_positions_y"),
                        py::arg("turbine_positions_z"),
                        py::arg("smearing_width"),
                        py::arg("level"),
                        py::arg("use_host_arrays"),
                        py::arg("number_of_blades") = 3,
                        py::arg("hub_config") = py::none(),
                        py::arg("tower_config") = py::none(),
                        py::arg("vawt_config") = py::none())
        .def_property_readonly("number_of_points_per_blade", &ActuatorFarm::getNumberOfPointsPerBlade)
        .def_property_readonly("number_of_turbines", &ActuatorFarm::getNumberOfTurbines)
        .def_property_readonly("number_of_blade_points_per_turbine", &ActuatorFarm::getNumberOfBladePointsPerTurbine)
        .def_property_readonly("number_of_blades_per_turbine", &ActuatorFarm::getNumberOfBladesPerTurbine)
        .def_property_readonly("number_of_blade_points", &ActuatorFarm::getNumberOfBladePoints)
        .def_property_readonly("number_of_indices", &ActuatorFarm::getNumberOfIndices)
        .def_property_readonly("total_number_of_points", &ActuatorFarm::getTotalNumberOfPoints)
        .def_property_readonly("number_of_hub_points", &ActuatorFarm::getNumberOfHubPoints)
        .def_property_readonly("number_of_tower_points", &ActuatorFarm::getNumberOfTowerPoints)

        .def("get_turbine_pos", [](ActuatorFarm& al, uint turbine){ 
            real position[3] = {al.getTurbinePosX(turbine), al.getTurbinePosY(turbine), al.getTurbinePosZ(turbine)}; return arr(3,  position);
            }, py::arg("turbine"))
        .def("get_all_turbine_pos_x", [](ActuatorFarm& al){ return arr(al.getNumberOfTurbines(), al.getAllTurbinePosX()); } )
        .def("get_all_turbine_pos_y", [](ActuatorFarm& al){ return arr(al.getNumberOfTurbines(), al.getAllTurbinePosY()); } )
        .def("get_all_turbine_pos_z", [](ActuatorFarm& al){ return arr(al.getNumberOfTurbines(), al.getAllTurbinePosZ()); } )
    
        .def("get_all_blade_coords_x", [](ActuatorFarm& al){ return arr({al.getNumberOfTurbines(), al.getNumberOfBladesPerTurbine(), al.getNumberOfPointsPerBlade()}, al.getAllBladeCoordsX()); } )
        .def("get_all_blade_coords_y", [](ActuatorFarm& al){ return arr({al.getNumberOfTurbines(), al.getNumberOfBladesPerTurbine(), al.getNumberOfPointsPerBlade()}, al.getAllBladeCoordsY()); } )
        .def("get_all_blade_coords_z", [](ActuatorFarm& al){ return arr({al.getNumberOfTurbines(), al.getNumberOfBladesPerTurbine(), al.getNumberOfPointsPerBlade()}, al.getAllBladeCoordsZ()); } )        
        .def("get_all_blade_velocities_x", [](ActuatorFarm& al){ return arr({al.getNumberOfTurbines(), al.getNumberOfBladesPerTurbine(), al.getNumberOfPointsPerBlade()}, al.getAllBladeVelocitiesX()); } )
        .def("get_all_blade_velocities_y", [](ActuatorFarm& al){ return arr({al.getNumberOfTurbines(), al.getNumberOfBladesPerTurbine(), al.getNumberOfPointsPerBlade()}, al.getAllBladeVelocitiesY()); } )
        .def("get_all_blade_velocities_z", [](ActuatorFarm& al){ return arr({al.getNumberOfTurbines(), al.getNumberOfBladesPerTurbine(), al.getNumberOfPointsPerBlade()}, al.getAllBladeVelocitiesZ()); } )
        .def("get_all_blade_forces_x", [](ActuatorFarm& al){ return arr({al.getNumberOfTurbines(), al.getNumberOfBladesPerTurbine(), al.getNumberOfPointsPerBlade()}, al.getAllBladeForcesX()); } )
        .def("get_all_blade_forces_y", [](ActuatorFarm& al){ return arr({al.getNumberOfTurbines(), al.getNumberOfBladesPerTurbine(), al.getNumberOfPointsPerBlade()}, al.getAllBladeForcesY()); } )
        .def("get_all_blade_forces_z", [](ActuatorFarm& al){ return arr({al.getNumberOfTurbines(), al.getNumberOfBladesPerTurbine(), al.getNumberOfPointsPerBlade()}, al.getAllBladeForcesZ()); } )

        .def("get_all_blade_coords_x_device", [](ActuatorFarm& al) -> intptr_t { return arr_to_cp(al.getAllBladeCoordsXDevice()); } )
        .def("get_all_blade_coords_y_device", [](ActuatorFarm& al) -> intptr_t { return arr_to_cp(al.getAllBladeCoordsYDevice()); } )
        .def("get_all_blade_coords_z_device", [](ActuatorFarm& al) -> intptr_t { return arr_to_cp(al.getAllBladeCoordsZDevice()); } )        
        .def("get_all_blade_velocities_x_device", [](ActuatorFarm& al) -> intptr_t { return arr_to_cp(al.getAllBladeVelocitiesXDevice()); } )
        .def("get_all_blade_velocities_y_device", [](ActuatorFarm& al) -> intptr_t { return arr_to_cp(al.getAllBladeVelocitiesYDevice()); } )
        .def("get_all_blade_velocities_z_device", [](ActuatorFarm& al) -> intptr_t { return arr_to_cp(al.getAllBladeVelocitiesZDevice()); } )
        .def("get_all_blade_forces_x_device", [](ActuatorFarm& al) -> intptr_t { return arr_to_cp(al.getAllBladeForcesXDevice()); } )
        .def("get_all_blade_forces_y_device", [](ActuatorFarm& al) -> intptr_t { return arr_to_cp(al.getAllBladeForcesYDevice()); } )
        .def("get_all_blade_forces_z_device", [](ActuatorFarm& al) -> intptr_t { return arr_to_cp(al.getAllBladeForcesZDevice()); } )

        .def("get_all_hub_coords_x_device", [](ActuatorFarm& al) -> intptr_t { return al.getNumberOfHubPoints() > 0 ? arr_to_cp(al.getAllHubCoordsXDevice()) : 0; })
        .def("get_all_hub_coords_y_device", [](ActuatorFarm& al) -> intptr_t { return al.getNumberOfHubPoints() > 0 ? arr_to_cp(al.getAllHubCoordsYDevice()) : 0; })
        .def("get_all_hub_coords_z_device", [](ActuatorFarm& al) -> intptr_t { return al.getNumberOfHubPoints() > 0 ? arr_to_cp(al.getAllHubCoordsZDevice()) : 0; })
        .def("get_all_hub_velocities_x_device", [](ActuatorFarm& al) -> intptr_t { return al.getNumberOfHubPoints() > 0 ? arr_to_cp(al.getAllHubVelocitiesXDevice()) : 0; })
        .def("get_all_hub_velocities_y_device", [](ActuatorFarm& al) -> intptr_t { return al.getNumberOfHubPoints() > 0 ? arr_to_cp(al.getAllHubVelocitiesYDevice()) : 0; })
        .def("get_all_hub_velocities_z_device", [](ActuatorFarm& al) -> intptr_t { return al.getNumberOfHubPoints() > 0 ? arr_to_cp(al.getAllHubVelocitiesZDevice()) : 0; })
        .def("get_all_hub_forces_x_device", [](ActuatorFarm& al) -> intptr_t { return al.getNumberOfHubPoints() > 0 ? arr_to_cp(al.getAllHubForcesXDevice()) : 0; })
        .def("get_all_hub_forces_y_device", [](ActuatorFarm& al) -> intptr_t { return al.getNumberOfHubPoints() > 0 ? arr_to_cp(al.getAllHubForcesYDevice()) : 0; })
        .def("get_all_hub_forces_z_device", [](ActuatorFarm& al) -> intptr_t { return al.getNumberOfHubPoints() > 0 ? arr_to_cp(al.getAllHubForcesZDevice()) : 0; })

        .def("get_all_tower_coords_x_device", [](ActuatorFarm& al) -> intptr_t { return al.getNumberOfTowerPoints() > 0 ? arr_to_cp(al.getAllTowerCoordsXDevice()) : 0; })
        .def("get_all_tower_coords_y_device", [](ActuatorFarm& al) -> intptr_t { return al.getNumberOfTowerPoints() > 0 ? arr_to_cp(al.getAllTowerCoordsYDevice()) : 0; })
        .def("get_all_tower_coords_z_device", [](ActuatorFarm& al) -> intptr_t { return al.getNumberOfTowerPoints() > 0 ? arr_to_cp(al.getAllTowerCoordsZDevice()) : 0; })
        .def("get_all_tower_velocities_x_device", [](ActuatorFarm& al) -> intptr_t { return al.getNumberOfTowerPoints() > 0 ? arr_to_cp(al.getAllTowerVelocitiesXDevice()) : 0; })
        .def("get_all_tower_velocities_y_device", [](ActuatorFarm& al) -> intptr_t { return al.getNumberOfTowerPoints() > 0 ? arr_to_cp(al.getAllTowerVelocitiesYDevice()) : 0; })
        .def("get_all_tower_velocities_z_device", [](ActuatorFarm& al) -> intptr_t { return al.getNumberOfTowerPoints() > 0 ? arr_to_cp(al.getAllTowerVelocitiesZDevice()) : 0; })
        .def("get_all_tower_forces_x_device", [](ActuatorFarm& al) -> intptr_t { return al.getNumberOfTowerPoints() > 0 ? arr_to_cp(al.getAllTowerForcesXDevice()) : 0; })
        .def("get_all_tower_forces_y_device", [](ActuatorFarm& al) -> intptr_t { return al.getNumberOfTowerPoints() > 0 ? arr_to_cp(al.getAllTowerForcesYDevice()) : 0; })
        .def("get_all_tower_forces_z_device", [](ActuatorFarm& al) -> intptr_t { return al.getNumberOfTowerPoints() > 0 ? arr_to_cp(al.getAllTowerForcesZDevice()) : 0; })

        .def("set_all_blade_coords", [](ActuatorFarm& al, arr coordsX, arr coordsY, arr coordsZ){ 
            al.setAllBladeCoords(np_to_arr(coordsX), np_to_arr(coordsY), np_to_arr(coordsZ)); 
        }, py::arg("blade_coords_x"), py::arg("blade_coords_y"), py::arg("blade_coords_z") )
        .def("set_all_blade_velocities", [](ActuatorFarm& al, arr velocitiesX, arr velocitiesY, arr velocitiesZ){ 
            al.setAllBladeVelocities(np_to_arr(velocitiesX), np_to_arr(velocitiesY), np_to_arr(velocitiesZ)); 
        }, py::arg("blade_velocities_x"), py::arg("blade_velocities_y"), py::arg("blade_velocities_z") )
        .def("set_all_blade_forces", [](ActuatorFarm& al, arr forcesX, arr forcesY, arr forcesZ){ 
            al.setAllBladeForces(np_to_arr(forcesX), np_to_arr(forcesY), np_to_arr(forcesZ));
        }, py::arg("blade_forces_x"), py::arg("blade_forces_y"), py::arg("blade_forces_z") )     
        .def("set_turbine_blade_coords", [](ActuatorFarm& al, uint turbine, arr coordsX, arr coordsY, arr coordsZ){ 
            al.setTurbineBladeCoords(turbine, np_to_arr(coordsX), np_to_arr(coordsY), np_to_arr(coordsZ)); 
        }, py::arg("turbine"), py::arg("blade_coords_x"), py::arg("blade_coords_y"), py::arg("blade_coords_z") )
        .def("set_turbine_blade_velocities", [](ActuatorFarm& al, uint turbine, arr velocitiesX, arr velocitiesY, arr velocitiesZ){
            al.setTurbineBladeVelocities(turbine, np_to_arr(velocitiesX), np_to_arr(velocitiesY), np_to_arr(velocitiesZ)); 
        }, py::arg("turbine"), py::arg("blade_velocities_x"), py::arg("blade_velocities_y"), py::arg("blade_velocities_z") )
        .def("set_turbine_blade_forces", [](ActuatorFarm& al, uint turbine, arr forcesX, arr forcesY, arr forcesZ){ 
            al.setTurbineBladeForces(turbine, np_to_arr(forcesX), np_to_arr(forcesY), np_to_arr(forcesZ)); 
        }, py::arg("turbine"), py::arg("blade_forces_x"), py::arg("blade_forces_y"), py::arg("blade_forces_z") )
        .def("update_forces_and_coordinates", &ActuatorFarm::updateForcesAndCoordinates, py::arg("time"), py::arg("deltaT"))
        .def("enable_output", &ActuatorFarm::enableOutput, py::arg("output_name"), py::arg("t_start_out"), py::arg("t_out"))
        .def("set_turbine_azimuth", &ActuatorFarm::setTurbineAzimuth, py::arg("turbine"), py::arg("azimuth"));

        py::class_<ActuatorFarmStandalone, ActuatorFarm, std::shared_ptr<ActuatorFarmStandalone>>(parentModule, "ActuatorFarmStandalone")
        .def(py::init<  SPtr<Parameter>,
                        SPtr<CudaMemoryManager>,
                        const real,
                        const uint,
                        const std::vector<real>&,
                        const std::vector<real>&,
                        const std::vector<real>&,
                        const std::vector<real>&,
                        const real,
                        const int,
                        const std::optional<ActuatorFarm::HubConfig>&,
                        const std::optional<ActuatorFarm::TowerConfig>&,
                        const std::optional<real>,
                        const std::optional<real>,
                        const std::optional<real>,
                        const uint,
                        const std::optional<std::vector<real>>&>(),
                        py::arg("para"),
                        py::arg("cuda_memory_manager"),
                        py::arg("diameter"),
                        py::arg("number_of_points_per_blade"),
                        py::arg("turbine_positions_x"),
                        py::arg("turbine_positions_y"),
                        py::arg("turbine_positions_z"),
                        py::arg("rotor_speeds"),
                        py::arg("smearing_width"),
                        py::arg("level"),
                        py::arg("hub_config") = py::none(),
                        py::arg("tower_config") = py::none(),
                        py::arg("hub_drag_coeff") = py::none(),
                        py::arg("hub_skin_friction_coeff") = py::none(),
                        py::arg("tower_drag_coeff") = py::none(),
                        py::arg("number_of_blades") = 3,
                        py::arg("blade_normal_coefficients") = py::none());

        py::class_<ActuatorFarmStandaloneVAWT, ActuatorFarm, std::shared_ptr<ActuatorFarmStandaloneVAWT>>(parentModule, "ActuatorFarmStandaloneVAWT")
        .def(py::init<  SPtr<Parameter>,
                         SPtr<CudaMemoryManager>,
                         real,
                         uint,
                         uint,
                         real,
                         const std::vector<real>&,
                         const std::vector<real>&,
                         const std::vector<real>&,
                         const std::vector<real>&,
                         real,
                         int,
                         const std::vector<real>&,
                         const std::vector<real>&,
                         const std::vector<real>&,
                         real,
                         real,
                         real,
                         real,
                         bool,
                         bool,
                         bool
                         >(),
                        py::arg("para"),
                        py::arg("cuda_memory_manager"),
                        py::arg("diameter"),
                        py::arg("number_of_blades"),
                        py::arg("number_of_points_per_blade"),
                        py::arg("rotorHeight"),
                        py::arg("turbine_positions_x"),
                        py::arg("turbine_positions_y"),
                        py::arg("turbine_positions_z"),
                        py::arg("rotor_speeds"),
                        py::arg("smearing_width"),
                        py::arg("level"),
                        py::arg("polar_angle_of_attack_deg"),
                        py::arg("polar_lift_coefficient"),
                        py::arg("polar_drag_coefficient"),
                        py::arg("chord"),
                        py::arg("pitch"),
                        py::arg("blade_mounting_point"),
                        py::arg("velocity_inlet"),
                        py::arg("flag_localized_smearing_width") = false,
                        py::arg("flag_flow_curvature") = false,
                        py::arg("flag_end_effects")=false);
    }
}

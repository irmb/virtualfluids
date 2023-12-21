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
#include <gpu/core/PreCollisionInteractor/Actuator/ActuatorFarm.h>
#include <gpu/core/PreCollisionInteractor/Actuator/ActuatorFarmStandalone.h>
#include <gpu/core/PreCollisionInteractor/PreCollisionInteractor.h>


class PyActuatorFarm : public ActuatorFarm 
{
public:
    using ActuatorFarm::ActuatorFarm; // Inherit constructors
    void updateForcesAndCoordinates() override 
    { 
        PYBIND11_OVERRIDE_PURE_NAME(void, ActuatorFarm, "update_forces_and_coordinates", updateForcesAndCoordinates); 
    }
};

namespace actuator_farm
{
    namespace py = pybind11;

    template<class dtype>
    dtype* np_to_arr(py::array_t<dtype, py::array::c_style> array){ return static_cast<dtype *>(array.request().ptr); };

    template<class dtype>
    intptr_t arr_to_cp(dtype* array){ return reinterpret_cast<intptr_t>(array); };

    void makeModule(py::module_ &parentModule)
    {
        using arr = py::array_t<real, py::array::c_style>;
        
        py::class_<ActuatorFarm, PreCollisionInteractor, PyActuatorFarm, std::shared_ptr<ActuatorFarm>>(parentModule, "ActuatorFarm", py::dynamic_attr())
        .def(py::init<  const real,
                        const std::vector<real>,
                        const std::vector<real>,
                        const std::vector<real>,
                        const std::vector<real>,
                        const real,
                        const real,
                        const int,
                        const real,
                        const real,
                        const bool>(), 
                        py::arg("diameter"),
                        py::arg("blade_radii"),
                        py::arg("turbine_positions_x"),
                        py::arg("turbine_positions_y"),
                        py::arg("turbine_positions_z"),
                        py::arg("density"), 
                        py::arg("smearing_width"),
                        py::arg("level"), 
                        py::arg("delta_t"), 
                        py::arg("delta_x"),
                        py::arg("use_host_arrays"))
        .def_property_readonly("number_of_nodes_per_blade", &ActuatorFarm::getNumberOfNodesPerBlade)
        .def_property_readonly("number_of_turbines", &ActuatorFarm::getNumberOfTurbines)
        .def_property_readonly("number_of_nodes_per_turbine", &ActuatorFarm::getNumberOfNodesPerTurbine)
        .def_property_readonly("number_of_blades_per_turbine", &ActuatorFarm::getNumberOfBladesPerTurbine)
        .def_property_readonly("number_of_grid_nodes", &ActuatorFarm::getNumberOfGridNodes)
        .def_property_readonly("number_of_indices", &ActuatorFarm::getNumberOfIndices)
        .def_property_readonly("density", &ActuatorFarm::getDensity)
        .def_property_readonly("delta_t", &ActuatorFarm::getDeltaT)
        .def_property_readonly("delta_x", &ActuatorFarm::getDeltaX)

        .def("get_turbine_pos", [](ActuatorFarm& al, uint turbine){ 
            real position[3] = {al.getTurbinePosX(turbine), al.getTurbinePosY(turbine), al.getTurbinePosZ(turbine)}; return arr(3,  position);
            }, py::arg("turbine"))
        .def("get_all_turbine_pos_x", [](ActuatorFarm& al){ return arr(al.getNumberOfTurbines(), al.getAllTurbinePosX()); } )
        .def("get_all_turbine_pos_y", [](ActuatorFarm& al){ return arr(al.getNumberOfTurbines(), al.getAllTurbinePosY()); } )
        .def("get_all_turbine_pos_z", [](ActuatorFarm& al){ return arr(al.getNumberOfTurbines(), al.getAllTurbinePosZ()); } )
    
        .def("get_all_blade_coords_x", [](ActuatorFarm& al){ return arr({al.getNumberOfTurbines(), al.getNumberOfBladesPerTurbine(), al.getNumberOfNodesPerBlade()}, al.getAllBladeCoordsX()); } )
        .def("get_all_blade_coords_y", [](ActuatorFarm& al){ return arr({al.getNumberOfTurbines(), al.getNumberOfBladesPerTurbine(), al.getNumberOfNodesPerBlade()}, al.getAllBladeCoordsY()); } )
        .def("get_all_blade_coords_z", [](ActuatorFarm& al){ return arr({al.getNumberOfTurbines(), al.getNumberOfBladesPerTurbine(), al.getNumberOfNodesPerBlade()}, al.getAllBladeCoordsZ()); } )        
        .def("get_all_blade_velocities_x", [](ActuatorFarm& al){ return arr({al.getNumberOfTurbines(), al.getNumberOfBladesPerTurbine(), al.getNumberOfNodesPerBlade()}, al.getAllBladeVelocitiesX()); } )
        .def("get_all_blade_velocities_y", [](ActuatorFarm& al){ return arr({al.getNumberOfTurbines(), al.getNumberOfBladesPerTurbine(), al.getNumberOfNodesPerBlade()}, al.getAllBladeVelocitiesY()); } )
        .def("get_all_blade_velocities_z", [](ActuatorFarm& al){ return arr({al.getNumberOfTurbines(), al.getNumberOfBladesPerTurbine(), al.getNumberOfNodesPerBlade()}, al.getAllBladeVelocitiesZ()); } )
        .def("get_all_blade_forces_x", [](ActuatorFarm& al){ return arr({al.getNumberOfTurbines(), al.getNumberOfBladesPerTurbine(), al.getNumberOfNodesPerBlade()}, al.getAllBladeForcesX()); } )
        .def("get_all_blade_forces_y", [](ActuatorFarm& al){ return arr({al.getNumberOfTurbines(), al.getNumberOfBladesPerTurbine(), al.getNumberOfNodesPerBlade()}, al.getAllBladeForcesY()); } )
        .def("get_all_blade_forces_z", [](ActuatorFarm& al){ return arr({al.getNumberOfTurbines(), al.getNumberOfBladesPerTurbine(), al.getNumberOfNodesPerBlade()}, al.getAllBladeForcesZ()); } )

        .def("get_turbine_blade_coords_x", [](ActuatorFarm& al, uint turbine){ return arr({al.getNumberOfBladesPerTurbine(), al.getNumberOfNodesPerBlade()}, al.getTurbineBladeCoordsXDevice(turbine)); }, py::arg("turbine") )
        .def("get_turbine_blade_coords_y", [](ActuatorFarm& al, uint turbine){ return arr({al.getNumberOfBladesPerTurbine(), al.getNumberOfNodesPerBlade()}, al.getTurbineBladeCoordsYDevice(turbine)); }, py::arg("turbine") )
        .def("get_turbine_blade_coords_z", [](ActuatorFarm& al, uint turbine){ return arr({al.getNumberOfBladesPerTurbine(), al.getNumberOfNodesPerBlade()}, al.getTurbineBladeCoordsZDevice(turbine)); }, py::arg("turbine") )        
        .def("get_turbine_blade_velocities_x", [](ActuatorFarm& al, uint turbine){ return arr({al.getNumberOfBladesPerTurbine(), al.getNumberOfNodesPerBlade()}, al.getTurbineBladeVelocitiesXDevice(turbine)); }, py::arg("turbine") )
        .def("get_turbine_blade_velocities_y", [](ActuatorFarm& al, uint turbine){ return arr({al.getNumberOfBladesPerTurbine(), al.getNumberOfNodesPerBlade()}, al.getTurbineBladeVelocitiesYDevice(turbine)); }, py::arg("turbine") )
        .def("get_turbine_blade_velocities_z", [](ActuatorFarm& al, uint turbine){ return arr({al.getNumberOfBladesPerTurbine(), al.getNumberOfNodesPerBlade()}, al.getTurbineBladeVelocitiesZDevice(turbine)); }, py::arg("turbine") )
        .def("get_turbine_blade_forces_x", [](ActuatorFarm& al, uint turbine){ return arr({al.getNumberOfBladesPerTurbine(), al.getNumberOfNodesPerBlade()}, al.getTurbineBladeForcesXDevice(turbine)); }, py::arg("turbine") )
        .def("get_turbine_blade_forces_y", [](ActuatorFarm& al, uint turbine){ return arr({al.getNumberOfBladesPerTurbine(), al.getNumberOfNodesPerBlade()}, al.getTurbineBladeForcesYDevice(turbine)); }, py::arg("turbine") )
        .def("get_turbine_blade_forces_z", [](ActuatorFarm& al, uint turbine){ return arr({al.getNumberOfBladesPerTurbine(), al.getNumberOfNodesPerBlade()}, al.getTurbineBladeForcesZDevice(turbine)); }, py::arg("turbine") )

        .def("get_all_blade_coords_x_device", [](ActuatorFarm& al) -> intptr_t { return arr_to_cp(al.getAllBladeCoordsXDevice()); } )
        .def("get_all_blade_coords_y_device", [](ActuatorFarm& al) -> intptr_t { return arr_to_cp(al.getAllBladeCoordsYDevice()); } )
        .def("get_all_blade_coords_z_device", [](ActuatorFarm& al) -> intptr_t { return arr_to_cp(al.getAllBladeCoordsZDevice()); } )        
        .def("get_all_blade_velocities_x_device", [](ActuatorFarm& al) -> intptr_t { return arr_to_cp(al.getAllBladeVelocitiesXDevice()); } )
        .def("get_all_blade_velocities_y_device", [](ActuatorFarm& al) -> intptr_t { return arr_to_cp(al.getAllBladeVelocitiesYDevice()); } )
        .def("get_all_blade_velocities_z_device", [](ActuatorFarm& al) -> intptr_t { return arr_to_cp(al.getAllBladeVelocitiesZDevice()); } )
        .def("get_all_blade_forces_x_device", [](ActuatorFarm& al) -> intptr_t { return arr_to_cp(al.getAllBladeForcesXDevice()); } )
        .def("get_all_blade_forces_y_device", [](ActuatorFarm& al) -> intptr_t { return arr_to_cp(al.getAllBladeForcesYDevice()); } )
        .def("get_all_blade_forces_z_device", [](ActuatorFarm& al) -> intptr_t { return arr_to_cp(al.getAllBladeForcesZDevice()); } )

        .def("get_turbine_blade_coords_x_device", [](ActuatorFarm& al, uint turbine) -> intptr_t { return arr_to_cp(al.getTurbineBladeCoordsXDevice(turbine)); }, py::arg("turbine") )
        .def("get_turbine_blade_coords_y_device", [](ActuatorFarm& al, uint turbine) -> intptr_t { return arr_to_cp(al.getTurbineBladeCoordsYDevice(turbine)); }, py::arg("turbine") )
        .def("get_turbine_blade_coords_z_device", [](ActuatorFarm& al, uint turbine) -> intptr_t { return arr_to_cp(al.getTurbineBladeCoordsZDevice(turbine)); }, py::arg("turbine") )        
        .def("get_turbine_blade_velocities_x_device", [](ActuatorFarm& al, uint turbine) -> intptr_t { return arr_to_cp(al.getTurbineBladeVelocitiesXDevice(turbine)); }, py::arg("turbine") )
        .def("get_turbine_blade_velocities_y_device", [](ActuatorFarm& al, uint turbine) -> intptr_t { return arr_to_cp(al.getTurbineBladeVelocitiesYDevice(turbine)); }, py::arg("turbine") )
        .def("get_turbine_blade_velocities_z_device", [](ActuatorFarm& al, uint turbine) -> intptr_t { return arr_to_cp(al.getTurbineBladeVelocitiesZDevice(turbine)); }, py::arg("turbine") )
        .def("get_turbine_blade_forces_x_device", [](ActuatorFarm& al, uint turbine) -> intptr_t { return arr_to_cp(al.getTurbineBladeForcesXDevice(turbine)); }, py::arg("turbine") )
        .def("get_turbine_blade_forces_y_device", [](ActuatorFarm& al, uint turbine) -> intptr_t { return arr_to_cp(al.getTurbineBladeForcesYDevice(turbine)); }, py::arg("turbine") )
        .def("get_turbine_blade_forces_z_device", [](ActuatorFarm& al, uint turbine) -> intptr_t { return arr_to_cp(al.getTurbineBladeForcesZDevice(turbine)); }, py::arg("turbine") )

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
        .def("update_forces_and_coordinates", &ActuatorFarm::updateForcesAndCoordinates)
        .def("enable_output", &ActuatorFarm::enableOutput, py::arg("output_name"), py::arg("t_start_out"), py::arg("t_out"))
        .def("set_turbine_azimuth", &ActuatorFarm::setTurbineAzimuth, py::arg("turbine"), py::arg("azimuth"));

        py::class_<ActuatorFarmStandalone, ActuatorFarm, std::shared_ptr<ActuatorFarmStandalone>>(parentModule, "ActuatorFarmStandalone")
        .def(py::init<  const real,
                        const uint,
                        const std::vector<real>,
                        const std::vector<real>,
                        const std::vector<real>,
                        const std::vector<real>,
                        const real,
                        const real,
                        const int,
                        const real,
                        const real>(), 
                        py::arg("diameter"),
                        py::arg("number_of_nodes_per_blade"), 
                        py::arg("turbine_positions_x"),
                        py::arg("turbine_positions_y"),
                        py::arg("turbine_positions_z"),
                        py::arg("rotor_speeds"),
                        py::arg("density"), 
                        py::arg("smearing_width"),
                        py::arg("level"), 
                        py::arg("delta_t"), 
                        py::arg("delta_x"));
    }


}
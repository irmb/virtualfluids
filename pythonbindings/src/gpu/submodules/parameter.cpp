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
#include <pybind11/functional.h>
#include <pybind11/stl.h>
#include <gpu/core/Parameter/Parameter.h>
#include "basics/constants/NumericConstants.h"
#include <basics/config/ConfigurationFile.h>
#include <gpu/core/PreCollisionInteractor/PreCollisionInteractor.h>


using namespace vf::basics::constant;

namespace parameter
{
    namespace py = pybind11;

    void makeModule(py::module_ &parentModule)
    {
        py::class_<Parameter, std::shared_ptr<Parameter>>(parentModule, "Parameter")
        .def(py::init<
                int,
                int,
                std::optional<const vf::basics::ConfigurationFile*>>(),
                py::arg("number_of_processes"),
                py::arg("my_ID"),
                py::arg("config_data"))
        .def(py::init<int, int>(),
                py::arg("number_of_processes"),
                py::arg("my_ID"))
        .def(py::init<const vf::basics::ConfigurationFile*>(), py::arg("config_data"))
        .def("set_forcing", &Parameter::setForcing, py::arg("forcing_x"), py::arg("forcing_y"), py::arg("forcing_z"))
        .def("set_quadric_limiters", &Parameter::setQuadricLimiters, py::arg("quadric_limiter_p"), py::arg("quadric_limiter_m"), py::arg("quadric_limiter_d"))
        .def("set_diff_on", &Parameter::setDiffOn, py::arg("is_diff"))
        .def("set_max_level", &Parameter::setMaxLevel, py::arg("number_of_levels"))
        .def("set_timestep_end", &Parameter::setTimestepEnd, py::arg("tend"))
        .def("set_timestep_out", &Parameter::setTimestepOut, py::arg("tout"))
        .def("set_timestep_start_out", &Parameter::setTimestepStartOut, py::arg("t_start_out"))
        .def("set_timestep_of_coarse_level", &Parameter::setTimestepOfCoarseLevel, py::arg("timestep"))
        .def("set_calc_turbulence_intensity", &Parameter::setCalcTurbulenceIntensity, py::arg("calc_velocity_and_fluctuations"))
        .def("set_output_path", &Parameter::setOutputPath, py::arg("o_path"))
        .def("set_output_prefix", &Parameter::setOutputPrefix, py::arg("o_prefix"))
        .def("set_print_files", &Parameter::setPrintFiles, py::arg("print_files"))
        .def("set_concentration_init", &Parameter::setConcentrationInit, py::arg("concentration"))
        .def("set_concentration_BC", &Parameter::setConcentrationBC, py::arg("concentrationBC"))
        .def("set_viscosity_LB", &Parameter::setViscosityLB, py::arg("viscosity"))
        .def("set_velocity_LB", &Parameter::setVelocityLB, py::arg("velocity"))
        .def("set_viscosity_ratio", &Parameter::setViscosityRatio, py::arg("viscosity_ratio"))
        .def("set_velocity_ratio", &Parameter::setVelocityRatio, py::arg("velocity_ratio"))
        .def("set_density_ratio", &Parameter::setDensityRatio, py::arg("density_ratio"))
        .def("set_devices", &Parameter::setDevices, py::arg("devices"))
        .def("set_max_dev", &Parameter::setMaxDev, py::arg("max_dev"))
        .def("set_is_body_force", &Parameter::setIsBodyForce, py::arg("is_body_force"))
        .def("set_use_streams", &Parameter::setUseStreams, py::arg("use_streams"))
        .def("configure_main_kernel", &Parameter::configureMainKernel, py::arg("kernel"))
        .def("set_AD_kernel", &Parameter::setADKernel, py::arg("ad_kernel"))
        .def("set_has_wall_model_monitor", &Parameter::setHasWallModelMonitor, py::arg("has_wall_monitor"))
        .def("set_outflow_pressure_correction_factor", &Parameter::setOutflowPressureCorrectionFactor, py::arg("correction_factor"))
        .def("set_initial_condition", [](Parameter &para, std::function<std::vector<float>(real, real, real)> &init_func)
        {
            para.setInitialCondition([init_func](real coordX, real coordY, real coordZ, real& rho, real& vx, real& vy, real& vz)
            {   
                std::vector<float> values = init_func(coordX, coordY, coordZ);
                rho = values[0];
                vx = values[1];
                vy = values[2];
                vz = values[3];
            });
        }, py::arg("init_func"))
        .def("set_initial_condition_uniform", [](Parameter &para, real velocity_x, real velocity_y, real velocity_z)
        {
            para.setInitialCondition([velocity_x, velocity_y, velocity_z](real coordX, real coordY, real coordZ, real& rho, real& vx, real& vy, real& vz) // must capture values explicitly!
            {
                rho = c0o1;
                vx = velocity_x;
                vy = velocity_y;
                vz = velocity_z;
            });
        }, py::arg("velocity_x"), py::arg("velocity_y"), py::arg("velocity_z"))
        .def("set_initial_condition_log_law", [](Parameter &para, real u_star, real z0, real velocityRatio)
        {
            para.setInitialCondition(
                [u_star, z0, velocityRatio](real coordX, real coordY, real coordZ, real& rho, real& vx, real& vy, real& vz)
                {
                    coordZ = coordZ > c0o1 ? coordZ : c0o1;

                    rho = c0o1;
                    vx  = u_star/c4o10 * log(coordZ/z0+c1o1) / velocityRatio;
                    vy = c0o1;
                    vz = c0o1;
                }
            );
        }, py::arg("u_star"), py::arg("z0"), py::arg("velocity_ratio"))
        .def("set_initial_condition_perturbed_log_law", [](Parameter &para, real u_star, real z0, real L_x, real L_z, real H, real velocityRatio)
        {
            para.setInitialCondition(
                [u_star, z0, L_x, L_z, H, velocityRatio](real coordX, real coordY, real coordZ, real& rho, real& vx, real& vy, real& vz)
                {
                    coordZ = coordZ > c0o1 ? coordZ : c0o1;
                    rho = c0o1;
                    vx  = (u_star/c4o10 * log(coordZ/z0+c1o1) + c2o1*sin(cPi*c16o1*coordX/L_x)*sin(cPi*c8o1*coordZ/H)/(pow(coordZ/H,c2o1)+c1o1)) / velocityRatio; 
                    vy  = c2o1*sin(cPi*c16o1*coordX/L_x)*sin(cPi*c8o1*coordZ/H)/(pow(coordZ/H,c2o1)+c1o1) / velocityRatio; 
                    vz  = c8o1*u_star/c4o10*(sin(cPi*c8o1*coordY/H)*sin(cPi*c8o1*coordZ/H)+sin(cPi*c8o1*coordX/L_x))/(pow(c1o2*L_z-coordZ, c2o1)+c1o1) / velocityRatio;
                }
            );
        }, py::arg("u_star"), py::arg("z0"), py::arg("length_x"), py::arg("length_z"), py::arg("height"), py::arg("velocity_ratio"))
        .def("add_actuator", &Parameter::addActuator, py::arg("actuator"))
        .def("add_probe", &Parameter::addProbe, py::arg("probe"))
        .def("get_output_path", &Parameter::getOutputPath)
        .def("get_output_prefix", &Parameter::getOutputPrefix)
        .def("get_velocity", &Parameter::getVelocity)
        .def("get_viscosity", &Parameter::getViscosity)
        .def("get_velocity_ratio", &Parameter::getVelocityRatio)
        .def("get_viscosity_ratio", &Parameter::getViscosityRatio)
        .def("get_density_ratio", &Parameter::getDensityRatio)
        .def("get_force_ratio", &Parameter::getForceRatio)
        .def("get_SGS_constant", &Parameter::getSGSConstant)
        .def("get_is_body_force", &Parameter::getIsBodyForce)
        ;

    }
}
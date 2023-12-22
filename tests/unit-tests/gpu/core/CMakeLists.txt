#=======================================================================================
# ____          ____    __    ______     __________   __      __       __        __
# \    \       |    |  |  |  |   _   \  |___    ___| |  |    |  |     /  \      |  |
#  \    \      |    |  |  |  |  |_)   |     |  |     |  |    |  |    /    \     |  |
#   \    \     |    |  |  |  |   _   /      |  |     |  |    |  |   /  /\  \    |  |
#    \    \    |    |  |  |  |  | \  \      |  |     |   \__/   |  /  ____  \   |  |____
#     \    \   |    |  |__|  |__|  \__\     |__|      \________/  /__/    \__\  |_______|
#      \    \  |    |   ________________________________________________________________
#       \    \ |    |  |  ______________________________________________________________|
#        \    \|    |  |  |         __          __     __     __     ______      _______
#         \         |  |  |_____   |  |        |  |   |  |   |  |   |   _  \    /  _____)
#          \        |  |   _____|  |  |        |  |   |  |   |  |   |  | \  \   \_______
#           \       |  |  |        |  |_____   |   \_/   |   |  |   |  |_/  /    _____  |
#            \ _____|  |__|        |________|   \_______/    |__|   |______/    (_______/
#
#  This file is part of VirtualFluids. VirtualFluids is free software: you can
#  redistribute it and/or modify it under the terms of the GNU General Public
#  License as published by the Free Software Foundation, either version 3 of
#  the License, or (at your option) any later version.
#
#  VirtualFluids is distributed in the hope that it will be useful, but WITHOUT
#  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
#  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
#  for more details.
#
#  You should have received a copy of the GNU General Public License along
#  with VirtualFluids (see COPYING.txt). If not, see <http://www.gnu.org/licenses/>.
#
#  SPDX-License-Identifier: GPL-3.0-or-later
#  SPDX-FileCopyrightText: Copyright © VirtualFluids Project contributors, see AUTHORS.md in root folder
#
#! \author Soeren Peters
#=======================================================================================
project(gpu_coreTests LANGUAGES CXX)

vf_add_tests(NAME gpu_coreTests PRIVATE_LINK gpu_core project_options project_warnings)

# set_target_properties(gpu_coreTests PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
target_include_directories(gpu_coreTests PRIVATE "${VF_THIRD_DIR}/cuda_samples/")
target_include_directories(gpu_coreTests PRIVATE "${VF_ROOT_DIR}/src/gpu/GridGenerator/")
set_source_files_properties(Communication/ExchangeData27Test.cpp PROPERTIES LANGUAGE CUDA)
set_source_files_properties(BoundaryConditions/BoundaryConditionFactoryTest.cpp PROPERTIES LANGUAGE CUDA)
set_source_files_properties(DataStructureInitializer/GridReaderGenerator/GridGeneratorTest.cpp PROPERTIES LANGUAGE CUDA)
set_source_files_properties(DataStructureInitializer/GridReaderGenerator/IndexRearrangementForStreamsTest.cpp PROPERTIES LANGUAGE CUDA)
set_source_files_properties(Kernel/Kernels/BasicKernels/FluidFlow/Compressible/CumulantK17/CumulantK17Test.cpp PROPERTIES LANGUAGE CUDA)
set_source_files_properties(Utilities/KernelUtilitiesTests.cpp PROPERTIES LANGUAGE CUDA)
set_source_files_properties(Output/DistributionDebugWriterTest.cpp PROPERTIES LANGUAGE CUDA)
set_source_files_properties(Parameter/ParameterTest.cpp PROPERTIES LANGUAGE CUDA)
set_source_files_properties(PreCollisionInteractor/Actuator/ActuatorFarmInlinesTest.cpp PROPERTIES LANGUAGE CUDA)

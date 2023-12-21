#  SPDX-License-Identifier: CC-BY-4.0
#  SPDX-FileCopyrightText: Copyright © VirtualFluids Project contributors, see AUTHORS.md in root folder
SET(CMAKE_CUDA_ARCHITECTURES "75")

list(APPEND USER_APPS "apps/gpu/ActuatorLine")
list(APPEND USER_APPS "apps/gpu/SphereMultiGPU")

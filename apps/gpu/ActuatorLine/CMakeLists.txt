PROJECT(ActuatorLine LANGUAGES CUDA CXX)

vf_add_library(BUILDTYPE binary PRIVATE_LINK basics gpu_core GridGenerator FILES ActuatorLine.cpp)

set_source_files_properties(ActuatorLine.cpp PROPERTIES LANGUAGE CUDA)

set_target_properties(ActuatorLine PROPERTIES CUDA_SEPARABLE_COMPILATION ON)

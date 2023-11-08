PROJECT(ChannelFlow LANGUAGES CUDA CXX)

vf_add_library(BUILDTYPE binary PRIVATE_LINK basics gpu_core GridGenerator FILES ChannelFlow.cpp)

set_source_files_properties(ChannelFlow.cpp PROPERTIES LANGUAGE CUDA)

set_target_properties(ChannelFlow PROPERTIES CUDA_SEPARABLE_COMPILATION ON)

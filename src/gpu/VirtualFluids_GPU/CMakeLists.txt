project(VirtualFluids_GPU LANGUAGES CUDA CXX)

set(additional_libraries "")
if(MSVC)
    set(additional_libraries ws2_32 Traffic) # ws_32 throws an error on Phoenix
endif()

vf_add_library(BUILDTYPE static PRIVATE_LINK ${additional_libraries} GridGenerator basics MPI::MPI_CXX PUBLIC_LINK lbmCuda)

linkBoost(COMPONENTS "serialization")

#SET(TPN_WIN32 "/EHsc")
#https://stackoverflow.com/questions/6832666/lnk2019-when-including-asio-headers-solution-generated-with-cmake
#https://stackoverflow.com/questions/27442885/syntax-error-with-stdnumeric-limitsmax

#set_target_properties(VirtualFluids_GPU PROPERTIES CUDA_RESOLVE_DEVICE_SYMBOLS ON)
set_target_properties(VirtualFluids_GPU PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
#set_target_properties(VirtualFluids_GPU PROPERTIES POSITION_INDEPENDENT_CODE ON)
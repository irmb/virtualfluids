project(VirtualFluids_GPU LANGUAGES CUDA CXX)

set(additional_libraries "")
if(MSVC)
    set(additional_libraries ws2_32 Traffic) # ws_32 throws an error on Phoenix
endif()

vf_add_library(PRIVATE_LINK ${additional_libraries} GridGenerator basics MPI::MPI_CXX)

linkBoost(COMPONENTS "serialization")

#SET(TPN_WIN32 "/EHsc")
#https://stackoverflow.com/questions/6832666/lnk2019-when-including-asio-headers-solution-generated-with-cmake
#https://stackoverflow.com/questions/27442885/syntax-error-with-stdnumeric-limitsmax

vf_add_tests()

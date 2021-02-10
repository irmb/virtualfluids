project(VirtualFluids_GPU LANGUAGES CUDA CXX)

IF(MSVC)
    set(libsToLink ws2_32 GridGenerator basics Traffic) # ws_32 throws an error on Phoenix
ELSE(MSVC)
    #set(libsToLink GridGenerator VirtualFluidsBasics Core Traffic)
    set(libsToLink GridGenerator basics)
ENDIF(MSVC)


vf_add_library(PRIVATE_LINK ${libsToLink})

linkMPI()
linkBoost(COMPONENTS "serialization")

#SET(TPN_WIN32 "/EHsc")
#https://stackoverflow.com/questions/6832666/lnk2019-when-including-asio-headers-solution-generated-with-cmake
#https://stackoverflow.com/questions/27442885/syntax-error-with-stdnumeric-limitsmax

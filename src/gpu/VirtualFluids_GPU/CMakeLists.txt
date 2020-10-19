
IF(MSVC)
    set(libsToLink ws2_32 GridGenerator basics Traffic) # ws_32 throws an error on Phoenix
ELSE(MSVC)
    #set(libsToLink GridGenerator VirtualFluidsBasics Core Traffic)
    set(libsToLink GridGenerator basics)
ENDIF(MSVC)


vf_add_library(BUILDTYPE shared PRIVATE_LINK ${libsToLink})

linkMPI()
linkCUDA()
linkBoost("serialization")

#SET(TPN_WIN32 "/EHsc")
#https://stackoverflow.com/questions/6832666/lnk2019-when-including-asio-headers-solution-generated-with-cmake
#https://stackoverflow.com/questions/27442885/syntax-error-with-stdnumeric-limitsmax

IF(MSVC)
    vf_get_library_name(library_name)
    set_target_properties(${library_name} PROPERTIES LINK_FLAGS "/ignore:4251")
ENDIF(MSVC)

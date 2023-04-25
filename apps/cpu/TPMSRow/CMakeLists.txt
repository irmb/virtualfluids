PROJECT(TPMSRow)

vf_add_library(BUILDTYPE binary PRIVATE_LINK VirtualFluidsCore muparser basics ${MPI_CXX_LIBRARIES} FILES TPMSRow.cpp )

vf_get_library_name (library_name)
#target_include_dires(${library_name} PRIVATE ${APPS_ROOT_CPU})
#target_include_dires(${library_name} PRIVATE "/cluster/lib/boost/1.63.0/gcc")
target_include_directories(${library_name} PRIVATE ${APPS_ROOT_CPU})
#target_include_directories(${library_name} PRIVATE "/cluster/lib/boost/1.63.0/gcc/include/")


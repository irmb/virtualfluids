#############################################################
###                  Options                              ###
#############################################################
SET(VF_CPU_ENABLE_VTK OFF CACHE BOOL "include VTK library support")
SET(VF_CP_ENABLE_CATALYST OFF CACHE BOOL "include Paraview Catalyst support")

IF(${VF_CPU_ENABLE_VTK})
    FIND_PACKAGE(VTK REQUIRED)
    INCLUDE_DIRECTORIES(${VTK_INCLUDE_DIRS})
    target_compile_definitions(project_options INTERFACE VF_VTK)
ENDIF()

IF(${VF_CP_ENABLE_CATALYST})
    find_package(ParaView 4.3 REQUIRED COMPONENTS vtkPVPythonCatalyst)
    include("${PARAVIEW_USE_FILE}")
    target_compile_definitions(project_options INTERFACE VF_CATALYST)
ENDIF()

target_compile_definitions(project_options INTERFACE VF_METIS)

#############################################################
###                  Libraries                            ###
#############################################################

add_subdirectory(${VF_THIRD_DIR}/MuParser)
add_subdirectory(${VF_THIRD_DIR}/metis/metis-5.1.0)

add_subdirectory(src/cpu/core)

if(VF_ENABLE_PYTHON_BINDINGS)
    add_subdirectory(src/cpu/simulationconfig)
endif()

#############################################################
###                  Apps                                 ###
#############################################################

add_subdirectory(apps/cpu/FlowAroundCylinder)
#add_subdirectory(apps/cpu/LidDrivenCavity)
#add_subdirectory(apps/cpu/LaminarPlaneFlow)
#add_subdirectory(apps/cpu/LaminarPipeFlow)
#add_subdirectory(apps/cpu/AcousticPulse)

if(VF_ENABLE_BOOST)
    add_subdirectory(apps/cpu/GyroidsRow)
endif()

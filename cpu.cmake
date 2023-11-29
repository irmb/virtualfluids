#############################################################
###                  Options                              ###
#############################################################
SET(VFCPU_USE_VTK OFF CACHE BOOL "include VTK library support")
SET(VFCPU_USE_CATALYST OFF CACHE BOOL "include Paraview Catalyst support")

IF(${VFCPU_USE_VTK})
    FIND_PACKAGE(VTK REQUIRED)
    INCLUDE_DIRECTORIES(${VTK_INCLUDE_DIRS})
ENDIF()

IF(${VFCPU_USE_CATALYST})
    find_package(ParaView 4.3 REQUIRED COMPONENTS vtkPVPythonCatalyst)
    include("${PARAVIEW_USE_FILE}")
ENDIF()

list(APPEND VF_COMPILER_DEFINITION VF_METIS)

IF(${VFCPU_USE_VTK})
    list(APPEND VF_COMPILER_DEFINITION VF_VTK)
ENDIF()
IF(${VFCPU_USE_CATALYST})
    list(APPEND VF_COMPILER_DEFINITION VF_CATALYST)
ENDIF()

#############################################################
###                  Libraries                            ###
#############################################################

add_subdirectory(${VF_THIRD_DIR}/MuParser)
add_subdirectory(${VF_THIRD_DIR}/metis/metis-5.1.0)

add_subdirectory(src/cpu/core)

if(BUILD_VF_PYTHON_BINDINGS)
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

if(BUILD_USE_BOOST)
    add_subdirectory(apps/cpu/GyroidsRow)
endif()

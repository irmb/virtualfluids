cmake_minimum_required(VERSION 2.4)

INCLUDE("../../../../CMakeCABMacros.txt")
INCLUDE("../../../../CMakeSetCompilerFlags.txt")

CHECK_FOR_VARIABLE(CAB_MACHINE "machine name, e.g. ALTIX, ARWEN")
SET(CMAKE_CONFIG_FILE "${SOURCE_ROOT}/cmake_config_files/${CAB_MACHINE}.config.cmake")

IF(CAB_MACHINE AND EXISTS ${CMAKE_CONFIG_FILE} )
   
   PROJECT(geo3dCellCutTests)

   SET(EXECUTABLE_NAME geo3dCellCutTests)

   #erst hier das config file einf�gen, ansonsten werden manche settings durch (Project) �berschrieben)  
   INCLUDE(${CMAKE_CONFIG_FILE})  
  
   #################################################################
   ###   PACKAGES						###
   #################################################################
   INCLUDE(${SOURCE_ROOT}/basics/utilities/CMakePackage.txt)
   INCLUDE(${SOURCE_ROOT}/basics/objects/CMakePackage.txt)
   INCLUDE(${SOURCE_ROOT}/geometry3d/CMakePackage.txt)
   INCLUDE(${SOURCE_ROOT}/geometry3d/creator/CMakePackage.txt)
   INCLUDE(${SOURCE_ROOT}/geometry3d/fem/CMakePackage.txt)
   
   #################################################################
   ###   OWN DEFINES 						###
   #################################################################
   FILE(GLOB SPECIFIC_FILES ${CMAKE_CURRENT_SOURCE_DIR}/*.h
                            ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp ) 


   SET(ALL_SOURCES ${ALL_SOURCES} ${SPECIFIC_FILES})
   SOURCE_GROUP(example FILES ${SPECIFIC_FILES})

   SET_COMPILER_SPECIFIC_FLAGS(${CAB_COMPILER} BINARY)
   
   #################################################################
   ###   ADDITIONAL_MAKE_CLEAN_FILES                             ###
   #################################################################
   SET_DIRECTORY_PROPERTIES(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES "${GENERATED_FILES}")

   #################################################################
   ###   EXCECUTABLE						###
   #################################################################
   ADD_EXECUTABLE(${EXECUTABLE_NAME} ${ALL_SOURCES} )

   #################################################################
   ###   ADDITIONAL LINK PROPERTIES    			        ###
   #################################################################
   IF(CAB_ADDITIONAL_LINK_FLAGS)
     SET_TARGET_PROPERTIES(${EXECUTABLE_NAME} PROPERTIES LINK_FLAGS ${CAB_ADDITIONAL_LINK_FLAGS})
   ENDIF(CAB_ADDITIONAL_LINK_FLAGS)
ELSE()
    IF(CAB_MACHINE)
      MESSAGE("CAB_MACHINE error - following file is missing: \n ${CMAKE_CONFIG_FILE}")
    ELSE()
      MESSAGE("check CAB_MACHINE!!!")
    ENDIF()  
ENDIF()
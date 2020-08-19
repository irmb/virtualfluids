SET(SOURCE_ROOT $ENV{CAB_DIR} CACHE PATH "(e.g. d:/temp/source)" )
STRING(REGEX REPLACE "\\\\" "/" SOURCE_ROOT ${SOURCE_ROOT}) # "\" --> "/"

IF(EXISTS ${SOURCE_ROOT})
 IF(EXISTS ${SOURCE_ROOT}/basics)
  
  INCLUDE(${SOURCE_ROOT}/CMakeCABMacros.txt) 
 
  PROJECT (STL2INP)
  SET(EXECUTABLE_NAME stl2inp)

  #################################################################
  ###   PACKAGES						###
  #################################################################
 
  INCLUDE(${SOURCE_ROOT}/basics/objects/CMakeLists.txt)
  INCLUDE(${SOURCE_ROOT}/basics/utilities/CMakeLists.txt)
  INCLUDE(${SOURCE_ROOT}/basics/relation/CMakeLists.txt)

 
  INCLUDE(${SOURCE_ROOT}/numerics/geometry2d/CMakeLists.txt)
  INCLUDE(${SOURCE_ROOT}/numerics/geometry2d/creator/CMakeLists.txt)
 
  INCLUDE(${SOURCE_ROOT}/numerics/geometry3d/CMakeLists.txt)
  INCLUDE(${SOURCE_ROOT}/numerics/geometry3d/creator/CMakeLists.txt)
 
  INCLUDE(${SOURCE_ROOT}/octree/CMakeLists.txt)
  INCLUDE(${SOURCE_ROOT}/octree/facette/CMakeLists.txt)
 
  INCLUDE(${SOURCE_ROOT}/quadtree/CMakeLists.txt)
  INCLUDE(${SOURCE_ROOT}/quadtree/nodeadaptation/CMakeLists.txt)

 
  #################################################################
  ###   OWN DEFINES                                             ###
  #################################################################
 # ADD_DEFINITIONS( -DCAB_QT )
 # ADD_DEFINITIONS( -DCAB_QTONLY )
  ADD_DEFINITIONS( -DMEMPOOL_A2PIIOCTNODEVD)
  SET(SPECIFIC_FILES main.cpp stl2inp.h stl2inp.cpp QDefineUniformMesh.h QDefineUniformMesh.cpp)
  SET(ALL_SOURCES ${ALL_SOURCES} ${SPECIFIC_FILES})
  SOURCE_GROUP(z_specific FILES ${SPECIFIC_FILES})

  IF(WIN32)
    ADD_CXX_FLAGS("/wd4996")
  ELSE(WIN32)
    ADD_CXX_FLAGS("-O3 -mcpu=athlon-4 -fomit-frame-pointer -finline-functions -funroll-all-loops")
  ENDIF(WIN32)

  

  #################################################################
  ### QT SPECIFIC (only has effects if a QT source is included)	###
  #################################################################
  #QT specific 
   SET(NEED_QT "YES")
  INCLUDE(${SOURCE_ROOT}/CMakeQtMacros.txt)

  IF(QT_FOUND)
    INCLUDE_DIRECTORIES( ${QT_INCLUDE_DIR} 
			 ${QT_INCLUDE_PATH} 
			 ${QT_QTCORE_INCLUDE_DIR}
			 ${QT_QTGUI_INCLUDE_DIR}  
                       )
    LINK_LIBRARIES ( ${QT_QTCORE_LIBRARY} 
                     ${QT_QTGUI_LIBRARY} 
                    )
    ADD_DEFINITIONS( ${QT_DEFINITIONS})
   ELSE(QT_FOUND)
    IF(${NEED_QT} MATCHES "YES") 
      MESSAGE("Ups\nAt least one package needs Qt!\nPlease check Qt settings\n(e.g. librarys within Advanced Values)")
    ENDIF(${NEED_QT} MATCHES "YES") 
   
    
  ENDIF(QT_FOUND)

 IF(QT_FOUND)
     SET(SUBDIRPATH numerics/geometry3d/examples/stl2inp) 
     SET(CURRENT_DIR ${SOURCE_ROOT}/${SUBDIRPATH})
 
     ################################################################
    ###         Qt4      UI FILES                                ###
    ################################################################
    FILE(GLOB UI_FILES ${CURRENT_DIR}/*.ui)           #collect ui files
    QT4_WRAP_UI(${CURRENT_DIR} OUTFILES ${UI_FILES})  #wrap ui files
    REMOVE(TEMP_FILES ${OUTFILES} )
    SET(TEMP_FILES ${TEMP_FILES} ${OUTFILES} )
    
    #make subfolders for VS with new files
    SOURCE_GROUP(${SUBDIRPATH} FILES ${OUTFILES})
    
    IF(WIN32)
      SET(ALL_SOURCES ${ALL_SOURCES} ${UI_FILES})
      SOURCE_GROUP(${SUBDIRPATH} FILES ${UI_FILES})
    ENDIF(WIN32)
 
   ################################################################
    ###       Qt4        HEADERS TO BE MOCED                     ###
    ################################################################
    MAKE_DIRECTORY(${CURRENT_DIR}${QTGEN_MOC})
    SET(MOC_FILES ) #empty MOC_FILES
    #SET(MOC_CLASSES  ${CURRENT_DIR}/stl2inp.h QDefineUniformMesh.ui )

    QT4_WRAP_CPP(${CURRENT_DIR}${QTGEN_MOC} MOC_FILES ${MOC_CLASSES})
    REMOVE(TEMP_FILES ${MOC_FILES})
    SET(TEMP_FILES ${TEMP_FILES} ${MOC_FILES})
    SOURCE_GROUP(${SUBDIRPATH}${QTGEN_MOC} FILES ${MOC_FILES})
   
    SET(ALL_SOURCES ${ALL_SOURCES} ${TEMP_FILES})
   
  ENDIF(QT_FOUND)


  #################################################################
  ###   ADDITIONAL_MAKE_CLEAN_FILES                             ###
  #################################################################
  SET_DIRECTORY_PROPERTIES(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES "${GENERATED_FILES}")
  
  #################################################################
  ###   EXCECUTABLE						###
  #################################################################
    ADD_EXECUTABLE( ${EXECUTABLE_NAME} 
                  ${ALL_SOURCES}
                )


  #################################################################
  ###   ADDITIONAL LINK PROPERTIES    			        ###
  #################################################################
  IF(CAB_ADDITIONAL_LINK_FLAGS)
    SET_TARGET_PROPERTIES(${EXECUTABLE_NAME} PROPERTIES LINK_FLAGS ${CAB_ADDITIONAL_LINK_FLAGS})
  ENDIF(CAB_ADDITIONAL_LINK_FLAGS)

 ELSE(EXISTS ${SOURCE_ROOT}/basics)
   MESSAGE("Set Path to \"source\" directory at\nSOURCE_ROOT\nis not correct")
 ENDIF(EXISTS ${SOURCE_ROOT}/basics)
ELSE(EXISTS ${SOURCE_ROOT})
 SET(SOURCE_ROOT "CAB_DIR NOT FOUND" CACHE PATH "(e.g. d:/temp/source)" FORCE)
 MESSAGE("Please Set Path to \"source\" directory at\nSOURCE_ROOT\n(e.g. D:/temp/source)")
ENDIF(EXISTS ${SOURCE_ROOT})

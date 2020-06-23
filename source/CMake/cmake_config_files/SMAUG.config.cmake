IF(NOT CMAKE_CXX_COMPILER)
   MESSAGE(FATAL_ERROR "before cmake-config-file can be included -> project must be extecuted")
ENDIF()
  

#################################################################################
#  BOOST  
#################################################################################
SET(BOOST_VERSION "1.44")
SET(BOOST_USE_MULTITHREAD ON)
SET(BOOST_USE_STATIC_LIBS ON)

IF(BOOST_VERSION MATCHES "1.44" )
  SET(BOOST_INCLUDEDIR "C:/boost/boost_1_44_0")
  SET(BOOST_LIBRARYDIR "C:/boost/boost_1_44_0/stage64/lib") 
ENDIF()
#################################################################################
#  ZOLTAN  
#################################################################################
  SET(ZOLTAN_INCLUDEDIR "c:/zoltan/include")
  SET(ZOLTAN_LIBRARY "c:/zoltan/lib/Debug/zoltan.lib") 

set(SOURCE_ROOT "g:/Kostja/pFluid/source") 



  

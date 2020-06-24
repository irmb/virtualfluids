#################################################################################
# COMPILER FLAG VAR
#################################################################################

SET(CAB_COMPILER "gcc41" CACHE STRING "gcc41, gcc43, intel9" FORCE) 
# on MacOS X __APPLE__ is defined but not __unix__
add_definitions(-D__unix__)
SET(CMAKE_OSX_ARCHITECTURES x86_64 CACHE STRING "ppc;i386;ppc64;x86_64" FORCE)

#ADD_CXX_FLAGS(-arch x86_64 -L/usr/local/lib/x86_64)# -Wall -ansi)
#ADD_CXX_BUILDTYPE_FLAGS(DEBUG -DDEBUG -D_DEBUG -gdwarf-2 -arch x86_64 -L/usr/local/lib/x86_64 -Wall -ansi)
#ADD_CXX_BUILDTYPE_FLAGS(RELEASE -O3 -DNDEBUG -fomit-frame-pointer -finline-functions -funroll-all-loops -arch x86_64 -L/usr/local/lib/x86_64 -Wall -ansi)
ADD_CXX_BUILDTYPE_FLAGS(DEBUG -DDEBUG -D_DEBUG -gdwarf-2 -Wall -ansi)
ADD_CXX_BUILDTYPE_FLAGS(RELEASE -O3 -DNDEBUG -fomit-frame-pointer -finline-functions -funroll-all-loops -Wall -ansi)

#ADD_C_FLAGS(-arch x86_64 -L/usr/local/lib/x86_64)# -Wall -ansi)


SET(BOOST_VERSION "1.47.0")
#IF(CAB_COMPILER MATCHES "gcc41")
# SET(BOOST_USE_MULTITHREAD ON)
# SET(BOOST_USE_STATIC_LIBS ON)
#
#   SET(BOOST_COMPILER_SUFFIX -xgcc40)
# 
# IF(BOOST_VERSION MATCHES "1.37.0" )
#		SET(BOOST_INCLUDEDIR "/Users/hg/bin/boost/boost_1_37_0")
#		SET(BOOST_LIBRARYDIR "/Users/hg/bin/boost/boost_1_37_0/lib")
# ENDIF()
## IF(BOOST_VERSION MATCHES "1.37.0" )
##		SET(BOOST_INCLUDEDIR "/opt/local/include")
##		SET(BOOST_LIBRARYDIR "/opt/local/lib")
## ENDIF()
#ENDIF()
#IF(BOOST_VERSION AND NOT BOOST_INCLUDEDIR)
# MESSAGE("${BOOST_VERSION} not found on ${CAB_MACHINE} for specified compiler")
#ENDIF()

SET(BOOST_USE_MULTITHREAD ON)
SET(BOOST_USE_STATIC_LIBS ON)

IF(BOOST_VERSION MATCHES "1.47.0" )
  SET(BOOST_INCLUDEDIR "/Users/hg/bin/boost/boost_1_47_0_bin/include")
  SET(BOOST_LIBRARYDIR "/Users/hg/bin/boost/boost_1_47_0_bin/lib") 
ENDIF()
#################################################################################
#  ZOLTAN  
#################################################################################
IF(${USE_ZOLTAN})
  SET(ZOLTAN_INCLUDEDIR "/Users/hg/bin/zoltan_distrib_v3.501_bin/include")
  SET(ZOLTAN_DEBUG_LIBRARY "/Users/hg/bin/zoltan_distrib_v3.501_bin/lib/libzoltan.a") 
  SET(ZOLTAN_RELEASE_LIBRARY "/Users/hg/bin/zoltan_distrib_v3.501_bin/lib/libzoltan.a") 
ENDIF()
set(SOURCE_ROOT "/Users/hg/hg+/Konstantin_Kucher/Patchcode/src")  
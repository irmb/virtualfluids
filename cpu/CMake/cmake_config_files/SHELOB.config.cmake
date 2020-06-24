LIST(APPEND CAB_ADDTIONAL_COMPILER_FLAGS -D__unix__)
LIST(APPEND CAB_ADDTIONAL_COMPILER_FLAGS -D__UNIX__)
LIST(APPEND CAB_ADDTIONAL_COMPILER_FLAGS -D__APPLE__)

#################################################################################
#  RUBY
#################################################################################
#IF(CMAKE_CL_64)
#   SET(RUBY_EXECUTABLE "C:/ruby191/bin/ruby.exe")
#   SET(RUBY_LIBRARY "C:/ruby191/lib/msvcr90-ruby191.lib")
#ELSE()
#   SET(RUBY_EXECUTABLE "C:/ruby186/bin/ruby.exe")
#   SET(RUBY_LIBRARY "C:/ruby186/lib/msvcrt-ruby18.lib")
#ENDIF()

#################################################################################
#  SWIG
#################################################################################
SET(SWIG_EXECUTABLE "/Users/freud/dev/swig/bin/swig")


########################################################################################
##                            BOOST ALLGMEINGUELTIG                                   ##
########################################################################################
#standard boost: 1.41.0
SET(BOOST_VERSION "1.41.0" CACHE STRING "std: 1.41.0")

SET(MPI_DIR  "/usr/local/mpich2/gcc401")

IF(CAB_COMPILER MATCHES "gcc42")
   SET(BOOST_USE_MULTITHREAD ON)

   SET(BOOST_COMPILER_SUFFIX -xgcc40)
   SET(BOOST_USE_STATIC_LIBS ON)

   IF(BOOST_VERSION MATCHES "1.41.0" )
      SET(BOOST_INCLUDEDIR "/scratch/shared/boost_1_41_0/include")
      SET(BOOST_LIBRARYDIR "/scratch/shared/boost_1_41_0/lib")
   ELSEIF(BOOST_VERSION MATCHES "1.39.0" )
      SET(BOOST_INCLUDEDIR "/usr/local/boost/gcc401/boost_1_39_0/include/boost-1_39")
      SET(BOOST_LIBRARYDIR "/usr/local/boost/gcc401/boost_1_39_0/lib")
   ENDIF()
ENDIF()

IF(BOOST_VERSION AND NOT BOOST_INCLUDEDIR)
  MESSAGE("${BOOST_VERSION} not found on ${CAB_MACHINE} for specified compiler")
ENDIF()



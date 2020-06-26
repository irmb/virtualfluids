#LIST(APPEND CAB_ADDTIONAL_COMPILER_FLAGS -D__unix__)
#LIST(APPEND CAB_ADDTIONAL_COMPILER_FLAGS -D__UNIX__)
LIST(APPEND CAB_ADDTIONAL_COMPILER_FLAGS -D__AIX__)


IF(NOT CMAKE_CXX_COMPILER)
   MESSAGE(FATAL_ERROR "before cmake-config-file can be included -> project must be extecuted")
ENDIF()

#################################################################################
# MPI
#################################################################################
SET(USE_MPI_CXX_SYNTAX OFF)
#################################################################################
#  BOOST
#################################################################################
SET(BOOST_VERSION "1.44")
SET(BOOST_USE_MULTITHREAD ON)
SET(BOOST_USE_STATIC_LIBS ON)

IF(BOOST_VERSION MATCHES "1.44" )
  SET(BOOST_INCLUDEDIR "/sp6/userdeisa/lrz102bk/tools/boost_1_44_0")
  SET(BOOST_LIBRARYDIR "/sp6/userdeisa/lrz102bk/tools/boost_1_44_0/stage/lib")
ENDIF()

set(SOURCE_ROOT "/sp6/userdeisa/lrz102bk/projects/pFluid/source")
set(pFluid_ROOT "/sp6/userdeisa/lrz102bk/projects/pFluid")
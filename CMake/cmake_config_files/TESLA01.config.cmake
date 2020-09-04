#Don't change:
SET(METIS_ROOT ${VF_THIRD_DIR}/metis/metis-5.1.0 CACHE PATH "METIS ROOT")
SET(GMOCK_ROOT ${VF_THIRD_DIR}/googletest CACHE PATH "GMOCK ROOT")
SET(JSONCPP_ROOT ${VF_THIRD_DIR}/jsoncpp CACHE PATH "JSONCPP ROOT")
SET(FFTW_ROOT ${VF_THIRD_DIR}/fftw/fftw-3.3.7 CACHE PATH "JSONCPP ROOT")


#SET TO CORRECT PATH:
SET(BOOST_ROOT  "C:\\Libraries\\boost_1_65_1"  CACHE PATH "BOOST_ROOT")
SET(BOOST_LIBRARYDIR  "C:\\Libraries\\boost_1_65_1\\lib" CACHE PATH "BOOST_LIBRARYDIR")

SET(VTK_DIR "C:/Libraries/VTK-8.0.1/build")

IF(${USE_METIS})
    SET(METIS_INCLUDEDIR "C:/Libraries/metis-5.1.0//include")
    SET(METIS_DEBUG_LIBRARY "C:/Libraries/metis-5.1.0/build/libmetis/Debug/metis.lib")
    SET(METIS_RELEASE_LIBRARY "C:/Libraries/metis-5.1.0/build/libmetis/Release/metis.lib")
ENDIF()

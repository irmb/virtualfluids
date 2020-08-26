#ifndef VIRTUAL_FLUIDS_DEFINITIONS_H
#define VIRTUAL_FLUIDS_DEFINITIONS_H

#cmakedefine BUILD_SHARED_LIBS

// disable warnings 
#pragma warning(disable: 4251)// occurs normally while exporting standard library: "needs to have dll-interface to be"
#pragma warning(disable: 4275) // on dll-interface class <classname> used as base for dll-interface class <classname>


// double or single precision
#cmakedefine VF_DOUBLE_ACCURACY

// External libraries
#cmakedefine VF_BUILD_WITH_MPI
#cmakedefine VF_BUILD_WITH_METIS
#cmakedefine VF_BUILD_WITH_CUDA

// Compiler
#cmakedefine VF_CXX_COMPILER_IS_GNU
#cmakedefine VF_CXX_COMPILER_IS_INTEL
#cmakedefine VF_CXX_COMPILER_IS_IBM
#cmakedefine VF_CXX_COMPILER_IS_MSVC
#cmakedefine VF_CXX_COMPILER_IS_CLANG

#endif

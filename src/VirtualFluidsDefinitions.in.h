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


// Generic helper definitions for shared library support
#if defined _WIN32 || defined __CYGWIN__
  #define VF_SHARED_LIB_IMPORT __declspec(dllimport)
  #define VF_SHARED_LIB_EXPORT __declspec(dllexport)
  #define VF_SHARED_LIB_LOCAL
#else
  #if __GNUC__ >= 4
    #define VF_SHARED_LIB_IMPORT __attribute__ ((visibility ("default")))
    #define VF_SHARED_LIB_EXPORT __attribute__ ((visibility ("default")))
    #define VF_SHARED_LIB_LOCAL  __attribute__ ((visibility ("hidden")))
  #else
    #define VF_SHARED_LIB_IMPORT
    #define VF_SHARED_LIB_EXPORT
    #define VF_SHARED_LIB_LOCAL
  #endif
#endif

// Now we use the generic helper definitions above to define VF_PUBLIC, VF_PROTECTED
// and VF_PRIVATE. VF_PUBLIC is for symbols part of the public application programming
// interface (API), VF_PROTECTED is for symbols used e.g. by public templated or
// inlined code. These symbols must also be publicly available when compiling the
// application. VF_PRIVATE are symbols for internal use inside the library only.

#ifdef BUILD_SHARED_LIBS
   // defined if VF is compiled as a shared library
   #ifdef VF_SHARED_LIB_SELECT_IMPORTS
      // defined if we are building the VF SHARED_LIB (instead of using it)
      #define VF_PUBLIC VF_SHARED_LIB_IMPORT
   #else
      #define VF_PUBLIC VF_SHARED_LIB_EXPORT
   #endif
   #define VF_PRIVATE VF_SHARED_LIB_LOCAL
#else
   // VF_SHARED_LIB is not defined: this means VF is a static library
   #define VF_PUBLIC
   #define VF_PRIVATE
#endif
#define VF_PROTECTED VF_PUBLIC

#endif

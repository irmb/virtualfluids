
#ifndef VirtualFluids_GPU_EXPORT_H
#define VirtualFluids_GPU_EXPORT_H

#ifdef VirtualFluids_GPU_BUILT_AS_STATIC
#  define VirtualFluids_GPU_EXPORT
#  define VIRTUALFLUIDS_GPU_NO_EXPORT
#else
#  ifndef VirtualFluids_GPU_EXPORT
#    ifdef VirtualFluids_GPU_EXPORTS
        /* We are building this library */
#      define VirtualFluids_GPU_EXPORT __declspec(dllexport)
#    else
        /* We are using this library */
#      define VirtualFluids_GPU_EXPORT __declspec(dllimport)
#    endif
#  endif

#  ifndef VIRTUALFLUIDS_GPU_NO_EXPORT
#    define VIRTUALFLUIDS_GPU_NO_EXPORT 
#  endif
#endif

#ifndef VIRTUALFLUIDS_GPU_DEPRECATED
#  define VIRTUALFLUIDS_GPU_DEPRECATED __declspec(deprecated)
#endif

#ifndef VIRTUALFLUIDS_GPU_DEPRECATED_EXPORT
#  define VIRTUALFLUIDS_GPU_DEPRECATED_EXPORT VirtualFluids_GPU_EXPORT VIRTUALFLUIDS_GPU_DEPRECATED
#endif

#ifndef VIRTUALFLUIDS_GPU_DEPRECATED_NO_EXPORT
#  define VIRTUALFLUIDS_GPU_DEPRECATED_NO_EXPORT VIRTUALFLUIDS_GPU_NO_EXPORT VIRTUALFLUIDS_GPU_DEPRECATED
#endif

#if 0 /* DEFINE_NO_DEPRECATED */
#  ifndef VIRTUALFLUIDS_GPU_NO_DEPRECATED
#    define VIRTUALFLUIDS_GPU_NO_DEPRECATED
#  endif
#endif

#endif

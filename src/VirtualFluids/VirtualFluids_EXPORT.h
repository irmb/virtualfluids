
#ifndef VirtualFluids_EXPORT_H
#define VirtualFluids_EXPORT_H

#ifdef VirtualFluids_BUILT_AS_STATIC
#  define VirtualFluids_EXPORT
#  define VIRTUALFLUIDS_NO_EXPORT
#else
#  ifndef VirtualFluids_EXPORT
#    ifdef VirtualFluids_EXPORTS
        /* We are building this library */
#      define VirtualFluids_EXPORT __declspec(dllexport)
#    else
        /* We are using this library */
#      define VirtualFluids_EXPORT __declspec(dllimport)
#    endif
#  endif

#  ifndef VIRTUALFLUIDS_NO_EXPORT
#    define VIRTUALFLUIDS_NO_EXPORT 
#  endif
#endif

#ifndef VIRTUALFLUIDS_DEPRECATED
#  define VIRTUALFLUIDS_DEPRECATED __declspec(deprecated)
#endif

#ifndef VIRTUALFLUIDS_DEPRECATED_EXPORT
#  define VIRTUALFLUIDS_DEPRECATED_EXPORT VirtualFluids_EXPORT VIRTUALFLUIDS_DEPRECATED
#endif

#ifndef VIRTUALFLUIDS_DEPRECATED_NO_EXPORT
#  define VIRTUALFLUIDS_DEPRECATED_NO_EXPORT VIRTUALFLUIDS_NO_EXPORT VIRTUALFLUIDS_DEPRECATED
#endif

#if 0 /* DEFINE_NO_DEPRECATED */
#  ifndef VIRTUALFLUIDS_NO_DEPRECATED
#    define VIRTUALFLUIDS_NO_DEPRECATED
#  endif
#endif

#endif


#ifndef GridGenerator_EXPORT_H
#define GridGenerator_EXPORT_H

#ifdef GridGenerator_BUILT_AS_STATIC
#  define GridGenerator_EXPORT
#  define GRIDGENERATOR_NO_EXPORT
#else
#  ifndef GridGenerator_EXPORT
#    ifdef GridGenerator_EXPORTS
        /* We are building this library */
#      define GridGenerator_EXPORT __declspec(dllexport)
#    else
        /* We are using this library */
#      define GridGenerator_EXPORT __declspec(dllimport)
#    endif
#  endif

#  ifndef GRIDGENERATOR_NO_EXPORT
#    define GRIDGENERATOR_NO_EXPORT 
#  endif
#endif

#ifndef GRIDGENERATOR_DEPRECATED
#  define GRIDGENERATOR_DEPRECATED __declspec(deprecated)
#endif

#ifndef GRIDGENERATOR_DEPRECATED_EXPORT
#  define GRIDGENERATOR_DEPRECATED_EXPORT GridGenerator_EXPORT GRIDGENERATOR_DEPRECATED
#endif

#ifndef GRIDGENERATOR_DEPRECATED_NO_EXPORT
#  define GRIDGENERATOR_DEPRECATED_NO_EXPORT GRIDGENERATOR_NO_EXPORT GRIDGENERATOR_DEPRECATED
#endif

#if 0 /* DEFINE_NO_DEPRECATED */
#  ifndef GRIDGENERATOR_NO_DEPRECATED
#    define GRIDGENERATOR_NO_DEPRECATED
#  endif
#endif

#endif

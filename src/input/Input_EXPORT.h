
#ifndef Input_EXPORT_H
#define Input_EXPORT_H

#ifdef Input_BUILT_AS_STATIC
#  define Input_EXPORT
#  define INPUT_NO_EXPORT
#else
#  ifndef Input_EXPORT
#    ifdef Input_EXPORTS
        /* We are building this library */
#      define Input_EXPORT __declspec(dllexport)
#    else
        /* We are using this library */
#      define Input_EXPORT __declspec(dllimport)
#    endif
#  endif

#  ifndef INPUT_NO_EXPORT
#    define INPUT_NO_EXPORT 
#  endif
#endif

#ifndef INPUT_DEPRECATED
#  define INPUT_DEPRECATED __declspec(deprecated)
#endif

#ifndef INPUT_DEPRECATED_EXPORT
#  define INPUT_DEPRECATED_EXPORT Input_EXPORT INPUT_DEPRECATED
#endif

#ifndef INPUT_DEPRECATED_NO_EXPORT
#  define INPUT_DEPRECATED_NO_EXPORT INPUT_NO_EXPORT INPUT_DEPRECATED
#endif

#if 0 /* DEFINE_NO_DEPRECATED */
#  ifndef INPUT_NO_DEPRECATED
#    define INPUT_NO_DEPRECATED
#  endif
#endif

#endif

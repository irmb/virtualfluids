
#ifndef MuParser_EXPORT_H
#define MuParser_EXPORT_H

#ifdef MuParser_BUILT_AS_STATIC
#  define MuParser_EXPORT
#  define MUPARSER_NO_EXPORT
#else
#  ifndef MuParser_EXPORT
#    ifdef MuParser_EXPORTS
        /* We are building this library */
#      define MuParser_EXPORT __declspec(dllexport)
#    else
        /* We are using this library */
#      define MuParser_EXPORT __declspec(dllimport)
#    endif
#  endif

#  ifndef MUPARSER_NO_EXPORT
#    define MUPARSER_NO_EXPORT 
#  endif
#endif

#ifndef MUPARSER_DEPRECATED
#  define MUPARSER_DEPRECATED __declspec(deprecated)
#endif

#ifndef MUPARSER_DEPRECATED_EXPORT
#  define MUPARSER_DEPRECATED_EXPORT MuParser_EXPORT MUPARSER_DEPRECATED
#endif

#ifndef MUPARSER_DEPRECATED_NO_EXPORT
#  define MUPARSER_DEPRECATED_NO_EXPORT MUPARSER_NO_EXPORT MUPARSER_DEPRECATED
#endif

#if 0 /* DEFINE_NO_DEPRECATED */
#  ifndef MUPARSER_NO_DEPRECATED
#    define MUPARSER_NO_DEPRECATED
#  endif
#endif

#endif


#ifndef MUPARSER_EXPORT_H
#define MUPARSER_EXPORT_H

#ifdef MUPARSER_STATIC_DEFINE
#  define MUPARSER_EXPORT
#  define MUPARSER_NO_EXPORT
#else
#  ifndef MUPARSER_EXPORT
#    ifdef MuParser_EXPORTS
        /* We are building this library */
#      define MUPARSER_EXPORT 
#    else
        /* We are using this library */
#      define MUPARSER_EXPORT 
#    endif
#  endif

#  ifndef MUPARSER_NO_EXPORT
#    define MUPARSER_NO_EXPORT 
#  endif
#endif

#ifndef MUPARSER_DEPRECATED
#  define MUPARSER_DEPRECATED __attribute__ ((__deprecated__))
#endif

#ifndef MUPARSER_DEPRECATED_EXPORT
#  define MUPARSER_DEPRECATED_EXPORT MUPARSER_EXPORT MUPARSER_DEPRECATED
#endif

#ifndef MUPARSER_DEPRECATED_NO_EXPORT
#  define MUPARSER_DEPRECATED_NO_EXPORT MUPARSER_NO_EXPORT MUPARSER_DEPRECATED
#endif

#if 0 /* DEFINE_NO_DEPRECATED */
#  ifndef MUPARSER_NO_DEPRECATED
#    define MUPARSER_NO_DEPRECATED
#  endif
#endif

#endif /* MUPARSER_EXPORT_H */

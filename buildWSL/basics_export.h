
#ifndef BASICS_EXPORT_H
#define BASICS_EXPORT_H

#ifdef BASICS_STATIC_DEFINE
#  define BASICS_EXPORT
#  define BASICS_NO_EXPORT
#else
#  ifndef BASICS_EXPORT
#    ifdef basics_EXPORTS
        /* We are building this library */
#      define BASICS_EXPORT 
#    else
        /* We are using this library */
#      define BASICS_EXPORT 
#    endif
#  endif

#  ifndef BASICS_NO_EXPORT
#    define BASICS_NO_EXPORT 
#  endif
#endif

#ifndef BASICS_DEPRECATED
#  define BASICS_DEPRECATED __attribute__ ((__deprecated__))
#endif

#ifndef BASICS_DEPRECATED_EXPORT
#  define BASICS_DEPRECATED_EXPORT BASICS_EXPORT BASICS_DEPRECATED
#endif

#ifndef BASICS_DEPRECATED_NO_EXPORT
#  define BASICS_DEPRECATED_NO_EXPORT BASICS_NO_EXPORT BASICS_DEPRECATED
#endif

#if 0 /* DEFINE_NO_DEPRECATED */
#  ifndef BASICS_NO_DEPRECATED
#    define BASICS_NO_DEPRECATED
#  endif
#endif

#endif /* BASICS_EXPORT_H */

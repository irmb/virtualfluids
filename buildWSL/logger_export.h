
#ifndef LOGGER_EXPORT_H
#define LOGGER_EXPORT_H

#ifdef LOGGER_STATIC_DEFINE
#  define LOGGER_EXPORT
#  define LOGGER_NO_EXPORT
#else
#  ifndef LOGGER_EXPORT
#    ifdef logger_EXPORTS
        /* We are building this library */
#      define LOGGER_EXPORT 
#    else
        /* We are using this library */
#      define LOGGER_EXPORT 
#    endif
#  endif

#  ifndef LOGGER_NO_EXPORT
#    define LOGGER_NO_EXPORT 
#  endif
#endif

#ifndef LOGGER_DEPRECATED
#  define LOGGER_DEPRECATED __attribute__ ((__deprecated__))
#endif

#ifndef LOGGER_DEPRECATED_EXPORT
#  define LOGGER_DEPRECATED_EXPORT LOGGER_EXPORT LOGGER_DEPRECATED
#endif

#ifndef LOGGER_DEPRECATED_NO_EXPORT
#  define LOGGER_DEPRECATED_NO_EXPORT LOGGER_NO_EXPORT LOGGER_DEPRECATED
#endif

#if 0 /* DEFINE_NO_DEPRECATED */
#  ifndef LOGGER_NO_DEPRECATED
#    define LOGGER_NO_DEPRECATED
#  endif
#endif

#endif /* LOGGER_EXPORT_H */


#ifndef LBM_EXPORT_H
#define LBM_EXPORT_H

#ifdef LBM_STATIC_DEFINE
#  define LBM_EXPORT
#  define LBM_NO_EXPORT
#else
#  ifndef LBM_EXPORT
#    ifdef lbm_EXPORTS
        /* We are building this library */
#      define LBM_EXPORT 
#    else
        /* We are using this library */
#      define LBM_EXPORT 
#    endif
#  endif

#  ifndef LBM_NO_EXPORT
#    define LBM_NO_EXPORT 
#  endif
#endif

#ifndef LBM_DEPRECATED
#  define LBM_DEPRECATED __attribute__ ((__deprecated__))
#endif

#ifndef LBM_DEPRECATED_EXPORT
#  define LBM_DEPRECATED_EXPORT LBM_EXPORT LBM_DEPRECATED
#endif

#ifndef LBM_DEPRECATED_NO_EXPORT
#  define LBM_DEPRECATED_NO_EXPORT LBM_NO_EXPORT LBM_DEPRECATED
#endif

#if 0 /* DEFINE_NO_DEPRECATED */
#  ifndef LBM_NO_DEPRECATED
#    define LBM_NO_DEPRECATED
#  endif
#endif

#endif /* LBM_EXPORT_H */


#ifndef MPI_EXPORT_H
#define MPI_EXPORT_H

#ifdef MPI_STATIC_DEFINE
#  define MPI_EXPORT
#  define MPI_NO_EXPORT
#else
#  ifndef MPI_EXPORT
#    ifdef mpi_EXPORTS
        /* We are building this library */
#      define MPI_EXPORT 
#    else
        /* We are using this library */
#      define MPI_EXPORT 
#    endif
#  endif

#  ifndef MPI_NO_EXPORT
#    define MPI_NO_EXPORT 
#  endif
#endif

#ifndef MPI_DEPRECATED
#  define MPI_DEPRECATED __attribute__ ((__deprecated__))
#endif

#ifndef MPI_DEPRECATED_EXPORT
#  define MPI_DEPRECATED_EXPORT MPI_EXPORT MPI_DEPRECATED
#endif

#ifndef MPI_DEPRECATED_NO_EXPORT
#  define MPI_DEPRECATED_NO_EXPORT MPI_NO_EXPORT MPI_DEPRECATED
#endif

#if 0 /* DEFINE_NO_DEPRECATED */
#  ifndef MPI_NO_DEPRECATED
#    define MPI_NO_DEPRECATED
#  endif
#endif

#endif /* MPI_EXPORT_H */

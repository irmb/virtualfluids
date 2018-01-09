
#ifndef dem_coupling_EXPORT_H
#define dem_coupling_EXPORT_H

#ifdef dem_coupling_BUILT_AS_STATIC
#  define dem_coupling_EXPORT
#  define DEM_COUPLING_NO_EXPORT
#else
#  ifndef dem_coupling_EXPORT
#    ifdef dem_coupling_EXPORTS
        /* We are building this library */
#      define dem_coupling_EXPORT __declspec(dllexport)
#    else
        /* We are using this library */
#      define dem_coupling_EXPORT __declspec(dllimport)
#    endif
#  endif

#  ifndef DEM_COUPLING_NO_EXPORT
#    define DEM_COUPLING_NO_EXPORT 
#  endif
#endif

#ifndef DEM_COUPLING_DEPRECATED
#  define DEM_COUPLING_DEPRECATED __declspec(deprecated)
#endif

#ifndef DEM_COUPLING_DEPRECATED_EXPORT
#  define DEM_COUPLING_DEPRECATED_EXPORT dem_coupling_EXPORT DEM_COUPLING_DEPRECATED
#endif

#ifndef DEM_COUPLING_DEPRECATED_NO_EXPORT
#  define DEM_COUPLING_DEPRECATED_NO_EXPORT DEM_COUPLING_NO_EXPORT DEM_COUPLING_DEPRECATED
#endif

#if 0 /* DEFINE_NO_DEPRECATED */
#  ifndef DEM_COUPLING_NO_DEPRECATED
#    define DEM_COUPLING_NO_DEPRECATED
#  endif
#endif

#endif

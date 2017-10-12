
#ifndef VirtualFluidsBasics_EXPORT_H
#define VirtualFluidsBasics_EXPORT_H

#ifdef VirtualFluidsBasics_BUILT_AS_STATIC
#  define VirtualFluidsBasics_EXPORT
#  define VIRTUALFLUIDSBASICS_NO_EXPORT
#else
#  ifndef VirtualFluidsBasics_EXPORT
#    ifdef VirtualFluidsBasics_EXPORTS
        /* We are building this library */
#      define VirtualFluidsBasics_EXPORT __declspec(dllexport)
#    else
        /* We are using this library */
#      define VirtualFluidsBasics_EXPORT __declspec(dllimport)
#    endif
#  endif

#  ifndef VIRTUALFLUIDSBASICS_NO_EXPORT
#    define VIRTUALFLUIDSBASICS_NO_EXPORT 
#  endif
#endif

#ifndef VIRTUALFLUIDSBASICS_DEPRECATED
#  define VIRTUALFLUIDSBASICS_DEPRECATED __declspec(deprecated)
#endif

#ifndef VIRTUALFLUIDSBASICS_DEPRECATED_EXPORT
#  define VIRTUALFLUIDSBASICS_DEPRECATED_EXPORT VirtualFluidsBasics_EXPORT VIRTUALFLUIDSBASICS_DEPRECATED
#endif

#ifndef VIRTUALFLUIDSBASICS_DEPRECATED_NO_EXPORT
#  define VIRTUALFLUIDSBASICS_DEPRECATED_NO_EXPORT VIRTUALFLUIDSBASICS_NO_EXPORT VIRTUALFLUIDSBASICS_DEPRECATED
#endif

#if 0 /* DEFINE_NO_DEPRECATED */
#  ifndef VIRTUALFLUIDSBASICS_NO_DEPRECATED
#    define VIRTUALFLUIDSBASICS_NO_DEPRECATED
#  endif
#endif

#endif

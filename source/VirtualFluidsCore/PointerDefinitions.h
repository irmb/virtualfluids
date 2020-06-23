/*
*  Author: S. Peters
*  mail: peters@irmb.tu-bs.de
*/
#ifndef SHARED_POINTER_H
#define SHARED_POINTER_H

#define useStdSmartPointer

#ifdef useStdSmartPointer
  #include <memory>
  #define smartPointerNamespace std
#endif


template <class T>
using SPtr = smartPointerNamespace::shared_ptr<T>;

template <class T>
using WPtr = smartPointerNamespace::weak_ptr<T>;

//template <class T>
//using UPtr = smartPointerNamespace::unique_ptr<T>;

template <class T>
using enableSharedFromThis = smartPointerNamespace::enable_shared_from_this<T>;

#define dynamicPointerCast smartPointerNamespace::dynamic_pointer_cast

template <class T>
using RPtr = T*;

#endif

//#ifndef VF_BOOST
//   #include <memory>
//   #define smartPointerNamespace std
//#else
//   #include <boost/enable_shared_from_this.hpp>
//   #include <boost/pointer_cast.hpp>
//   #include <boost/shared_ptr.hpp>
//   #define smartPointerNamespace boost
//#endif
//
//#define SPtr smartPointerNamespace::shared_ptr 
//
//#define WPtr smartPointerNamespace::weak_ptr
//
//#define enableSharedFromThis smartPointerNamespace::enable_shared_from_this
//
//#define dynamicPointerCast smartPointerNamespace::dynamic_pointer_cast

//#endif
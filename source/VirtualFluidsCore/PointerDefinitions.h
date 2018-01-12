/*
*  Author: S. Peters
*  mail: peters@irmb.tu-bs.de
*/
#ifndef SHARED_POINTER_H
#define SHARED_POINTER_H

#include <memory>

#define useStdSmartPointer

#ifdef useStdSmartPointer
  #define smartPointerNamespace std
#else
  #define smartPointerNamespace boost
#endif


template <class T>
using SPtr = smartPointerNamespace::shared_ptr<T>;

template <class T>
using WPtr = smartPointerNamespace::weak_ptr<T>;

template <class T>
using UPtr = smartPointerNamespace::unique_ptr<T>;

template <class T>
using enableSharedFromThis = smartPointerNamespace::enable_shared_from_this<T>;

#define dynamicPointerCast smartPointerNamespace::dynamic_pointer_cast

template <class T>
using RPtr = T*;

#endif

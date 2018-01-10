/*
*  Author: S. Peters
*  mail: peters@irmb.tu-bs.de
*/
#ifndef SHARED_POINTER_H
#define SHARED_POINTER_H

#include <memory>

template <class T>
using SPtr = std::shared_ptr<T>;

template <class T>
using WPtr = std::weak_ptr<T>;

template <class T>
using UPtr = std::unique_ptr<T>;

template <class T>
using RPtr = T*;

#endif

//=======================================================================================
// ____          ____    __    ______     __________   __      __       __        __
// \    \       |    |  |  |  |   _   \  |___    ___| |  |    |  |     /  \      |  |
//  \    \      |    |  |  |  |  |_)   |     |  |     |  |    |  |    /    \     |  |
//   \    \     |    |  |  |  |   _   /      |  |     |  |    |  |   /  /\  \    |  |
//    \    \    |    |  |  |  |  | \  \      |  |     |   \__/   |  /  ____  \   |  |____
//     \    \   |    |  |__|  |__|  \__\     |__|      \________/  /__/    \__\  |_______|
//      \    \  |    |   ________________________________________________________________
//       \    \ |    |  |  ______________________________________________________________|
//        \    \|    |  |  |         __          __     __     __     ______      _______
//         \         |  |  |_____   |  |        |  |   |  |   |  |   |   _  \    /  _____)
//          \        |  |   _____|  |  |        |  |   |  |   |  |   |  | \  \   \_______
//           \       |  |  |        |  |_____   |   \_/   |   |  |   |  |_/  /    _____  |
//            \ _____|  |__|        |________|   \_______/    |__|   |______/    (_______/
//
//  This file is part of VirtualFluids. VirtualFluids is free software: you can
//  redistribute it and/or modify it under the terms of the GNU General Public
//  License as published by the Free Software Foundation, either version 3 of
//  the License, or (at your option) any later version.
//
//  VirtualFluids is distributed in the hope that it will be useful, but WITHOUT
//  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//  FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
//  for more details.
//
//  SPDX-License-Identifier: GPL-3.0-or-later
//  SPDX-FileCopyrightText: Copyright © VirtualFluids Project contributors, see AUTHORS.md in root folder
//
//! \author Henry Korb
//=======================================================================================
#include <basics/DataTypes.h>
#include <pybind11/detail/using_smart_holder.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <stdexcept>
#include <string>

namespace numpy_utils
{
namespace py = pybind11;

template <typename T>
void checkDim(const py::array_t<T>& arr, py::ssize_t ndim)
{
    if (arr.ndim() == ndim)
        return;
    throw std::runtime_error("array needs to have ndim " + std::to_string(ndim));
}
template <typename T>
void checkSize(const py::array_t<T>& arr, py::ssize_t dim, py::ssize_t size)
{
    if (size == 0 || arr.shape(dim) == size)
        return;
    throw std::runtime_error("array needs to have size " + std::to_string(size) + " in dimension" + std::to_string(dim));
}
template <size_t ndim, typename T>
void checkShape(const py::array_t<T>& arr, std::array<py::ssize_t, ndim> shape)
{
    checkDim(arr, ndim);
    for (size_t dim = 0; dim < ndim; dim++)
        checkSize(arr, dim, shape[dim]);
}

template <size_t N, typename T>
inline std::array<T, N> convertArrayToArray(const py::array_t<T>& arr)
{
    checkShape<1>(arr, { N });
    std::array<T, N> cpparr;
    for (size_t i; i < N; i++)
        cpparr[i] = arr.at(i);
    return cpparr;
}
inline real3 convertArrayToReal3(const py::array_t<real>& arr)
{
    checkShape<1>(arr, { 3 });
    return { arr.at(0), arr.at(1), arr.at(2) };
}
inline std::vector<real3> convertArrayToVecReal3(const py::array_t<real>& arr)
{
    checkShape<2>(arr, { 0, 3 });
    std::vector<real3> vec;

    for (py::ssize_t i = 0; i < arr.shape(0); i++)
        vec.emplace_back(real3 { arr.at(i, 0), arr.at(i, 1), arr.at(i, 2) });
    return vec;
}
inline std::vector<real> convertArrayToVecReal(const py::array_t<real>& arr)
{
    std::vector<real> vec;
    checkShape<1>(arr, { 0 });
    for (py::ssize_t i = 0; i < arr.shape(0); i++) {
        vec.emplace_back(arr.at(i));
    }
    return vec;
}

template <typename T>
inline std::vector<T> convertArrayToVec(const py::array_t<real>& arr);

template <>
inline std::vector<real> convertArrayToVec(const py::array_t<real>& arr)
{
    return convertArrayToVecReal(arr);
}

template <>
inline std::vector<real3> convertArrayToVec(const py::array_t<real>& arr)
{
    return convertArrayToVecReal3(arr);
}
} // namespace numpy_utils
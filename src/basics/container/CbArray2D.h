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
//! \addtogroup container
//! \ingroup basics
//! \{
//! \author Soeren Freudiger, Sebastian Geller
//=======================================================================================
#ifndef CBARRAY2D_H
#define CBARRAY2D_H

#include <iomanip>

#include <algorithm>
#include <basics/utilities/UbEqual.h>
#include <basics/utilities/UbException.h>
#include <typeinfo>

//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
// IndexClasses

// IndexerX2X1:
//        4 5 6
// Array  1 2 3  -->  vector 1 2 3 4 5 6
// optimaler schleifendurchlauf
// for(alle X2)
//  for(alle X1)
class IndexerX2X1
{
public:
    using size_type = int;

public:
    inline std::size_t getIndex(const size_type &x1, const size_type &x2, const size_type &nx1,
                                const size_type & /*nx2*/) const
    {
        return nx1 * x2 + x1;
    }
    inline std::size_t getStartIndexOfSortedArray(const size_type & /*x1*/, const size_type &x2, const size_type &nx1,
                                                  const size_type & /*nx2*/) const
    {
        return nx1 * x2;
    }
};

// IndexerX1X2:
//        4 5 6
// Array  1 2 3  -->  vector 1 4 2 5 3 6
// optimaler schleifendurchlauf
// for(alle X1)
//  for(alle X2)
class IndexerX1X2
{
public:
    using size_type = int;

public:
    inline std::size_t getIndex(const size_type &x1, const size_type &x2, const size_type & /*nx1*/,
                                const size_type &nx2) const
    {
        return nx2 * x1 + x2;
    }
    inline std::size_t getStartIndexOfSortedArray(const size_type &x1, const size_type & /*x2*/,
                                                  const size_type & /*nx1*/, const size_type &nx2) const
    {
        return nx2 * x1;
    }
};

//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
// CbArray2D
//////////////////////////////////////////////////////////////////////////
//! \brief 2D Array
//! \details the data is stored in a vector
//!
//! Rangecheck active, if:
//!
//! -debug  : not defined "NO_CB_RANGECHECK"
//!
//! -release: not defined "NO_CB_RANGECHECK" && defined "CB_RANGECHECK"
//////////////////////////////////////////////////////////////////////////
template <typename T, typename IndexClass = IndexerX2X1>
class CbArray2D
{
public:
    using value_type      = T;
    using indexer_type    = IndexClass;
    using size_type       = typename IndexClass::size_type;
    using reference       = typename std::vector<value_type>::reference;
    using const_reference = typename std::vector<value_type>::const_reference;
    using pointer         = typename std::vector<value_type>::pointer;
    using const_pointer   = typename std::vector<value_type>::const_pointer;

private:
    template <typename value_type2, typename IndexClass2>
    friend class CbArray2D;

public:
    /*=======================================================================*/
    CbArray2D() { this->resize(0, 0); }
    /*=======================================================================*/
    CbArray2D(const size_type &nx2, const size_type &nx1) { this->resize(nx2, nx1); }
    /*=======================================================================*/
    CbArray2D(const size_type &nx2, const size_type &nx1, const value_type &val) { this->resize(nx2, nx1, val); }
    /*=======================================================================*/
    CbArray2D(const size_type &uniformDimensionSize /*nx1==nx2*/)
    {
        this->resize(uniformDimensionSize, uniformDimensionSize);
    }
    /*=======================================================================*/
    // ssbernimmt vector als daten vector! (erstellt KEINE kopie!!!, vec ist anschlieueend leer, da swap verwendet wird)
    CbArray2D(std::vector<value_type> &vec, const size_type &nx1, const size_type &nx2)
    {
        assert((nx1 * nx2) == vec.size());
        this->data.swap(vec);
        this->resize(nx1, nx2);
    }
    /*=======================================================================*/
    CbArray2D(const CbArray2D &src) : nx1(src.nx1), nx2(src.nx2), data(src.data) {}
    /*=======================================================================*/
    template <typename value_type2>
    CbArray2D(const CbArray2D<value_type2> &src) : nx1(src.nx1), nx2(src.nx2)
    {
        // Sourcedaten kopieren
        this->data.resize(src.data.size());
        for (std::size_t i = 0; i < data.size(); ++i)
            this->data[i] = src.data[i];
    }
    /*=======================================================================*/
    virtual ~CbArray2D() = default;
    /*=======================================================================*/
    CbArray2D &operator=(const CbArray2D &rhs)
    {
        if (this == &rhs)
            return *this;

        this->nx1 = rhs.nx1;
        this->nx2 = rhs.nx2;

        // Laenge anpassen
        this->data.resize(rhs.data.size());
        // gespeicherte Datenelemente loeschen
        this->data.clear();

        // Sourcedaten kopieren
        this->data = rhs.data;

        return *this;
    }
    /*=======================================================================*/
    // durch value_type2 kann man z.B. ein float array einem double array zuweisen!
    template <typename value_type2, typename IndexClass2>
    CbArray2D &operator=(const CbArray2D<value_type2, IndexClass2> &rhs)
    {
        this->nx1 = rhs.nx1;
        this->nx2 = rhs.nx2;

        // gespeicherte Datenelemente loeschen
        this->data.clear();
        // Laenge anpassen
        this->data.resize(rhs.data.size());

        // Sourcedaten kopieren (!! koennte anderen Indexer besitzen!!! -> operator() benutzen)
        // ACHTUNG: fuer diese Konvertierung muss bei Klassen der demenstrechende operator
        //         implementiert sein, e.g.: class value_type2 {public: inline operator value_type2() const { return
        //         value_type2(); }
        for (int x1 = 0; x1 < this->nx1; x1++)
            for (int x2 = 0; x2 < this->nx2; x2++)
                this->operator()(x1, x2) = static_cast<value_type>(rhs.operator()(x1, x2));

        return *this;
    }
    /*=======================================================================*/
    bool operator==(const CbArray2D &rhs) const
    {
        if (this == &rhs)
            return true;

        if (this->nx1 != rhs.nx1 || this->nx2 != rhs.nx2 || this->data.size() != rhs.data.size()) {
            return false;
        }

        return std::equal(this->data.begin(), this->data.end(), rhs.data.begin(), UbEqual<value_type, value_type>());
    }
    /*=======================================================================*/
    template <typename value_type2, typename IndexClass2>
    bool operator==(const CbArray2D<value_type2, IndexClass2> &rhs) const
    {
        if (this->data.size() != rhs.data.size())
            return false;

        // Sourcedaten einzeln checken (!! koennte anderen Indexer besitzen!!! -> operator() benutzen)
        for (int x1 = 0; x1 < this->nx1; x1++)
            for (int x2 = 0; x2 < this->nx2; x2++)
                if (!isUbEqual(this->operator()(x1, x2), rhs.operator()(x1, x2)))
                    return false;

        return true;
    }
    /*=======================================================================*/
    bool operator!=(const CbArray2D &rhs) const { return !(*this == rhs); }
    /*=======================================================================*/
    template <typename value_type2, typename IndexClass2>
    bool operator!=(const CbArray2D<value_type2, IndexClass2> &rhs) const
    {
        return !(*this == rhs);
    }
    /*=======================================================================*/
    reference operator()(const size_type &x1, const size_type &x2)
    {
#if !defined(NO_CB_RANGECHECK) && (defined(_DEBUG) || defined(CB_RANGECHECK))
        if (!this->indicesInRange(x1, x2))
            UB_THROW(UbException(UB_EXARGS, getExceptionErrorString(x1, x2)));
#endif

        return this->data[indexer.getIndex(x1, x2, nx1, nx2)];
    }
    /*=======================================================================*/
    const_reference operator()(const size_type &x1, const size_type &x2) const
    {
#if !defined(NO_CB_RANGECHECK) && (defined(_DEBUG) || defined(CB_RANGECHECK))
        if (!this->indicesInRange(x1, x2))
            UB_THROW(UbException(UB_EXARGS, getExceptionErrorString(x1, x2)));
#endif

        return this->data[indexer.getIndex(x1, x2, nx1, nx2)];
    }
    /*=======================================================================*/
    pointer getStartAdressOfSortedArray(const size_type &x1, const size_type &x2)
    {
#if !defined(NO_CB_RANGECHECK) && (defined(_DEBUG) || defined(CB_RANGECHECK))
        if (!this->indicesInRange(x1, x2))
            UB_THROW(UbException(UB_EXARGS, getExceptionErrorString(x1, x2)));
#endif
        return &this->data[indexer.getStartIndexOfSortedArray(x1, x2, nx1, nx2)];
    }
    /*=======================================================================*/
    const_pointer getStartAdressOfSortedArray(const size_type &x1, const size_type &x2) const
    {
#if !defined(NO_CB_RANGECHECK) && (defined(_DEBUG) || defined(CB_RANGECHECK))
        if (!this->indicesInRange(x1, x2))
            UB_THROW(UbException(UB_EXARGS, getExceptionErrorString(x1, x2)));
#endif
        return &this->data[indexer.getStartIndexOfSortedArray(x1, x2, nx1, nx2)];
    }
    /*=======================================================================*/
    void setObject(const size_type &x1, const size_type &x2, const value_type &value)
    {
#if !defined(NO_CB_RANGECHECK) && (defined(_DEBUG) || defined(CB_RANGECHECK))
        if (!this->indicesInRange(x1, x2))
            UB_THROW(UbException(UB_EXARGS, getExceptionErrorString(x1, x2)));
#endif
        this->data[indexer.getIndex(x1, x2, nx1, nx2)] = value;
    }
    /*=======================================================================*/
    reference getObject(const size_type &x1, const size_type &x2)
    {
#if !defined(NO_CB_RANGECHECK) && (defined(_DEBUG) || defined(CB_RANGECHECK))
        if (!this->indicesInRange(x1, x2))
            UB_THROW(UbException(UB_EXARGS, getExceptionErrorString(x1, x2)));
#endif
        return this->data[indexer.getIndex(x1, x2, nx1, nx2)];
    }
    /*=======================================================================*/
    typename std::vector<value_type>::const_reference getObject(const size_type &x1, const size_type &x2) const
    {
        return this->operator()(x1, x2);
    }
    /*=======================================================================*/
    bool isEmpty() const { return data.empty(); }
    size_type getNX1() const { return this->nx1; }
    size_type getNX2() const { return this->nx2; }
    /*=======================================================================*/
    void reset(const T &val) { std::fill(this->data.begin(), this->data.end(), val); }
    /*=======================================================================*/
    std::string toString() const
    {
        std::stringstream text;
        for (size_type x2 = 0; x2 < this->nx2; x2++) {
            for (size_type x1 = 0; x1 < this->nx1; x1++) {
                // hier kommts zum Konflikt ab  und an ...
                text << this->getObject(x1, x2) << ", ";
            }
            text << "\n";
        }

        return text.str();
    }
    /*=======================================================================*/
    std::string getInfo() const
    {
        std::stringstream text;
        text << "CbArray2D< storageType=" << typeid(T).name() << ", indexer=" << typeid(IndexClass).name() << " >";
        text << "( nx1=" << this->nx1 << ", nx2=" << this->nx2 << ")";
        return text.str();
    }
    /*=======================================================================*/
    void resize(const size_type &uniformDimensionSize) { this->resize(uniformDimensionSize, uniformDimensionSize); }
    /*=======================================================================*/
    void resize(const size_type &nx1, const size_type &nx2)
    {
        this->nx1 = nx1;
        this->nx2 = nx2;
        this->data.resize(nx1 * nx2);
    }
    /*=======================================================================*/
    void resize(const size_type &nx1, const size_type &nx2, const value_type &initVal)
    {
        this->nx1 = nx1;
        this->nx2 = nx2;
        this->data.resize(nx1 * nx2, initVal);
    }
    /*=======================================================================*/
    void clear()
    {
        this->nx1 = 0;
        this->nx2 = 0;
        this->data.clear();
    }
    /*=======================================================================*/
    std::vector<value_type> &getDataVector() { return this->data; }
    /*=======================================================================*/
    const std::vector<value_type> &getDataVector() const { return this->data; }
    /*=======================================================================*/
    inline size_type getDataVectorIndex(const size_type &x1, const size_type &x2) const
    {
#if !defined(NO_CB_RANGECHECK) && (defined(_DEBUG) || defined(CB_RANGECHECK))
        if (!this->indicesInRange(x1, x2))
            UB_THROW(UbException(UB_EXARGS, getExceptionErrorString(x1, x2)));
#endif

        return indexer.getIndex(x1, x2, nx1, nx2);
    }

protected:
    /*=======================================================================*/
    // success -> true
    // else    -> false
    inline bool indicesInRange(const size_type &x1, const size_type &x2) const
    {
        if (x1 < 0 || x1 >= this->nx1 || x2 < 0 || x2 >= this->nx2) {
            return false;
        }
        return true;
    }
    /*=======================================================================*/
    std::string getExceptionErrorString(const size_type &x1, const size_type &x2) const
    {
        std::stringstream out("index out of range - ");
        out << "(" << x1 << "," << x2 << ") not in (" << nx1 << "," << nx2 << ")";
        return out.str();
    }
    /*=======================================================================*/

protected:
    size_type nx1;
    size_type nx2;
    indexer_type indexer;
    std::vector<value_type> data;
};

#endif // CBARRAY2D_H

//! \}

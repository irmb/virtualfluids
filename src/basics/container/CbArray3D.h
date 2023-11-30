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
//  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
//  for more details.
//
//  You should have received a copy of the GNU General Public License along
//  with VirtualFluids (see COPYING.txt). If not, see <http://www.gnu.org/licenses/>.
//
//! \file CbArray3D.h
//! \ingroup container
//! \author Soeren Freudiger, Sebastian Geller
//=======================================================================================
#ifndef CBARRAY3D_H
#define CBARRAY3D_H

#include <iomanip>

#include <basics/PointerDefinitions.h>
#include <algorithm>
#include <basics/utilities/UbEqual.h>
#include <basics/utilities/UbException.h>
#include <typeinfo>

//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
// IndexClasses

// IndexerX3X2X1:
//                4 5 6          10 11 12
// Array  ebene A 1 2 3  ebene B  7  8  9 -->  vector 1 2 3 4 5 6 7 8 9 10 11 12
// x1-reihen "liegen am stueck" im speicher
// optimaler schleifendurchlauf
// for(alle X3)
//  for(alle X2)
//    for(alle X1)
class IndexerX3X2X1 // FunctorX1SortedForX1X2Plane
{
public:
    using size_type = size_t;

public:
    inline std::size_t getIndex(const size_type &x1, const size_type &x2, const size_type &x3, const size_type &nx1,
                                const size_type &nx2, const size_type & /*nx3*/) const
    {
        return nx1 * (nx2 * x3 + x2) + x1;
    }
    inline std::size_t getStartIndexOfSortedArray(const size_type & /*x1*/, const size_type &x2, const size_type &x3,
                                                  const size_type &nx1, const size_type &nx2,
                                                  const size_type & /*nx3*/) const
    {
        return nx1 * (nx2 * x3 + x2);
    }
};

// IndexerX1X2X3:
//                4 5 6          10 11 12
// Array  ebene A 1 2 3  ebene B  7  8  9 -->
// optimaler schleifendurchlauf
// for(alle X1)
//  for(alle X2)
//    for(alle X3)
class IndexerX1X2X3 // FunctorX3SortedForX3X2Plane
{
public:
    using size_type = size_t;

public:
    inline std::size_t getIndex(const size_type &x1, const size_type &x2, const size_type &x3,
                                const size_type & /*nx1*/, const size_type &nx2, const size_type &nx3) const
    {
        return nx3 * (nx2 * x1 + x2) + x3;
    }
    inline std::size_t getStartIndexOfSortedArray(const size_type &x1, const size_type &x2, const size_type & /*x3*/
                                                  ,
                                                  const size_type & /*nx1*/, const size_type &nx2,
                                                  const size_type &nx3) const
    {
        return nx3 * (nx2 * x1 + x2);
    }
};

// IndexerX2X1X3:
//                4 5 6          10 11 12
// Array  ebene A 1 2 3  ebene B  7  8  9 -->  vector 1 7 2 8 3 9 4 10 5 11 6 12
// optimaler schleifendurchlauf
// for(alle X2)
//  for(alle X1)
//    for(alle X3)
class IndexerX2X1X3
{
public:
    using size_type = size_t;

public:
    inline std::size_t getIndex(const size_type &x1, const size_type &x2, const size_type &x3, const size_type &nx1,
                                const size_type & /*nx2*/, const size_type &nx3) const
    {
        return nx3 * (nx1 * x2 + x1) + x3;
    }
    inline std::size_t getStartIndexOfSortedArray(const size_type &x1, const size_type &x2, const size_type & /*x3*/
                                                  ,
                                                  const size_type &nx1, const size_type & /*nx2*/,
                                                  const size_type &nx3) const
    {
        return nx3 * (nx1 * x2 + x1);
    }
};

//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
// CbArray3D
//////////////////////////////////////////////////////////////////////////
//! \brief 3D Array
//! \details The data is stored in a vector
//!
//! Rangecheck active, if:
//!
//! -debug  : not defined "NO_CB_RANGECHECK"
//!
//! -release: not defined "NO_CB_RANGECHECK" && defined "CB_RANGECHECK"
//////////////////////////////////////////////////////////////////////////
template <typename T, typename IndexClass = IndexerX3X2X1>
class CbArray3D
{
public:
    using CbArray3DPtr = SPtr<CbArray3D<T, IndexClass>>;

    using value_type      = T;
    using indexer_type    = IndexClass;
    using size_type       = typename IndexClass::size_type;
    using reference       = typename std::vector<value_type>::reference;
    using const_reference = typename std::vector<value_type>::const_reference;
    using pointer         = typename std::vector<value_type>::pointer;
    using const_pointer   = typename std::vector<value_type>::const_pointer;

private:
    template <typename value_type2, typename IndexClass2>
    friend class CbArray3D;

public:
    /*=======================================================================*/
    CbArray3D() { this->resize(0, 0, 0); }
    /*=======================================================================*/
    CbArray3D(const size_type &nx1, const size_type &nx2, const size_type &nx3, const value_type &val)
    {
        this->resize(nx1, nx2, nx3, val);
    }
    /*=======================================================================*/
    CbArray3D(const size_type &nx1, const size_type &nx2, const size_type &nx3) { this->resize(nx1, nx2, nx3); }
    /*=======================================================================*/
    CbArray3D(const size_type &uniformDimensionSize /*nx1==nx2==nx3*/)
    {
        this->resize(uniformDimensionSize, uniformDimensionSize, uniformDimensionSize);
    }
    /*=======================================================================*/
    // ssbernimmt vector als daten vector! (erstellt KEINE kopie!!!, vec ist anschliessend leer, da swap verwendet wird)
    CbArray3D(std::vector<value_type> &vec, const size_type &nx1, const size_type &nx2, const size_type &nx3)
    {
        assert((nx1 * nx2 * nx3) == vec.size());
        this->data.swap(vec);
        this->resize(nx1, nx2, nx3);
    }
    /*=======================================================================*/
    CbArray3D(const CbArray3D &src) : nx1(src.nx1), nx2(src.nx2), nx3(src.nx3), data(src.data) {}
    /*=======================================================================*/
    template <typename value_type2>
    CbArray3D(const CbArray3D<value_type2> &src) : nx1(src.nx1), nx2(src.nx2), nx3(src.nx3)
    {
        // Sourcedaten kopieren
        this->data.resize(src.data.size());
        for (std::size_t i = 0; i < data.size(); ++i)
            this->data[i] = src.data[i];
    }
    /*=======================================================================*/
    virtual ~CbArray3D() = default;
    /*=======================================================================*/
    CbArray3D &operator=(const CbArray3D &rhs)
    {
        if (this == &rhs)
            return *this;

        this->nx1 = rhs.nx1;
        this->nx2 = rhs.nx2;
        this->nx3 = rhs.nx3;

        // gespeicherte Datenelemente loeschen
        // Laenge anpassen
        this->data.resize(rhs.data.size());
        // gespeicherte Datenelemente loeschen
        this->data.clear();

        // Sourcedaten kopieren
        this->data = rhs.data;

        return *this;
    }
    /*=======================================================================*/
    // durch value_type2 kann man z.B. ein float array einer double array zuweisen!
    template <typename value_type2, typename IndexClass2>
    CbArray3D &operator=(const CbArray3D<value_type2, IndexClass2> &rhs)
    {
        this->nx1 = rhs.nx1;
        this->nx2 = rhs.nx2;
        this->nx3 = rhs.nx3;

        // gespeicherte Datenelemente loeschen
        this->data.clear();
        // Laenge anpassen
        this->data.resize(rhs.data.size());

        // Sourcedaten kopieren (!! koennte anderen Indexer besitzen!!! -> operator() benutzen)
        for (int x3 = 0; x3 < this->nx3; x3++)
            for (int x2 = 0; x2 < this->nx2; x2++)
                for (int x1 = 0; x1 < this->nx1; x1++)
                    this->operator()(x1, x2, x3) = static_cast<value_type>(rhs.operator()(x1, x2, x3));

        return *this;
    }
    /*=======================================================================*/
    bool operator==(const CbArray3D &rhs) const
    {
        if (this == &rhs)
            return true;

        if (this->nx1 != rhs.nx1 || this->nx2 != rhs.nx2 || this->nx3 != rhs.nx3 ||
            this->data.size() != rhs.data.size()) {
            return false;
        }

        return std::equal(this->data.begin(), this->data.end(), rhs.data.begin(), UbEqual<value_type, value_type>());
    }
    /*=======================================================================*/
    template <typename value_type2, typename IndexClass2>
    bool operator==(const CbArray3D<value_type2, IndexClass2> &rhs) const
    {
        if (this->data.size() != rhs.data.size())
            return false;

        // Sourcedaten einzeln checken (!! koennte anderen Indexer besitzen!!! -> operator() benutzen)
        for (int x3 = 0; x3 < this->nx3; x3++)
            for (int x2 = 0; x2 < this->nx2; x2++)
                for (int x1 = 0; x1 < this->nx1; x1++)
                    if (!isUbEqual(this->operator()(x1, x2, x3), rhs.operator()(x1, x2, x3)))
                        return false;

        return true;
    }
    /*=======================================================================*/
    bool operator!=(const CbArray3D &src) const { return !(*this == src); }
    /*=======================================================================*/
    template <typename value_type2, typename IndexClass2>
    bool operator!=(const CbArray3D<value_type2, IndexClass2> &rhs) const
    {
        return !(*this == rhs);
    }
    /*=======================================================================*/
    reference operator()(const size_type &x1, const size_type &x2, const size_type &x3)
    {
#if !defined(NO_CB_RANGECHECK) && (defined(_DEBUG) || defined(CB_RANGECHECK))
        if (!this->indicesInRange(x1, x2, x3))
            UB_THROW(UbException(UB_EXARGS, getExceptionErrorString(x1, x2, x3)));
#endif

        return this->data[indexer.getIndex(x1, x2, x3, nx1, nx2, nx3)];
    }
    /*=======================================================================*/
    const_reference operator()(const size_type &x1, const size_type &x2, const size_type &x3) const
    {
#if !defined(NO_CB_RANGECHECK) && (defined(_DEBUG) || defined(CB_RANGECHECK))
        if (!this->indicesInRange(x1, x2, x3))
            UB_THROW(UbException(UB_EXARGS, getExceptionErrorString(x1, x2, x3)));
#endif

        return this->data[indexer.getIndex(x1, x2, x3, nx1, nx2, nx3)];
    }
    /*=======================================================================*/
    pointer getStartAdressOfSortedArray(const size_type &x1, const size_type &x2, const size_type &x3)
    {
#if !defined(NO_CB_RANGECHECK) && (defined(_DEBUG) || defined(CB_RANGECHECK))
        if (!this->indicesInRange(x1, x2, x3))
            UB_THROW(UbException(UB_EXARGS, getExceptionErrorString(x1, x2, x3)));
#endif

        return &this->data[indexer.getStartIndexOfSortedArray(x1, x2, x3, nx1, nx2, nx3)];
    }
    /*=======================================================================*/
    const_pointer getStartAdressOfSortedArray(const size_type &x1, const size_type &x2, const size_type &x3) const
    {
#if !defined(NO_CB_RANGECHECK) && (defined(_DEBUG) || defined(CB_RANGECHECK))
        if (!this->indicesInRange(x1, x2, x3))
            UB_THROW(UbException(UB_EXARGS, getExceptionErrorString(x1, x2, x3)));
#endif

        return &this->data[indexer.getStartIndexOfSortedArray(x1, x2, x3, nx1, nx2, nx3)];
    }
    /*=======================================================================*/
    void setObject(const size_type &x1, const size_type &x2, const size_type &x3, const value_type &value)
    {
#if !defined(NO_CB_RANGECHECK) && (defined(_DEBUG) || defined(CB_RANGECHECK))
        if (!this->indicesInRange(x1, x2, x3))
            UB_THROW(UbException(UB_EXARGS, getExceptionErrorString(x1, x2, x3)));
#endif

        this->data[indexer.getIndex(x1, x2, x3, nx1, nx2, nx3)] = value;
    }
    /*=======================================================================*/
    reference getObject(const size_type &x1, const size_type &x2, const size_type &x3)
    {
#if !defined(NO_CB_RANGECHECK) && (defined(_DEBUG) || defined(CB_RANGECHECK))
        if (!this->indicesInRange(x1, x2, x3))
            UB_THROW(UbException(UB_EXARGS, getExceptionErrorString(x1, x2, x3)));
#endif

        return this->data[indexer.getIndex(x1, x2, x3, nx1, nx2, nx3)];
    }
    /*=======================================================================*/
    const_reference getObject(const size_type &x1, const size_type &x2, const size_type &x3) const
    {
        return (*this)(x1, x2, x3);
    }
    /*=======================================================================*/
    bool isEmpty() const { return data.empty(); }
    size_type getNX1() const { return this->nx1; }
    size_type getNX2() const { return this->nx2; }
    size_type getNX3() const { return this->nx3; }
    /*=======================================================================*/
    void reset(const value_type &val) { std::fill(this->data.begin(), this->data.end(), val); }
    /*=======================================================================*/
    std::string toString() const
    {
        std::stringstream text;
        for (size_type x1 = 0; x1 < this->nx1; x1++) {
            for (size_type x2 = 0; x2 < this->nx2; x2++) {
                for (size_type x3 = 0; x3 < this->nx3; x3++) {
                    text << (*this)(x1, x2, x3) << ", ";
                }
                text << std::endl;
            }
            text << std::endl << std::endl;
        }

        return text.str();
    }
    /*=======================================================================*/
    std::string getInfo() const
    {
        std::stringstream text;
        text << "CbArray3D< storageType=" << typeid(T).name() << ", indexer=" << typeid(IndexClass).name() << " >";
        text << "( nx1=" << this->nx1 << ", nx2=" << this->nx2 << ", nx3=" << this->nx3 << ")";
        return text.str();
    }
    /*=======================================================================*/
    void resize(const int &uniformDimensionSize)
    {
        this->resize(uniformDimensionSize, uniformDimensionSize, uniformDimensionSize);
    }
    /*=======================================================================*/
    void resize(const size_type &nx1, const size_type &nx2, const size_type &nx3)
    {
        this->nx1 = nx1;
        this->nx2 = nx2;
        this->nx3 = nx3;
        this->data.resize(nx1 * nx2 * nx3);
    }
    /*=======================================================================*/
    void resize(const size_type &nx1, const size_type &nx2, const size_type &nx3, const value_type &val)
    {
        this->nx1 = nx1;
        this->nx2 = nx2;
        this->nx3 = nx3;
        this->data.resize(nx1 * nx2 * nx3, val);
    }
    /*=======================================================================*/
    void clear()
    {
        this->nx1 = 0;
        this->nx2 = 0;
        this->nx3 = 0;
        this->data.clear();
    }
    /*=======================================================================*/
    std::vector<value_type> &getDataVector() { return this->data; }
    /*=======================================================================*/
    const std::vector<value_type> &getDataVector() const { return this->data; }
    /*=======================================================================*/
    inline std::size_t getDataVectorIndex(const size_type &x1, const size_type &x2, const size_type &x3) const
    {
#if !defined(NO_CB_RANGECHECK) && (defined(_DEBUG) || defined(CB_RANGECHECK))
        if (!this->indicesInRange(x1, x2, x3))
            UB_THROW(UbException(UB_EXARGS, getExceptionErrorString(x1, x2, x3)));
#endif

        return indexer.getIndex(x1, x2, x3, nx1, nx2, nx3);
    }

    /*=======================================================================*/
    // success -> true
    // else    -> false
    inline bool indicesInRange(const size_type &x1, const size_type &x2, const size_type &x3) const
    {
        if (x1 < 0 || x1 >= this->nx1 || x2 < 0 || x2 >= this->nx2 || x3 < 0 || x3 >= this->nx3) {
            return false;
        }
        return true;
    }

protected:
    /*=======================================================================*/
    std::string getExceptionErrorString(const size_type &x1, const size_type &x2, const size_type &x3) const
    {
        std::stringstream out("index out of range - ");
        out << "(" << x1 << "," << x2 << "," << x3 << ") not in (" << nx1 << "," << nx2 << "," << nx3 << ")";
        return out.str();
    }
    /*=======================================================================*/

protected:
    size_type nx1;
    size_type nx2;
    size_type nx3;
    indexer_type indexer;
    std::vector<value_type> data;
};

#endif // CBARRAY3D_H

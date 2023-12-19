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
#ifndef CBVECTOR_H
#define CBVECTOR_H

#include <PointerDefinitions.h> //for memcopy
#include <algorithm>            //for std::swap
#include <typeinfo>             //for typeid
#include <vector>

#include <basics/utilities/UbEqual.h>
#include <basics/utilities/UbSystem.h>

template <typename T>
class CbVectorAllocator;
template <typename T>
class CbVectorAllocatorStd;

//=========================================================================
//! \brief A class implements a container like a vector
//! \details
//! For this class it was required to ave only the type as template argument.
//! Hence, the allocator must be an abstract class. With out this requirement,
//! an allocator as second template argument would have been possible, as in the
//! STL vector. This would lead to the compiler generating two different classes
//! for the same data type with different allocators during compile time. Here it
//! is required that the same class can have different allocators.
//!
//! Rangecheck active, if:
//! -debug  : not defined "NO_CB_RANGECHECK"
//! -release: not defined "NO_CB_RANGECHECK" && defined "CB_RANGECHECK"
//=========================================================================
//////////////////////////////////////////////////////////////////////////
template <typename T>
class CbVector
{
public:
    using value_type = T;
    using pointer    = value_type *;
    using size_type  = std::size_t;

    friend class CbVectorAllocator<value_type>; // um auf ptrData und dataSize zugreifen zu koennen!

    CbVector<value_type>(const CbVector<value_type> &src) = delete;

public:
    /*==========================================================*/
    CbVector(CbVectorAllocator<value_type> *const &allocator = new CbVectorAllocatorStd<value_type>)
        : ptrData(NULL), allocator(allocator)
    {
        this->allocator->alloc(*this, 0, value_type());
    }
    /*==========================================================*/
    CbVector(const size_type size,
             CbVectorAllocator<value_type> *const &allocator = new CbVectorAllocatorStd<value_type>,
             const value_type &value                         = value_type())
        : ptrData(NULL), allocator(allocator)
    {
        this->allocator->alloc(*this, size, value);
    }
    /*==========================================================*/
    virtual ~CbVector()
    {
        if (allocator) {
            this->allocator->dealloc(*this);
            delete allocator;
            allocator = NULL;
        }
    }
    /*=======================================================================*/
    CbVector &operator=(const CbVector &src)
    {
        if (this == &src)
            return *this;

        // gespeicherte Datenelemente loeschen
        // Laenge anpassen
        this->allocator->resize(*this, src.size());

        // gespeicherte Datenelemente kopieren
        if (!src.empty()) {
            memcpy((char *)ptrData, (char *)&src[0], src.size() * sizeof(value_type));
            // for(size_type i=0; i<src.size(); i++)
            //   (*this)[i] = src[i];
        }

        return *this;
    }
    /*=======================================================================*/
    CbVector &operator=(const std::vector<value_type> &src)
    {
        // gespeicherte Datenelemente loeschen
        // Laenge anpassen
        this->allocator->resize(*this, src.size());

        // gespeicherte Datenelemente kopieren
        if (!src.empty()) {
            memcpy((char *)ptrData, (char *)&src[0], src.size() * sizeof(value_type));
            // for(size_type i=0; i<src.size(); i++)
            //   (*this)[i] = src[i];
        }

        return *this;
    }
    /*=======================================================================*/
    bool operator==(const CbVector &rhs) const
    {
        if (this == &rhs)
            return true;
        if (this->dataSize != rhs.dataSize)
            return false;

        for (size_type i = 0; i < rhs.size(); i++)
            if (!isUbEqual(this->operator[](i), rhs.operator[](i)))
                return false;

        return true;
    }
    /*==========================================================*/
    void setAllocator(CbVectorAllocator<value_type> *const &allocator)
    {
        if (this->allocator) {
            if (this->allocator == allocator)
                return;
            this->allocator->dealloc(*this);
            delete this->allocator;
        }
        this->allocator = allocator;
        this->allocator->alloc(*this, 0);
    }
    /*==========================================================*/
    size_type size() const { return dataSize; }
    /*==========================================================*/
    bool empty() const { return dataSize == 0; }
    /*==========================================================*/
    bool resize(const size_type &dataSize) { return allocator->resize(*this, dataSize); }
    /*==========================================================*/
    bool resize(const size_type &dataSize, const value_type &value)
    {
        return allocator->resize(*this, dataSize, value);
    }
    /*==========================================================*/
    void swap(CbVector &rhs)
    {
        if (this == &rhs)
            return;

        std::swap(this->ptrData, rhs.ptrData);
        std::swap(this->dataSize, rhs.dataSize);
        std::swap(this->allocator, rhs.allocator);
    }
    /*==========================================================*/
    value_type &operator[](const size_type &i)
    {
#if !defined(NO_CB_RANGECHECK) && (defined(_DEBUG) || defined(CB_RANGECHECK))
        if (i >= dataSize)
            UB_THROW(UbException(UB_EXARGS, "T=" + (std::string) typeid(*this).name() + UbSystem::toString(i) +
                                                " out of range (size=" + UbSystem::toString(dataSize) + ")"));
#endif // _DEBUG

        return ptrData[i];
    }
    /*==========================================================*/
    const value_type &operator[](const size_type &i) const
    {
#if !defined(NO_CB_RANGECHECK) && (defined(_DEBUG) || defined(CB_RANGECHECK))
        if (i >= dataSize)
            UB_THROW(UbException(UB_EXARGS, "T=" + (std::string) typeid(*this).name() + UbSystem::toString(i) +
                                                " out of range (size=" + UbSystem::toString(dataSize) + ")"));
#endif // _DEBUG

        return ptrData[i];
    }
    /*==========================================================*/
    CbVectorAllocator<value_type> *getAllocator() const { return allocator; }
    /*==========================================================*/

private:
    value_type *ptrData;
    size_type dataSize{ 0 };
    CbVectorAllocator<value_type> *allocator;
    // CbVector<value_type>& operator=(const CbVector<value_type>& src);
};

//////////////////////////////////////////////////////////////////////////
// CbVectorAllocator-Interface
//////////////////////////////////////////////////////////////////////////
template <typename T>
class CbVectorAllocator
{
public:
    using value_type = typename CbVector<T>::value_type;
    using size_type  = typename CbVector<value_type>::size_type;

public:
    CbVectorAllocator()          = default;
    virtual ~CbVectorAllocator() = default;

    virtual bool alloc(CbVector<value_type> &vec, const size_type &dataSize,
                       const value_type &value = value_type())  = 0;
    virtual bool resize(CbVector<value_type> &vec, const size_type &dataSize,
                        const value_type &value = value_type()) = 0;
    virtual bool dealloc(CbVector<value_type> &vec)             = 0;

protected:
    // folgende Methoden ersparen eine friend Deklaierung aller moeglichen Allocatoren
    // denn durch diese beiden Methoden haben sie exklusive Zugriffsrechte!
    //**********************************************************************************//
    inline value_type *&ptrDataOf(CbVector<value_type> &vec)
    {
        if (vec.getAllocator() != this)
            UB_THROW(UbException(UB_EXARGS, "allocator is not member of vec!"));
        return vec.ptrData;
    }
    //**********************************************************************************//
    inline size_type &dataSizeOf(CbVector<value_type> &vec)
    {
        if (vec.getAllocator() != this)
            UB_THROW(UbException(UB_EXARGS, "allocator is not member of vec!"));
        return vec.dataSize;
    }
};

//////////////////////////////////////////////////////////////////////////
// CbVectorAllocatorStd
//////////////////////////////////////////////////////////////////////////
template <typename T>
class CbVectorAllocatorStd : public CbVectorAllocator<T>
{
public:
    // typedefs wiederholen, da Basisklasse = template -> "Dependent-Base"-Problem
    using value_type = typename CbVector<T>::value_type;
    using size_type  = typename CbVector<value_type>::size_type;

public:
    CbVectorAllocatorStd() : CbVectorAllocator<value_type>() {}
    /*==========================================================*/
    bool alloc(CbVector<value_type> &src, const size_type &dataSize, const value_type &value = value_type()) override
    {
        return this->resize(src, dataSize, value);
    }
    /*==========================================================*/
    bool resize(CbVector<value_type> &vec, const size_type &dataSize, const value_type &value = value_type()) override
    {
        if (CbVectorAllocatorStd<value_type>::dataSizeOf(vec) == dataSize)
            return false;

        // new array
        value_type *new_data = new value_type[dataSize];
        // copy existing data to array
        if (this->ptrDataOf(vec)) {
            for (size_type i = 0; (i < vec.size() && i < dataSize); ++i)
                new_data[i] = CbVectorAllocatorStd<value_type>::ptrDataOf(vec)[i];
            delete[] this->ptrDataOf(vec);
        }
        this->ptrDataOf(vec) = new_data;
        // new value for new items
        for (size_type i = this->dataSizeOf(vec); i < dataSize; ++i)
            this->ptrDataOf(vec)[i] = value;
        // assign new dataSize
        this->dataSizeOf(vec) = dataSize;

        return true;
    }
    /*==========================================================*/
    bool dealloc(CbVector<value_type> &vec) override
    {
        if (this->ptrDataOf(vec)) {
            delete[] this->ptrDataOf(vec);
            this->ptrDataOf(vec) = NULL;
        }
        this->dataSizeOf(vec) = 0;
        return true;
    }
    /*==========================================================*/

private:
};

#endif // CBVECTOR_H

//! \}

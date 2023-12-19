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
//! \addtogroup memory
//! \ingroup basics
//! \{
//! \author Soeren Freudiger, Sebastian Geller
//=======================================================================================
#ifndef MBSMARTPTR_H
#define MBSMARTPTR_H

#include <basics/memory/MbSmartPtrBase.h>

#ifdef __clang__
#pragma clang system_header
#endif

//=====================================================
// Globale Funktion, um das Loeschen des referenzierten
// Objektes flexibler zu gestalten.
//
template <class ObjType>
void deleteRefPtr(ObjType *ptr)
{
    delete ptr;
}

//======================================================
// Die Reference-Pointer Klasse:
//
// Beim Referenzieren eines Objektes ueber einen SmartPointer wird ein Zaehler fuer die referezierte Objekt-
// adresse inkrementiert. Wird der Pointer wieder einem anderen Objekt zugewiesen, so wird der Zaehler fuer das
// urspruenglich referenzierte Objekt wieder dekremtiert, ebenfalls beim Destruktor des Reference-Pointers.
// Tatsaechlich geloescht wird das referenzierte Objekt erst, wenn der zugehoerige Zaehler auf 0 ist. Dies geschieht
// ueber die globale Template-Funktion deleteRefPtr(), die bei Bedarf ueberschrieben werden kann.
// Der Reference-Pointer verfuegt also sozusagen ueber eine automatische Garbage Collection

template <class ObjType>
class MbSmartPtr : public MbSmartPtrBase
{
public:
    // Konstruktoren //bei explicit geht der implizite cast nicht mehr, aber um keinen stress zu verursachen
    /*explicit*/ MbSmartPtr<ObjType>(const ObjType *pPtr = NULL) : MbSmartPtrBase(), mpPtr(NULL) { init(pPtr); }
    template <class ParamType>
    MbSmartPtr<ObjType>(const MbSmartPtr<ParamType> &ptr) : MbSmartPtrBase(), mpPtr(NULL)
    {
        init(ptr.get());
    }
    // Destruktor
    ~MbSmartPtr<ObjType>() override { init(NULL); }
    //---------------------------------------------------
    // Kopierkonstruktor
    MbSmartPtr<ObjType>(const MbSmartPtr<ObjType> &ptr) : MbSmartPtrBase(), mpPtr(NULL) { init(ptr.get()); }
    //---------------------------------------------------
    // Zuweisungsoperatoren
    template <class ParamType>
    const MbSmartPtr<ObjType> &operator=(const MbSmartPtr<ParamType> &ptr)
    {
        init(ptr.get());
        return *this;
    }
    const MbSmartPtr<ObjType> &operator=(const MbSmartPtr<ObjType> &ptr)
    {
        init(ptr.get());
        return *this;
    }

    const MbSmartPtr<ObjType> &operator=(const ObjType *pPtr)
    {
        init(pPtr);
        return *this;
    }
    //---------------------------------------------------
    // Dereferenzierung-Operatoren
    ObjType &operator*() const { return *mpPtr; }
    ObjType *operator->() const { return mpPtr; }
    bool operator!() const { return !mpPtr; }
    operator ObjType *() const { return mpPtr; }
    //---------------------------------------------------
    // Methoden
    ObjType *get() const { return mpPtr; }
    //---------------------------------------------------
    int ref_count() const { return MbSmartPtrBase::ref_count(mpPtr); }
    //---------------------------------------------------
    bool release() const { return MbSmartPtrBase::removeFromGC(mpPtr); }

private:
    void init(const ObjType *pPtr)
    {
        // Nur was tun, wenn wirklich noetig
        if (pPtr == mpPtr)
            return;

        // Aktuell referenziertes Objekt freigeben, dabei ueberpruefen, ob letztes Release
        if (mpPtr && releaseRef(mpPtr)) {
            // referenziertes Objekt loeschen
            deleteRefPtr(mpPtr);
        }

        // Wenn pPtr ein neues Objekt ist, Zugriffszaehler auf neues Objekt erhoehen
        mpPtr = const_cast<ObjType *>(pPtr);
        if (mpPtr)
            addRef(mpPtr);
    }

private:
    ObjType *mpPtr;
};

#endif // MBSMARTPTR_H

//! \}

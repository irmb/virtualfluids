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
#ifndef MBSMARTPTRBASE_H
#define MBSMARTPTRBASE_H

#include <iostream>
#include <map>

//============================================================
// Klasse MbSmartPtrBase
//
// Basisklasse, speziell fuer MbSmartPtr, die das eigentliche
// Reference-Counting uebernimmt.
//
class MbSmartPtrBase
{
    // Ursprung:
    // mpCntrMap ist ein Pointer, weil sichergestellt sein muss, dass die
    // Map existiert, wenn das erste mal darauf zugegriffen wird.
    // Ein Zugriff zwischen zwei statischen Objekten kann zum Fehler fuehren, da
    // die Reihenfolge der Konstruktorenaufrufe dann vom Linker bestimmt wird.

    // Anpassung a la UbWriter mit SingletonMap
    class MbSmartPtrBaseMap
    {
    private:
        MbSmartPtrBaseMap() = default;

        std::map<void *, int> mpCntrMap;

    public:
        MbSmartPtrBaseMap(const MbSmartPtrBaseMap &) = delete;
        const MbSmartPtrBaseMap &operator=(const MbSmartPtrBaseMap &) = delete;

        static MbSmartPtrBaseMap *getInstance()
        {
            static MbSmartPtrBaseMap instance;
            return &instance;
        }
        std::map<void *, int> &getMap() { return mpCntrMap; }
    };

protected:
    MbSmartPtrBase()          = default;
    virtual ~MbSmartPtrBase() = default;
    bool addRef(void *p);
    bool releaseRef(void *p);
    bool removeFromGC(void *ptr) const;
    int ref_count(void *ptr) const;
};

#endif // MBSMARTPTRBASE_H

//! \}

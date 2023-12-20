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
//! \addtogroup utilities
//! \ingroup basics
//! \{
//! \author Soeren Freudiger, Sebastian Geller
//=======================================================================================
#ifndef UBOBSERVER_H
#define UBOBSERVER_H

class UbObservable;

//////////////////////////////////////////////////////////////////////////
//!
//! \brief Observer
//! \details This interface must be implemented by classes which want to
//! observe other objects.
//! IMPORTANT: if you delete an observer, ensure to remove Observer from
//!            all his observed observable objects before!!!
//! example: see end of UbObservable.h-file
//!
//////////////////////////////////////////////////////////////////////////

class UbObserver
{
protected:
    UbObserver() = default;

public:
    virtual ~UbObserver() = default;

    /*======================================================================*/
    /*  Methods                                                           */
    /*                                                                      */
    /*!
    This function is called when the observable indicated that an object
    has changed.
    \param changedObject Object which has changed
    */
    virtual void objectChanged(UbObservable *changedObject) = 0;
    /*!
    This function is called when the observable indicated that an object
    should be deleted.
    \param objectForDeletion Object which should be deleted
    */
    virtual void objectWillBeDeleted(UbObservable *objectForDeletion) = 0;
};

#endif

//! \}

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
#ifndef UBOBSERVABLE_H
#define UBOBSERVABLE_H

#include <iostream>
#include <list>

#include <basics/utilities/UbObserver.h>

class UbObserver;

//////////////////////////////////////////////////////////////////////////
//!
//! \brief Observable object
//! \details This class provides Observables. The Observeres which observe this
//!  Observable are stored in an observerlist.
//!  IMPORTANT: objectWillBeDeleted is called at UbObserver::~UbObserver
//!             this destructor is called AFTER the destructor of the
//!             child classes. if you down_cast the pointer sent with the
//!             objectWillBeDeleted(UbObserver* objpointer) then have to go this:
//!
//!               if(dynamic_cast<UbObserver*>(observedObj)==objpointer)
//!                     (e.g.) observedObj=NULL;
//!   example: see end of file
//!
//!   a copy of an UbservableObject will NOT copy the observerList
//!  <UL>
//!    <LI><B>Extending:</B> This UbObservable is the observable object itself. Extending should be used
//!    where object types can be extended from UbObservable.
//!    <LI><B>Associating:</B> Initialization is done via the constructor <tt>UbObservable(ObservableObject)</tt>.
//!    Associating may be used, where object types to be observed could not be extended from UbObservable.
//!  </UL>
//!
//! see UbObserver
//!
//////////////////////////////////////////////////////////////////////////

class UbObservable
{
protected:
    /*======================================================================*/
    /*  Konstruktoren                                                       */
    /*                                                                      */
    /**
      Creates a UbObservable itself to be the object to be observed.
      Usually this constructor is used in extended classes.
    */
    UbObservable() = default;

    UbObservable(const UbObservable & /*src*/)
    {
        // no copy of observers !!!
    }

    // falls irgendein schlaumeier den =operator von UbObservable aufrufen sollte,
    // dann macht diesr auch keine kopie! (Allg: zuweisungsoperatoren werden nie vererbt
    UbObservable &operator=(const UbObservable & /*src*/) { return *this; }

    //   /**
    //     Creates a UbObservable for the specified Object to be observed.
    //     Usually this constructor is used in associating UbObservable.
    //     @param object Object to be observed
    //   */
public:
    /*======================================================================*/
    /*  Destruktor                                                          */
    /*                                                                      */
    /**
     */
    virtual ~UbObservable() { this->notifyObserversObjectWillBeDeleted(); }

    /*======================================================================*/
    /*  methods                                                            */
    /*                                                                      */

    /**
    Adds an UbObserver to the observerlist.
    @param observer the observer to add to this observable (note that an observer may observe c1 observable more than
    once)
    */
    virtual void addObserver(UbObserver *observer)
    {
        if (!observer)
            return;
        for (std::list<UbObserver *>::iterator pos = mObserverList.begin(); pos != mObserverList.end(); ++pos) {
            if (*pos == observer)
                return;
        }
        this->mObserverList.push_back(observer);
    }
    /**
    Deletes an UbObserver from the observerlist.
    @param observer the observer to remove from this observable (note that all observers identical are deleted)
    ( delete means delete Heap... but here we're only removing a pointer)
    */
    virtual void removeObserver(UbObserver *observer)
    {
        if (!observer)
            return;
        this->mObserverList.remove(observer);
    }
    /**
    Deletes all Observers from the observerlist.
    ( delete means delete Heap... but here we're only removing a pointer)
    */
    virtual void removeAllObservers() { this->mObserverList.clear(); }

    /**
      Checks whether the specified UbObserver observes this observable.
      @param observer the observer to remove from this observable (note that all observers identical are deleted)
      @return true if the specified observer observes this observable
    */
    virtual bool isObservedBy(UbObserver *observer)
    {
        if (!observer)
            return false;
        for (std::list<UbObserver *>::iterator pos = mObserverList.begin(); pos != mObserverList.end(); ++pos) {
            if (*pos == observer)
                return true;
        }
        return false;
    }
    /**
      Notifies all of its observers that something happened. Does nothing, if the observed object is null.
      Calls the Method UbObserver.objectChanged(Object) with the object of this observable as parameter.
      The Method UbObserver.objectChanged(Object) must be defined
      by each class implementing the interface TiObserver
    */
    virtual void notifyObserversObjectChanged()
    {
        std::list<UbObserver *>::iterator tmp_pos; // es kann sein, dass der aktuelle observer waehrend
                                                   // objectChanged() removed wird...
        for (std::list<UbObserver *>::iterator pos = mObserverList.begin(); pos != mObserverList.end();) {
            // cout<<"in notifyObserversObjectChanged\n";
            // cout<<this->mObserverList.size()<<endl;

            tmp_pos = pos++; // erst tmp_pos=pos und dann pos++
            (*tmp_pos)->objectChanged(this);
        }
    }

    std::list<UbObserver *> *getObserverList() { return &mObserverList; }

    virtual std::string toString() { return "UbObservable - toString()"; }

private:
    /**
      Notifies all of its observers that something happened. Does nothing, if the observed object is null.
      Calls the Method UbObserver.objectChanged(Object) with the object of this observable as parameter.
      The Method UbObserver.objectChanged(Object) must be defined
      by each class implementing the interface TiObserver
    */
    virtual void notifyObserversObjectWillBeDeleted()
    {
        std::list<UbObserver *>::iterator tmp_pos; // es kann sein, dass der aktuelle observer waehrend
                                                   // objectWillBeDeleted() removed wird...
        for (std::list<UbObserver *>::iterator pos = mObserverList.begin(); pos != mObserverList.end();) {
            // cout<<"in notifyObserversObjectWillBeDeleted\n";
            // cout<<this->mObserverList.size()<<endl;

            tmp_pos = pos++;
            (*tmp_pos)->objectWillBeDeleted(this);
        }
    }

    std::list<UbObserver *> mObserverList;
};
/*=========================================================================*/

#ifdef RCF_USE_SF_SERIALIZATION
SF_NO_CTOR(UbObservable);
#endif // RCF_USE_SF_SERIALIZATION

#endif

////  E X A M P L E
////===================
// class Point : public UbObservable
//{
// public:
//   Point(){x=y=0;}
//   ~Point(){}
//   void setXCorrdinates(int x,int y)
//   {
//     this->x = x; this->y = y;
//     this->notifyObserverObjectChanged();
//   }
// private:
//   int x,y;
//};
// class VisPoint : public UbObserver
//{
// public:
//   VisPoint(Point* point)
//   {
//      this->point = point;
//      this->point->addObserver(this);
//   }
//   ~VisPoint()
//   {
//      if(this->point) point->removeObserver(this);
//   }
//   void update() { /* do some actualisation stuff */ }
//   void objectChanged(UbObservable* changedObject)
//   {
//      Point* point = dynamic_cast<Point*>(changedObject);
//      if( !this->point || this->point != point ) return;
//      this->repaint();
//   }
//   void objectWillBeDeleted(UbObservable* objectForDeletion)
//   {
//      if(!this->point) return;
//      UbObservable* obsobjet = dynamic_cast<UbObservable*>(this->point);
//      if(obsobjet == objectForDeletion) this->point = NULL;
//      ///////////////////////////////////////////////////////////////////
//      //*********************************************************************//
//      //INGEGEN erster annahmen nicht verwenden, da es nicht immer funktioniert
//      //z.B. bei mehrfachvererbung haut es nicht hin!
//      ////      Point* point = reinterpret_cast<point*>(objectForDeletion);
//      ////if(!this->point || objectForDeletion != this->point) return;
//      ////this->point = NULL;
//      //*********************************************************************//
//      //was hingegen immer moeglich sein sollte:
//      //if(dynamic_cast<void*>(objectForDeletion)==dynamic_cast<void*>(this->point))
//   }
// private:
//   Point* point;
//};

//! \}

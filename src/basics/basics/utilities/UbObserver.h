//  _    ___      __              __________      _     __
// | |  / (_)____/ /___  ______ _/ / ____/ /_  __(_)___/ /____
// | | / / / ___/ __/ / / / __ `/ / /_  / / / / / / __  / ___/
// | |/ / / /  / /_/ /_/ / /_/ / / __/ / / /_/ / / /_/ (__  )
// |___/_/_/   \__/\__,_/\__,_/_/_/   /_/\__,_/_/\__,_/____/
//
#ifndef UBOBSERVER_H
#define UBOBSERVER_H

class UbObservable;
/*=========================================================================*/
/*  Observer                                                               */
/*                                                                         */
/**
This interface must be implemented by classes which want to
observe other objects.
IMPORTANT: if you delete an observer, ensure to remove Observer from
           all his oberved observable objects before!!!
example: see end of UbObservable.h-file
<BR><BR><HR>
@author <A HREF="mailto:geller@cab.bau.tu-bs.de">S. Geller</A>
@author <A HREF="mailto:muffmolch@gmx.de">S. Freudiger</A>
@version 1.1 - 20.11.04
*/

class UbObserver 
{
protected:

   UbObserver(){}

public:

   virtual ~UbObserver(){}

   /*======================================================================*/
   /*  Methoden                                                            */
   /*                                                                      */
   /**
   This function is called when the observable indicated that an object
   has changed.
   @param changedObject Object which has changed
   */
   virtual void objectChanged(UbObservable* changedObject)=0;
   /**
   This function is called when the observable indicated that an object
   should be deleted.
   @param objectForDeletion Object which should be deleted
   */
   virtual void objectWillBeDeleted(UbObservable* objectForDeletion)=0;
};

#endif



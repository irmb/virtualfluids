//  _    ___      __              __________      _     __
// | |  / (_)____/ /___  ______ _/ / ____/ /_  __(_)___/ /____
// | | / / / ___/ __/ / / / __ `/ / /_  / / / / / / __  / ___/
// | |/ / / /  / /_/ /_/ / /_/ / / __/ / / /_/ / / /_/ (__  )
// |___/_/_/   \__/\__,_/\__,_/_/_/   /_/\__,_/_/\__,_/____/
//
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
   //Ursprung:
   // mpCntrMap ist ein Pointer, weil sichergestellt sein muss, dass die
   // Map existiert, wenn das erste mal darauf zugegriffen wird.
   // Ein Zugriff zwischen zwei statischen Objekten kann zum Fehler fuehren, da
   // die Reihenfolge der Konstruktorenaufrufe dann vom Linker bestimmt wird.

   //Anpassung a la UbWriter mit SingletonMap
   class MbSmartPtrBaseMap
   {
   private:
      MbSmartPtrBaseMap() = default;
      MbSmartPtrBaseMap( const MbSmartPtrBaseMap& );                  //no copy allowed
      const MbSmartPtrBaseMap& operator=( const MbSmartPtrBaseMap& ); //no copy allowed

      std::map<void*,int> mpCntrMap;
   public:
      static MbSmartPtrBaseMap* getInstance() { static MbSmartPtrBaseMap instance; return &instance; }
      std::map<void*,int>& getMap()           { return mpCntrMap;                                    }
   };

protected:
   MbSmartPtrBase() = default;
   virtual ~MbSmartPtrBase() = default;
   bool addRef(void* p);
	bool releaseRef(void* p);
   bool removeFromGC(void* ptr) const;
   int  ref_count(void* ptr) const;
};

#endif //MBSMARTPTRBASE_H

//  _    ___      __              __________      _     __
// | |  / (_)____/ /___  ______ _/ / ____/ /_  __(_)___/ /____
// | | / / / ___/ __/ / / / __ `/ / /_  / / / / / / __  / ___/
// | |/ / / /  / /_/ /_/ / /_/ / / __/ / / /_/ / / /_/ (__  )
// |___/_/_/   \__/\__,_/\__,_/_/_/   /_/\__,_/_/\__,_/____/
//
#ifndef UBPOINTERWRAPPER_H
#define UBPOINTERWRAPPER_H

//kappselt dynamische Objekte zur remote uebetragung
//bei RCF werden z.B. aufgrund GC alle lokalen Objekte und 
//"nackte" Pointer die automatisch als shared_ptr initialisert 
//werde nach Methoden-Aufruf zerstoert
//hierfuer kann man dann den UbPointerWrapper verwenden

template<typename T>
class UbPointerWrapper
{
public:
	UbPointerWrapper() : pointer(NULL) {}
	
	UbPointerWrapper(T* pointer) : pointer(pointer) {}

   T* get() { return pointer; }

   template<class Archive>
	void SF_SERIALIZE(Archive & ar) 
   {
		ar & pointer;
	}

private:
   T* pointer;
};

#endif //UBPOINTERWRAPPER_H

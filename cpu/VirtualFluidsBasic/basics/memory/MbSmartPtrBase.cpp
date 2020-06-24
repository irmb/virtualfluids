#include <basics/memory/MbSmartPtrBase.h>

using namespace std;

bool MbSmartPtrBase::addRef(void* ptr)
{
   MbSmartPtrBaseMap::getInstance()->getMap()[ptr]++;
	return true;
}
//-------------------------------------------------
bool MbSmartPtrBase::releaseRef(void* ptr)
{
   map<void*,int>& ptrMap = MbSmartPtrBaseMap::getInstance()->getMap();
   map<void*,int>::iterator pos=ptrMap.find(ptr);
	
   if( pos!=ptrMap.end() )
	{
		pos->second--;
		
      if(pos->second==0)
		{
			ptrMap.erase(pos);
			return true;
		}
	}
	return false;
}
//-------------------------------------------------
bool MbSmartPtrBase::removeFromGC(void* ptr) const 
{
   if( MbSmartPtrBaseMap::getInstance()->getMap().erase(ptr) ) return true;
   return false;
}
//-------------------------------------------------
int MbSmartPtrBase::ref_count(void* ptr) const 
{
   map<void*,int>& ptrMap = MbSmartPtrBaseMap::getInstance()->getMap();
   map<void*,int>::iterator pos=ptrMap.find(ptr);

   if( pos!=ptrMap.end() ) return pos->second;
   else                    return 0;
}



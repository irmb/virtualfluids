//  _    ___      __              __________      _     __
// | |  / (_)____/ /___  ______ _/ / ____/ /_  __(_)___/ /____
// | | / / / ___/ __/ / / / __ `/ / /_  / / / / / / __  / ___/
// | |/ / / /  / /_/ /_/ / /_/ / / __/ / / /_/ / / /_/ (__  )
// |___/_/_/   \__/\__,_/\__,_/_/_/   /_/\__,_/_/\__,_/____/
//
#ifndef MBMEMPOOL_H
#define MBMEMPOOL_H

#include <queue>
#include <list>


template <typename TData, int asize>
class MbMemPool
{

protected:
   MbMemPool(){}
   // Alle 3 Attribute sind Singleton-Objekte !
   // Allokiere Blocke der Groesse m_asize
   static int m_asize;       
   // Halte alle freien Pointer (jedes einzelne Element)  in eine FIFO Liste
   static std::queue<void*> m_queue;
   // Pointer auf Bloecke zum Loeschen !
   static std::list<TData*> m_list;

public:

   
   ~MbMemPool(){}

   // Daten-Bloecke Loeschen, damit wird der gesamte Speicher freigegeben,
   // erst aufrufen, wenn die objekte nicht mehr gebraucht werden!
   static void	deallocatePool();

   void* operator new(std::size_t size)
   {
      void*  pNew;
      TData* feld;	
      int i;

      //i=m_queue.size();
      //pNew = m_queue.front();
      if(m_queue.size()==0) 
      {
         //Wenn kein freier Speicher mehr vorhanden, Block anlegen
         feld = new TData[m_asize];
         m_list.push_back(feld);
         for(i=0 ; i<m_asize ; i++)
         {
            pNew = (void*) &(feld[i]);
            m_queue.push( pNew );
         }
      }
      pNew = m_queue.front();
      m_queue.pop();
      return pNew;

   }

   void  operator delete(void* p)
   {
      m_queue.push(p);
   }
};


template <typename TData, int asize> 
std::queue<void*>  MbMemPool<TData,asize>::m_queue;

template <typename TData, int asize> 
std::list<TData*>  MbMemPool<TData,asize>::m_list;

template <typename TData, int asize> 
int  MbMemPool<TData,asize>::m_asize=asize;

template <typename TData, int asize> 
void MbMemPool<TData,asize>::deallocatePool()
{	
   for(typename std::list<TData*>::iterator pos=m_list.begin() ; pos!=m_list.end(); ++pos)
   {
      delete[] pos;
   }
}

#endif //MBMEMPOOL_H

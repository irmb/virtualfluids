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
#ifndef CBVECTORPOOL_H
#define CBVECTORPOOL_H

#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <typeinfo>
#include <vector>

#include <basics/container/CbVector.h>
#include <basics/utilities/UbException.h>
#include <basics/utilities/UbLogger.h>
#include <basics/utilities/UbTuple.h>

//#include "MPICommunicator.h"
//
//#include <execinfo.h>
//#include <stdio.h>
//#include <stdlib.h>
//#include <unistd.h>

/*=========================================================================*/
/*  CbVectorPool                                                               */
/*                                                                         */
/**
<BR><BR>
@author <A HREF="mailto:muffmolch@gmx.de">S. Freudiger</A>
@version 1.0 - 08.11.07
@version 1.1 - 09.02.08
*/

/*
Durch Verwendung eines CbVectors in Verbindung mit einem CbVectorAllocatorPool
wird der Datenvector nicht direkt im CbVector gehalten, sondern ist ein Teil
des Datenvectors des Uebergabe-CbVectorPools.
Die Methoden der von CbVectors funktionieren fehlerfrei
Es mss einem jedoch bewusst sein, dass die "resize"-Methoden l�nger ben�tigen, da
u.U. viele Elemente im Speicher verschoeben werden muessen.
Der Poolvector enthaelt KEINE gaps, so dass er z.B. gut zur �bertragung via MPI
geeignet ist...

Verhaltensweise bei Zerstoeren des Pools:
wird der Pool zerst�rt bevor man die CbVectoren zerst�rt, so wird beim n�chsten
Datenzugriffsversuch eine entsprechende Exception geworfen, denn alle DatenElemente
des CbVEctors werden restet und der Pool dort zu NULL gesetzt.

Verhaltensweise bei Zerstoeren eines CbVectors:
hier ganz normal der Datenspeicher wieder freigegen und der Poolvektor verk�rzt
*/

//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

template <typename T>
class CbVectorAllocatorPool;

/*==================================================================*/
template <typename T>
class CbVectorPool
{
public:
    using value_type = typename CbVector<T>::value_type;
    using size_type  = typename CbVector<T>::size_type;
    using Pool       = std::vector<value_type>;

    using CbVectorKey     = std::string;
    using CbVectorMap     = std::map<CbVectorKey, CbVector<value_type> *>;
    using CbVectorMapIter = typename CbVectorMap::iterator;

public:
    //////////////////////////////////////////////////////////////////////////
    CbVectorPool(const size_type &startPoolSize = 20000) // startPoolSize*sizeof(T)/1024/1024 [MB]
        : poolStartAdress(NULL), nextCbVectorStartIndexInPool(0), nextCbVectorKey()
    {
        pool.reserve(startPoolSize);
    }
    /*==================================================================*/
    virtual ~CbVectorPool()
    {
        // hier werden lediglich ihre datenvektoren "resetet"
        for (CbVectorMapIter it = cbVectorMap.begin(); it != cbVectorMap.end(); ++it) {
            CbVector<value_type> &vec = *it->second;
            CbVectorAllocatorPool<value_type> &allocator =
                dynamic_cast<CbVectorAllocatorPool<value_type> &>(*vec.getAllocator());
            // FIXME: //if(allocator.ptrVectorPool != this) UB_THROW( UbException(UB_EXARGS,"CbVectorAllocator is part
            // of different Pool") );

            // allocator daten reseten
            allocator.ptrVectorPool    = NULL;
            allocator.key              = CbVectorKey();
            allocator.startIndexInPool = 0;

            // Datenzeiger/-groessen reseten
            allocator.ptrDataOf(vec)  = NULL;
            allocator.dataSizeOf(vec) = 0;
        }
    }
    /*==========================================================*/
    CbVectorKey getNextCbVectorKey() const { return this->nextCbVectorKey; }
    /*==================================================================*/
    bool allocVectorData(CbVector<value_type> &vec, const size_type &dataSize, const value_type &value = value_type())
    {
        // pool-allocator holen
        CbVectorAllocatorPool<value_type> &allocator =
            dynamic_cast<CbVectorAllocatorPool<value_type> &>(*vec.getAllocator());
        if (allocator.ptrVectorPool != this)
            UB_THROW(UbException(UB_EXARGS, "CbVectorAllocator is part of different Pool"));

        // alloc nur wenn cbVector noch kein Element von Pool!
        if (cbVectorMap.find(allocator.key) == cbVectorMap.end()) {
            return this->allocData(allocator, vec, dataSize, value);
        }

        UB_THROW(UbException(UB_EXARGS, "vector-key=" + UbSystem::toString(allocator.key) + " already taken! (e.g. SetConnectorBlockVisitor was called several times"));
    }
    /*==================================================================*/
    bool resizeVectorData(CbVector<value_type> &vec, const size_type &dataSize, const value_type &value = value_type())
    {
        CbVectorAllocatorPool<value_type> &allocator =
            dynamic_cast<CbVectorAllocatorPool<value_type> &>(*vec.getAllocator());
        if (allocator.ptrVectorPool != this)
            UB_THROW(UbException(UB_EXARGS, "CbVectorAllocator is part of different Pool"));

        // cbVector noch nicht in map?
        CbVectorMapIter pos = cbVectorMap.find(allocator.key);

        if (pos != cbVectorMap.end()) // cbVector vorhanden
        {
            // wenn bei alloc keine Laenge zugewiesen wurde, so erfolgt das nun
            if (allocator.startIndexInPool == 0 && allocator.ptrDataOf(vec) == NULL)
                return this->allocData(allocator, vec, dataSize, value);
            else
                return this->resizeData(allocator, vec, dataSize, value);
        }

        UB_THROW(UbException(UB_EXARGS,
                             "vector gehoert laut allocator zum pool aber den key gibt s nicht... wie kann das sein?"));
    }
    /*==================================================================*/
    bool deallocVectorData(CbVector<value_type> &vec)
    {
        CbVectorAllocatorPool<value_type> &allocator =
            dynamic_cast<CbVectorAllocatorPool<value_type> &>(*vec.getAllocator());
        if (allocator.ptrVectorPool != this)
            UB_THROW(UbException(UB_EXARGS, "CbVectorAllocator is part of different Pool"));

        // nur wenn vector auch teil des
        if (cbVectorMap.erase(allocator.key) > 0) {
            if (this->resizeData(allocator, vec, 0, 0)) {
                allocator.ptrVectorPool    = NULL;
                allocator.key              = CbVectorKey();
                allocator.startIndexInPool = 0;

                // das Datenzeiger/-groessen reseten wird bereits in resize durchgefuehrt
                return true;
            } else
                UB_THROW(UbException(UB_EXARGS, "unknown error"));
        }

        // SPtr<Communicator> comm = MPICommunicator::getInstance();
        // int myid = comm->getProcessID();

        //      // Get the name of the processor
        //      char machinename[MPI_MAX_PROCESSOR_NAME];
        //      int name_len;
        //      MPI_Get_processor_name(machinename, &name_len);
        //      UBLOG(logINFO, "PID = " << myid << " host name: " << machinename);
        //
        //      int j, nptrs;
        //#define SIZE 100
        //      void *buffer[100];
        //      char **strings;
        //
        //      nptrs = backtrace(buffer, SIZE);
        //      printf("backtrace() returned %d addresses\n", nptrs);
        //
        //      /* The call backtrace_symbols_fd(buffer, nptrs, STDOUT_FILENO)
        //      would produce similar output to the following: */
        //
        //      strings = backtrace_symbols(buffer, nptrs);
        //      if (strings == NULL)
        //      {
        //         perror("backtrace_symbols");
        //         exit(EXIT_FAILURE);
        //      }
        //
        //      for (j = 0; j < nptrs; j++)
        //         printf("%s\n", strings[j]);
        //
        //      free(strings);

        UB_THROW(UbException(UB_EXARGS,
                             "vector gehoert laut allocator zum pool aber den key gibt s nicht... wie kann das sein?"));
    }
    /*==================================================================*/
    friend std::ostream &operator<<(std::ostream &os, const CbVectorPool &cbPool)
    {
        os << "map" << std::endl;
        for (CbVectorMapIter pos = cbPool.cbVectorMap.begin(); pos != cbPool.cbVectorMap.end(); ++pos) {
            CbVectorAllocatorPool<value_type> &tmpAllocator =
                dynamic_cast<CbVectorAllocatorPool<value_type> &>(*pos->second->getAllocator());
            os << "vector-size=" << pos->second->size() << "vector-Adress=" << tmpAllocator->ptrDataOf(*pos->second)
               << ", allocator(key=" << tmpAllocator.key << ", startIndex=" << tmpAllocator.startIndexInPool << ")"
               << std::endl;
            for (size_type i = 0; i < pos->second->size(); i++)
                os << (*pos->second)[i] << ",";
            os << std::endl;
        }
        os << "pool" << std::endl;
        for (size_type i = 0; i < cbPool.pool.size(); i++) {
            os << cbPool.pool[i] << ",";
            os << std::endl;
        }

        return os;
    }
    /*==================================================================*/
    typename CbVectorMap::size_type getNofStoredVectors() const { return this->cbVectorMap.size(); }
    /*==================================================================*/
    typename Pool::size_type getPoolSize() const { return this->pool.size(); }
    /*==================================================================*/
    // checks if all vectors have one to one pool-entries
    bool consistencyCheck()
    {
        std::vector<int> pool2(pool.size(), 0);
        for (CbVectorMapIter it = cbVectorMap.begin(); it != cbVectorMap.end(); ++it) {
            CbVector<value_type> &tmpVec = *it->second;
            CbVectorAllocatorPool<value_type> &tmpAllocator =
                dynamic_cast<CbVectorAllocatorPool<value_type> &>(*tmpVec.getAllocator());
            for (size_type i = tmpAllocator.startIndexInPool; i < tmpAllocator.startIndexInPool + tmpVec.size(); ++i) {
                pool2.at(i)++;
            }
        }
        for (size_type i = 0; i < pool2.size(); ++i) {
            if (pool2.at(i) > 1) {
                UBLOG(logERROR, UB_FUNCTION << " - test failed typo 1");
                return false;
            }
            if (pool2.at(i) < 1) {
                UBLOG(logERROR, UB_FUNCTION << " - test failed typo 2");
                return false;
            }
        }
        return true;
    }

protected:
    /*==================================================================*/
    inline bool allocData(CbVectorAllocatorPool<value_type> &allocator, CbVector<value_type> &vec,
                          const size_type &dataSize, const value_type &value)
    {
        // safety checks
        if (allocator.startIndexInPool != 0 || allocator.ptrDataOf(vec) != NULL || allocator.dataSizeOf(vec) != 0) {
            UB_THROW(UbException(UB_EXARGS, "zu allokierender vector ist nicht ganz sauber!!"));
        }

        // poolVector vergroessern
        if (dataSize > 0) {
            pool.resize(pool.size() + dataSize, value);

            // Zeiger der vorhandenen CbVectoren neu setzen, wenn Pool im Speicher verschoben wurde
            if (poolStartAdress != &pool.front()) {
                poolStartAdress = &pool.front();
                for (CbVectorMapIter it = cbVectorMap.begin(); it != cbVectorMap.end(); ++it) {
                    CbVector<value_type> &tmpVec = *it->second;
                    CbVectorAllocatorPool<value_type> &tmpAllocator =
                        dynamic_cast<CbVectorAllocatorPool<value_type> &>(*tmpVec.getAllocator());

                    if (!tmpAllocator.ptrDataOf(tmpVec))
                        continue; // Fall: CbVector hat noch keinen Datenbereich (data zeigt auf NULL)
                    tmpAllocator.ptrDataOf(tmpVec) = &pool[tmpAllocator.startIndexInPool];
                }
                // std::cout<<"CbVectorPoolMpi::allocVectorData vector wurde im speicher verschoben - adressen
                // angepasst!!!"<<std::endl;
            }

            // aktuellem element adresse zuweisen (wurde evtl schon inder schleife zuvor gemacht)
            allocator.ptrDataOf(vec)   = &pool.at(nextCbVectorStartIndexInPool);
            allocator.startIndexInPool = nextCbVectorStartIndexInPool;

            // neuen StartIndex fuer naechstes Element berechnen
            nextCbVectorStartIndexInPool += dataSize;
            if (nextCbVectorStartIndexInPool != pool.size())
                UB_THROW(UbException(UB_EXARGS, "index Problem... Annahme falsch?"));
        }

        // vector zu map hinzuf�gen (speicher wird dann anschliessend zugwiesen)
        cbVectorMap.insert(
            std::make_pair(allocator.key, &vec)); // ist angeblich performanter als  cbVectorMap[ allocator.key ] =
                                                  // cbVector; //aus Effective STL von Scott Meyer
        allocator.dataSizeOf(vec) = dataSize;

        // dummDoof nextKey-Generung...
        if (allocator.key >= this->nextCbVectorKey)
            this->nextCbVectorKey = allocator.key + "1";

        return true;
    }
    /*==========================================================*/
    bool resizeData(CbVectorAllocatorPool<value_type> &allocator, CbVector<value_type> &vec, const size_type &dataSize,
                    const value_type &value)
    {
        // datenvector verlaengern/-kuerzen
        typename Pool::iterator startPos =
            pool.begin() + allocator.startIndexInPool; // startPosition der cbVector-Daten im Pool
        if (vec.size() > dataSize)
            pool.erase(startPos + dataSize, startPos + vec.size());
        else
            pool.insert(startPos + vec.size(), dataSize - vec.size(), value);

        //////////////////////////////////////////////////////////////////////////
        // adressen und laengen der einzelnen vectoren anpassen
        if (!pool.empty()) {
            bool poolMoved  = (poolStartAdress != &pool.front());
            poolStartAdress = &pool.front();

            for (CbVectorMapIter it = cbVectorMap.begin(); it != cbVectorMap.end(); ++it) {
                CbVector<value_type> &tmpVec = *it->second;

                if (tmpVec.size() > 0) {
                    CbVectorAllocatorPool<value_type> &tmpAllocator =
                        dynamic_cast<CbVectorAllocatorPool<value_type> &>(*tmpVec.getAllocator());
                    // liegt CbVector VOR ver�ndertem CbVector?
                    if (tmpAllocator.startIndexInPool <=
                        allocator.startIndexInPool) // ja: anpassung NUR wenn pool verschoben wurde!
                    {
                        if (poolMoved && tmpVec.size() > 0)
                            tmpAllocator.ptrDataOf(tmpVec) = &pool[tmpAllocator.startIndexInPool];
                    } else // nein: -> Adresse + Index MUSS immer angepasst werden
                    {
                        tmpAllocator.startIndexInPool += dataSize - vec.size();
                        tmpAllocator.ptrDataOf(tmpVec) = &pool[tmpAllocator.startIndexInPool];
                    }
                }
            }
        } else // Sonderfall: alle Elemente haben Laenge 0 -> kein pool -> alle Feld-Adressen auf NULL setzen!
        {
            poolStartAdress = NULL;
            for (CbVectorMapIter it = cbVectorMap.begin(); it != cbVectorMap.end(); ++it) {
                CbVector<value_type> &tmpVec = *it->second;
                CbVectorAllocatorPool<value_type> &tmpAllocator =
                    dynamic_cast<CbVectorAllocatorPool<value_type> &>(*tmpVec.getAllocator());
                tmpAllocator.startIndexInPool = 0;
            }
        }

        // restliche Daten von cbVector + allocator aktualisieren
        allocator.dataSizeOf(vec) = dataSize;
        if (dataSize == 0) {
            allocator.ptrDataOf(vec)   = NULL;
            allocator.startIndexInPool = 0;
        }

        nextCbVectorStartIndexInPool = pool.size();

        return true;
    }

protected:
    /*==================================================================*/
    void getCbVectorData(const CbVector<value_type> &vec, CbVectorKey &vectorKey, size_type &startIndexInPool,
                         size_type &dataSize)
    {
        CbVectorAllocatorPool<value_type> &allocator =
            dynamic_cast<CbVectorAllocatorPool<value_type> &>(*vec.getAllocator());

        startIndexInPool = allocator.startIndexInPool;
        vectorKey        = allocator.key;
        dataSize         = vec.size();
    }
    /*==================================================================*/
    void setCbVectorData(CbVector<value_type> &vec, const CbVectorKey &vectorKey, const size_type &startIndexInPool,
                         const size_type &dataSize)
    {
        CbVectorAllocatorPool<value_type> &allocator =
            dynamic_cast<CbVectorAllocatorPool<value_type> &>(*vec.getAllocator());

        allocator.startIndexInPool = startIndexInPool;
        allocator.key              = vectorKey;
        allocator.dataSizeOf(vec)  = dataSize;
        allocator.ptrDataOf(vec)   = &this->pool[startIndexInPool];
    }
    /*==================================================================*/

    CbVectorMap cbVectorMap; // informationsmap fuer MPIData und zugewiesener vector

    Pool pool;                                             // globaler Datenvector
    typename Pool::pointer poolStartAdress;                // StartAdresse des aktuellen Datenvektors
    typename Pool::size_type nextCbVectorStartIndexInPool; // StartIndex fuer den naechsten CbVector

    // key - erstmal dummdoof
    CbVectorKey nextCbVectorKey;
};

//////////////////////////////////////////////////////////////////////////
//  CbVectorAllocatorPool
//////////////////////////////////////////////////////////////////////////
template <typename T>
class CbVectorAllocatorPool : public CbVectorAllocator<T>
{
public:
    // typedefs wiederholen, da Basisklasse = template -> "Dependent-Base"-Problem
    using value_type = typename CbVector<T>::value_type;
    using size_type  = typename CbVector<value_type>::size_type;

    friend class CbVectorPool<value_type>;

    CbVectorAllocatorPool(const CbVectorAllocatorPool &) = delete;
    const CbVectorAllocatorPool &operator=(const CbVectorAllocatorPool &) = delete;

public:
    /*==========================================================*/
    CbVectorAllocatorPool(const typename CbVectorPool<value_type>::CbVectorKey &key,
                          CbVectorPool<value_type> *const &ptrVectorPool)
        : CbVectorAllocator<value_type>(), key(key), startIndexInPool(0), ptrVectorPool(ptrVectorPool)
    {
        if (!ptrVectorPool)
            UB_THROW(UbException(UB_EXARGS, "ptrVectorPool==NULL"));
    }
    /*==========================================================*/
    // hier wird der key automatisch erzeugt!
    CbVectorAllocatorPool(CbVectorPool<value_type> *const &ptrVectorPool)
        : CbVectorAllocator<value_type>(), startIndexInPool(0), ptrVectorPool(ptrVectorPool)
    {
        if (!ptrVectorPool)
            UB_THROW(UbException(UB_EXARGS, "ptrVectorPool==NULL"));
        key = ptrVectorPool->getNextCbVectorKey();
    }
    /*==========================================================*/
    bool alloc(CbVector<value_type> &vec, const size_type &dataSize, const value_type &value = value_type()) override
    {
        if (!ptrVectorPool)
            UB_THROW(UbException(UB_EXARGS, "vectorPool seems to be destroyed, ptrVectorPool==NULL"));
        return ptrVectorPool->allocVectorData(vec, dataSize, value);
    }
    /*==========================================================*/
    bool resize(CbVector<value_type> &vec, const size_type &dataSize, const value_type &value = value_type()) override
    {
        if (!ptrVectorPool)
            UB_THROW(UbException(UB_EXARGS, "vectorPool seems to be destroyed, ptrVectorPool==NULL"));
        return ptrVectorPool->resizeVectorData(vec, dataSize, value);
    }
    /*==========================================================*/
    bool dealloc(CbVector<value_type> &vec) override
    {
        if (ptrVectorPool)
            return this->ptrVectorPool->deallocVectorData(vec);
        // wenn kein ptrVectorPool -> wurde bereits deallokiert
        return true;
    }
    /*==========================================================*/
    const CbVectorPool<value_type> &getCbVectorPool()
    {
        if (!ptrVectorPool)
            UB_THROW(UbException(UB_EXARGS, "vectorPool seems to be destroyed, ptrVectorPool==NULL"));
        return *ptrVectorPool;
    }
    /*==========================================================*/

private:
    typename CbVectorPool<value_type>::CbVectorKey key;
    typename CbVectorPool<value_type>::Pool::size_type startIndexInPool;

    CbVectorPool<value_type> *ptrVectorPool;
};

#endif // CBVECTORPOOL_H

//! \}

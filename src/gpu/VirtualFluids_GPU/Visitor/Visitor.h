#ifndef Visitor_H
#define Visitor_H

#include <vector>
#include <string>

#include "PointerDefinitions.h"
#include "Core/DataTypes.h"
#include "VirtualFluids_GPU_export.h"
#include "DataStructureInitializer/GridProvider.h"
#include "GPU/CudaMemoryManager.h"

#include <cassert>

class Parameter;


class VIRTUALFLUIDS_GPU_EXPORT Visitor
{
private: 
	SPtr<Parameter> para;


protected:
    Visitor()
    {
        this->updateInterval = 1;
    }

public:
    virtual ~Visitor() {}

    virtual void visit(Parameter* para, CudaMemoryManager* cudaManager, int level, unsigned int t)=0;
    virtual void init(Parameter* para, GridProvider* gridProvider, CudaMemoryManager* cudaManager)=0;

    ////////////////////////////////////////////////////////////////////////////
    /// \brief setUpdateInterval
    /// \param _updateInterval
    /// \note usage: setUpdateInterval(this->tout);
    /// \note usage: setUpdateInterval(div_ru(this->tout,10U));
    ///
    void setUpdateInterval(const uint &_updateInterval)
    {
        assert(_updateInterval>0);
        this->updateInterval = _updateInterval;
    }

    bool isDue(const uint &tLB) const
    {
        return (tLB%this->updateInterval==0);
    }

protected:
    uint updateInterval;                                                        ///< update interval in number of timesteps of the coars patch (1 = each time step)

};

#endif
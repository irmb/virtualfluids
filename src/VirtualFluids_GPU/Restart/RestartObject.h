#ifndef RestartObject_H
#define RestartObject_H

#include "LBM/LB.h"
#include "LBM/D3Q27.h"
#include "Parameter/Parameter.h"
#include "basics/utilities/UbSystem.h"
#include "boost/serialization/serialization.hpp"
#include "boost/serialization/vector.hpp"

class RestartObject
{
public:
   RestartObject();
   ~RestartObject();
   void load(Parameter* para);
   void safe(Parameter* para);
   void clear(Parameter* para);
protected:

private:
	//////////////////////////////////////////////////////////////////////////
	//Restart
	friend class boost::serialization::access;
	template<class Archive>
	void serialize(Archive & ar, const unsigned int version)
	{
		ar & fs;
	}
	//////////////////////////////////////////////////////////////////////////
	std::vector< std::vector<real> > fs;
	//////////////////////////////////////////////////////////////////////////

};

#endif

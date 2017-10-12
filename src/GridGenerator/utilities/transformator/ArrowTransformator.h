#ifndef ArrowTransformator_h
#define ArrowTransformator_h

#include <memory>
#include "GridGenerator/global.h"
#include "GridGenerator_EXPORT.h"

class Arrow;

class ArrowTransformator
{
public:
    static GridGenerator_EXPORT std::shared_ptr<ArrowTransformator> makeTransformator(doubflo delta, doubflo dx, doubflo dy, doubflo dz);
	virtual ~ArrowTransformator() {}

protected:
	ArrowTransformator() {}
	
public:
	virtual void transformGridToWorld(std::shared_ptr<Arrow> arrow) const = 0;
};


#endif

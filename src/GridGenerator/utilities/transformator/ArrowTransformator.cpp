#include "ArrowTransformator.h"

#include "TransformatorImp.h"


std::shared_ptr<ArrowTransformator> ArrowTransformator::makeTransformator(doubflo delta, doubflo dx, doubflo dy, doubflo dz)
{
    return std::shared_ptr<ArrowTransformator>(new TransformatorImp(delta, dx, dy, dz));
}

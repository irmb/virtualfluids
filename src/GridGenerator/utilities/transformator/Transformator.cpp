#include "Transformator.h"

#include "TransformatorImp.h"

std::shared_ptr<Transformator> Transformator::makeTransformator(doubflo delta, doubflo dx, doubflo dy, doubflo dz)
{
    return std::shared_ptr<Transformator>(new TransformatorImp(delta, dx, dy, dz));
}
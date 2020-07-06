#include "Transformator.h"

#include "utilities/transformator/TransformatorImp.h"

std::shared_ptr<Transformator> Transformator::makeTransformator(real delta, real dx, real dy, real dz)
{
    return std::shared_ptr<Transformator>(new TransformatorImp(delta, dx, dy, dz));
}
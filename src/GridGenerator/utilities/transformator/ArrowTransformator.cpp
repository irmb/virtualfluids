#include "ArrowTransformator.h"

#include "TransformatorImp.h"


std::shared_ptr<ArrowTransformator> ArrowTransformator::makeTransformator(real delta, real dx, real dy, real dz)
{
    return std::shared_ptr<ArrowTransformator>(new TransformatorImp(delta, dx, dy, dz));
}

#include "SparseMatrix4D.h"


size_t SparseMatrix4D::nx1 = 0;
size_t SparseMatrix4D::nx2 = 0;
size_t SparseMatrix4D::nx3 = 0;
size_t SparseMatrix4D::nx4 = 0;
//////////////////////////////////////////////////////////////////////////
SparseMatrix4D::SparseMatrix4D()
{

}
//////////////////////////////////////////////////////////////////////////
SparseMatrix4D::~SparseMatrix4D()
{

}
//////////////////////////////////////////////////////////////////////////
SparseMatrix4DPtr SparseMatrix4D::getInstance()
{
   if( !SparseMatrix4D::instance )
      SparseMatrix4D::instance = SparseMatrix4DPtr(new SparseMatrix4D());
   return SparseMatrix4D::instance;
}
//////////////////////////////////////////////////////////////////////////
void SparseMatrix4D::setDimensions(size_t nx1, size_t nx2, size_t nx3, size_t nx4) 
{
   SparseMatrix4D::nx1 = nx1;
   SparseMatrix4D::nx2 = nx2;
   SparseMatrix4D::nx3 = nx3;
   SparseMatrix4D::nx4 = nx4;
}
//////////////////////////////////////////////////////////////////////////
void SparseMatrix4D::getDimensions(size_t& nx1, size_t& nx2, size_t& nx3, size_t& nx4) 
{
   nx1 = SparseMatrix4D::nx1;
   nx2 = SparseMatrix4D::nx2;
   nx3 = SparseMatrix4D::nx3;
   nx4 = SparseMatrix4D::nx4;
}
/////////////////////////////////////////////////////////////////////////


#include "SparseMatrix3D.h"


size_t SparseMatrix3D::nx1 = 0;
size_t SparseMatrix3D::nx2 = 0;
size_t SparseMatrix3D::nx3 = 0;
//////////////////////////////////////////////////////////////////////////
SparseMatrix3D::SparseMatrix3D()
{

}
//////////////////////////////////////////////////////////////////////////
SparseMatrix3D::~SparseMatrix3D()
{

}
//////////////////////////////////////////////////////////////////////////
SparseMatrix3DPtr SparseMatrix3D::getInstance()
{
   if( !SparseMatrix3D::instance )
      SparseMatrix3D::instance = SparseMatrix3DPtr(new SparseMatrix3D());
   return SparseMatrix3D::instance;
}
//////////////////////////////////////////////////////////////////////////
void SparseMatrix3D::setDimensions(size_t nx1, size_t nx2, size_t nx3) 
{
   SparseMatrix3D::nx1 = nx1;
   SparseMatrix3D::nx2 = nx2;
   SparseMatrix3D::nx3 = nx3;
}
//////////////////////////////////////////////////////////////////////////
void SparseMatrix3D::getDimensions(size_t& nx1, size_t& nx2, size_t& nx3) 
{
   nx1 = SparseMatrix3D::nx1;
   nx2 = SparseMatrix3D::nx2;
   nx3 = SparseMatrix3D::nx3;
}
/////////////////////////////////////////////////////////////////////////


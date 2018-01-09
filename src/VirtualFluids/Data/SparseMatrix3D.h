#ifndef SparseMatrix3D_h
#define SparseMatrix3D_h

#include <boost/unordered_map.hpp>

#include <memory>
class SparseMatrix3D;
typedef std::shared_ptr<SparseMatrix3D> SparseMatrix3DPtr;

class SparseMatrix3D
{
public:
   static SparseMatrix3DPtr getInstance();
   static void setDimensions(size_t nx1, size_t nx2, size_t nx3);
   static void getDimensions(size_t& nx1, size_t& nx2, size_t& nx3);
   virtual ~SparseMatrix3D();
   //////////////////////////////////////////////////////////////////////////
   inline size_t index(size_t x1, size_t x2, size_t x3)
   {
      return  nx1 * ( nx2 * x3 + x2) + x1 ;
   }
protected:
private:
   SparseMatrix3D();
   static SparseMatrix3DPtr instance;
   static size_t nx1, nx2, nx3;
};

#endif

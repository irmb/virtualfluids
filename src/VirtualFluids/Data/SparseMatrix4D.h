#ifndef SparseMatrix4D_h
#define SparseMatrix4D_h

#include <boost/unordered_map.hpp>

#include <boost/shared_ptr.hpp>
class SparseMatrix4D;
typedef boost::shared_ptr<SparseMatrix4D> SparseMatrix4DPtr;

class SparseMatrix4D
{
public:
   static SparseMatrix4DPtr getInstance();
   static void setDimensions(size_t nx1, size_t nx2, size_t nx3, size_t x4);
   static void getDimensions(size_t& nx1, size_t& nx2, size_t& nx3, size_t& x4);
   virtual ~SparseMatrix4D();
   //////////////////////////////////////////////////////////////////////////
   inline size_t index(size_t x1, size_t x2, size_t x3, size_t x4)
   {
      return nx4*(nx3*(nx2*x1+x2)+x3)+x4;
   }
protected:
private:
   SparseMatrix4D();
   static SparseMatrix4DPtr instance;
   static size_t nx1, nx2, nx3, nx4;
};

#endif

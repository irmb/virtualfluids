
#ifndef MBSHAREDPOINTERDEFINES_H
#define MBSHAREDPOINTERDEFINES_H


// Boost includes
#include <boost/shared_ptr.hpp>
#include <boost/weak_ptr.hpp>
#include <boost/enable_shared_from_this.hpp>

#define VFSharedFromThis boost::enable_shared_from_this
#define VFSharedPtr boost::shared_ptr
#define VFWeakPtr   boost::weak_ptr
#define VFDynamicPtrCast boost::dynamic_pointer_cast

template<typename T>
class VFPtrDeleter
{
public:
   void operator()(T* p) { delete p; }
};



// std includes
#include <vector>

//#ifdef WIN32
//#  include <memory>
//#else
//#  include<tr1/memory>
//#endif

//#  define DCSharedFromThis std::tr1::enable_shared_from_this
//#  define DCSharedPtr std::tr1::shared_ptr
//#  define DCWeakPtr   std::tr1::weak_ptr
//#  define DCDynamicPtrCast std::tr1::dynamic_pointer_cast

#endif

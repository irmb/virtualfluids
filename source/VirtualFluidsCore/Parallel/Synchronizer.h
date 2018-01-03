#ifndef SYNCHRONIZER_H
#define SYNCHRONIZER_H

#include <memory>
#include <boost/thread/thread.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/thread/condition.hpp>


class Synchronizer;
typedef std::shared_ptr<Synchronizer> SynchronizerPtr;

class Synchronizer
{
public:
   Synchronizer()
      : threshold(0), count(0), generation(0)
   {}

   Synchronizer(int count)
      : threshold(count), count(count), generation(0)
   {}

   ~Synchronizer() {}

   void incNumberOfThreads()
   {
      threshold++;
      count++;
   }

   void setNumberOfThreads(int count)
   {
      this->threshold = count;
      this->count     = count;
   }

   bool wait()
   {
      boost::mutex::scoped_lock lock(mutex);
      unsigned int gen = generation;

      if (--count == 0)
      {
         generation++;
         count = threshold;
         cond.notify_all();
         return true;
      }

      while (gen == generation)
         cond.wait(lock);
      return false;
   }

   boost::mutex               gmtx;
private:
   boost::mutex              mutex;
   boost::condition_variable cond;

   unsigned int threshold;
   unsigned int count;
   unsigned int generation;
};


#endif 

#include <basics/utilities/UbLogger.h>

#if defined(CAB_BOOST) && !defined(NO_THREADSAFE_LOGGING)

boost::mutex Output2Stream::mtx;

#endif // CAB_BOOST


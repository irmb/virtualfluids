#ifdef VF_PYTHON

#include <boost/python.hpp>
#include "vfluids.h"
using namespace boost::python;

BOOST_PYTHON_MODULE(virtualfluids)
{
   //bool (Configuration::*get1)(const std::string&, std::string&)const = &Configuration::get;
   //bool (Configuration::*get2)(const std::string&, int&)const = &Configuration::get;
   //bool (Configuration::*get3)(const std::string&, long&)const = &Configuration::get;
   //bool (Configuration::*get4)(const std::string&, double&)const = &Configuration::get;
   //bool (Configuration::*get5)(const std::string&, bool&)const = &Configuration::get;
   //class_<Configuration>("Configuration")

   //   .def("load", &Configuration::load)
   //   .def("get", get1)
   //   .def("get", get2)
   //   .def("get", get3)
   //   .def("get", get4)
   //   .def("get", get5)
   //   ;
   class_<Configuration>("Configuration")

      .def("load", &Configuration::load)
      .def("getString", &Configuration::getString)
      .def("getInt", &Configuration::getInt)
      .def("getLong", &Configuration::getLong)
      .def("getDouble", &Configuration::getDouble)
      .def("getFloat", &Configuration::getFloat)
      .def("getBool", &Configuration::getBool)
      ;



}

#endif

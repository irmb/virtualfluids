#ifndef BASICS_CONFIGURATIONFILE_H
#define BASICS_CONFIGURATIONFILE_H

#include "Logger.h"
#include <filesystem>
#include <map>
#include <vector>
#include <sstream>
#include <string>
#include <fstream>
#include <iostream>
#include <stdlib.h>

#include <basics/utilities/UbException.h>

//! \brief  Simple configuration file
//! \details The Configuration class presented here can read and keep values of any configuration file written in a format like this:
//!#
//!# Simulation parameters
//!#
//!
//!nbDimensions    = 2
//!temperature     = 25.001
//!epsilon         = 1.013e-14
//!writeLogFile    = false      # NOTE: Set to "true" in debug mode only.
//!                             #       Logging slows down the program.
//!errorMessage    = the simulation failed
//!origin          = 0.0 0.0 0.0 # x, y, z of origin
//!
//!Example how to use it:
//!
//!ConfigurationFile   config;
//!config.load(configname);
//!
//!int            nbDimensions = config.getValue<int>("nbDimensions");
//!float          temperature  = config.getValue<float>("temperature");
//!double         epsilon      = config.getValue<double>("epsilon");
//!bool           writeLogFile = config.getValue<bool>("writeLogFile");
//!string         errorMessage = config.getValue<string>("errorMessage");
//!vector<double> origin       = config.getVector<double>("origin");
//!            
//! \author  Konstantin Kutscher


namespace vf::basics
{


class ConfigurationFile
{
public:
   //! clear all values
   void clear();

   //! load a configuration file
   bool load(const std::string& File);

   //! check if value associated with given key exists
   bool contains(const std::string& key) const;

   //! get vector with key
   template<class T>
   std::vector<T> getVector(const std::string& key) const;

   //! get value with key
   template<class T>
   T getValue(const std::string& key) const;

   //! get value with key and default value
   template<class T>
   T getValue(const std::string& key, T defaultValue) const;

   static ConfigurationFile loadConfig(int argc, char *argv[], std::filesystem::path configPath = "./config.txt")
   {
      // the config file's default path can be replaced by passing a command line argument

      if (argc > 1) 
      {
         configPath = argv[1];
         VF_LOG_INFO("Using command line argument for config path: {}", configPath.string());
      } else {
         VF_LOG_INFO("Using default config path: {}", configPath.string());
      }

      vf::basics::ConfigurationFile config;
      config.load(configPath.string());
      return config;
   }

private:
   //! the container
   std::map<std::string, std::string> data;

   //! get string with key
   std::string  getString(const std::string& key) const;

   //! remove leading and trailing tabs and spaces
   static std::string trim(const std::string& str);

   //! convert string to data type T
   template<class T>
   T fromString(const std::string& str) const;

   void split(std::vector<std::string>& lst, const std::string& input, const std::string& separators, bool remove_empty = true) const;
};


//////////////////////////////////////////////////////////////////////////
template<class T>
std::vector<T> ConfigurationFile::getVector(const std::string& key) const
{
   std::string str = getString(key);
   std::vector<T> v;
   std::vector<std::string> strings;
   split(strings, str, "\t\n\r;, ");
   for (std::vector<std::string>::iterator it = strings.begin(); it != strings.end(); ++it)
   {
      if (*it != "")
      {
         v.push_back(fromString<T>(*it));
      }
   }
   return v;
}
//////////////////////////////////////////////////////////////////////////
template<class T>
T ConfigurationFile::fromString(const std::string& str) const
{
   std::istringstream stream(str);
   T t;
   stream >> t;
   return t;
}

template<>
bool ConfigurationFile::fromString<bool>(const std::string& str) const;

//////////////////////////////////////////////////////////////////////////
template<class T>
T ConfigurationFile::getValue(const std::string& key) const
{
   std::string str = getString(key);
   bool bFlag = false;
   if ((std::string)typeid(T).name() == (std::string)typeid(bool).name()) 
   {
      bFlag = true;
   }

   std::istringstream iss(str);
   T x;
   iss >> x;
   if (!iss && !bFlag)
      UB_THROW(UbException(UB_EXARGS, " cannot convert \"" + str + "\" to type <" + static_cast<std::string>(typeid(x).name()) + ">"));

   if (bFlag)
   {
      bool value = (str == "true");
      x = value;
   }

   return x;
}

template<class T>
T ConfigurationFile::getValue(const std::string& key, T defaultValue) const
{
   if (contains(key))
   {
      return getValue<T>(key);
   }
   else
   {
      return defaultValue;
   }
}

}

#endif

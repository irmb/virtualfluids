//=======================================================================================
// ____          ____    __    ______     __________   __      __       __        __
// \    \       |    |  |  |  |   _   \  |___    ___| |  |    |  |     /  \      |  |
//  \    \      |    |  |  |  |  |_)   |     |  |     |  |    |  |    /    \     |  |
//   \    \     |    |  |  |  |   _   /      |  |     |  |    |  |   /  /\  \    |  |
//    \    \    |    |  |  |  |  | \  \      |  |     |   \__/   |  /  ____  \   |  |____
//     \    \   |    |  |__|  |__|  \__\     |__|      \________/  /__/    \__\  |_______|
//      \    \  |    |   ________________________________________________________________
//       \    \ |    |  |  ______________________________________________________________|
//        \    \|    |  |  |         __          __     __     __     ______      _______
//         \         |  |  |_____   |  |        |  |   |  |   |  |   |   _  \    /  _____)
//          \        |  |   _____|  |  |        |  |   |  |   |  |   |  | \  \   \_______
//           \       |  |  |        |  |_____   |   \_/   |   |  |   |  |_/  /    _____  |
//            \ _____|  |__|        |________|   \_______/    |__|   |______/    (_______/
//
//  This file is part of VirtualFluids. VirtualFluids is free software: you can
//  redistribute it and/or modify it under the terms of the GNU General Public
//  License as published by the Free Software Foundation, either version 3 of
//  the License, or (at your option) any later version.
//
//  VirtualFluids is distributed in the hope that it will be useful, but WITHOUT
//  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//  FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
//  for more details.
//
//  SPDX-License-Identifier: GPL-3.0-or-later
//  SPDX-FileCopyrightText: Copyright © VirtualFluids Project contributors, see AUTHORS.md in root folder
//
//! \addtogroup utilities
//! \ingroup basics
//! \{
//! \author Soeren Freudiger, Sebastian Geller
//=======================================================================================
#ifndef UBLOGGER_H
#define UBLOGGER_H

#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(_WIN64) || defined(__WIN64__)
#include <windows.h>
#else
#include <sys/time.h>
#endif

#if defined(CAB_BOOST) && !defined(NO_THREADSAFE_LOGGING)
#include <boost/thread.hpp>
#endif // CAB_BOOST

enum LogLevel { logERROR, logWARNING, logINFO, logDEBUG, logDEBUG1, logDEBUG2, logDEBUG3, logDEBUG4, logDEBUG5 };

//////////////////////////////////////////////////////////////////////////
// template <typename OutputPolicy> class Log  - declaration
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//!
//! \brief
//! C++ Logger
//! \details Functionality:
//! Per logentry an object of type UbLogger is generated, the log string is passed to this object and
//! upon destruction of the object the string is written to a file or the screen depending on
//! the policy (=template paramter). Multiple log level are supported.
//!
//! helpermakro:  UBLOG
//!
//! Example 1:
//! \code
//! UBLOG(logINFO) << "Klasse::foo entered"; //endl is not required
//! \endcode
//!
//! Example 2:
//! \code
//! try
//! {
//!    UbLog::reportingLevel() = UbLog::logLevelFromString("DEBUG3");
//!    //UbLog::output_policy::setStream(&std::cerr); //<- clog is stdandard
//!    UbLog::output_policy::setStream("c:/temp/out.txt");  //you can not open these -> error message -> log is output
//!    in cerr
//!
//!    int count = 3;
//!    UBLOG(logINFO, "A loop with " << count << " iterations");
//!    for (int i = 0; i != count; ++i)
//!    {
//!        UBLOG(logERROR , "error  - the counter i = " << i );
//!        UBLOG(logDEBUG1, "debug1 - the counter i = " << i );
//!        UBLOG(logDEBUG2, "debug2 - the counter i = " << i );
//!        UBLOG(logDEBUG3, "debug3 - the counter i = " << i );
//!        //for MultiLine entries: -> formatting in logfile
//!        UBLOGML(logDEBUG3, "debug3 - the counter i = "<<endl<<" 2 zeile "<< "3. Zeile" << i);
//!        UBLOGML(logDEBUG3, "debug3 - the counter i = "<<endl<<" 2 zeile "<< "3. Zeile" << i);
//!        UBLOG2ML(logDEBUG3,std:cout,"debug3 - the counter i = "<<endl<<" 2 zeile "<< "3. Zeile" << i);
//!    }
//!    return 0;
//! }
//! catch(const std::exception& e)
//! {
//!    UBLOG(logERROR) << e.what();
//! }
//! \endcode
//! Idee based on:
//! Paper by Dr. Dobbs Portal,
//! September 05, 2007,
//! Logging In C++
//!
//////////////////////////////////////////////////////////////////////////
template <typename OutputPolicy>
class UbLogger
{
public:
    using output_policy = OutputPolicy;

public:
    UbLogger();
    virtual ~UbLogger();
    std::ostringstream &get(const LogLevel &level = logINFO);

public:
    // static, weil man so spaeter die ObjErstellunge ersparen kann,
    // falls level kleiner als Level
    static LogLevel &reportingLevel();

    static std::string logLevelToString(const LogLevel &level);
    static LogLevel logLevelFromString(const std::string &level);

    static std::string logTimeString();

protected:
    std::ostringstream os;

private:
    UbLogger(const UbLogger &);
    UbLogger &operator=(const UbLogger &);
};

//////////////////////////////////////////////////////////////////////////
// template <typename OutputPolicy> class Log  - implementation
//////////////////////////////////////////////////////////////////////////
template <typename OutputPolicy>
UbLogger<OutputPolicy>::UbLogger() = default;
/*==========================================================*/
template <typename OutputPolicy>
std::ostringstream &UbLogger<OutputPolicy>::get(const LogLevel &level)
{
    os << logTimeString() << " " << std::setw(6)
#if defined(CAB_BOOST) && !defined(NO_MT_LOGGING)
       << boost::this_thread::get_id() << " "
#endif
       << std::setw(8) << std::left << UbLogger<OutputPolicy>::logLevelToString(level) << ": "
       << std::string(level > logDEBUG ? 3 * (level - logDEBUG) : 0, ' '); //<baumartiger output :D

    return os;
}
/*==========================================================*/
template <typename OutputPolicy>
UbLogger<OutputPolicy>::~UbLogger()
{
    os << std::endl;
    OutputPolicy::output(os.str());
}
/*==========================================================*/
template <typename OutputPolicy>
LogLevel &UbLogger<OutputPolicy>::reportingLevel()
{
    static LogLevel reportLevel = logINFO;
    return reportLevel;
}
/*==========================================================*/
template <typename OutputPolicy>
std::string UbLogger<OutputPolicy>::logLevelToString(const LogLevel &level)
{
    static std::string const buffer[] = { "ERROR",  "WARNING", "INFO",   "DEBUG", "DEBUG1",
                                          "DEBUG2", "DEBUG3",  "DEBUG4", "DEBUG5" };
    return buffer[level];
}
/*==========================================================*/
template <typename OutputPolicy>
LogLevel UbLogger<OutputPolicy>::logLevelFromString(const std::string &level)
{
    if (level == "DEBUG5")
        return logDEBUG5;
    if (level == "DEBUG4")
        return logDEBUG4;
    if (level == "DEBUG3")
        return logDEBUG3;
    if (level == "DEBUG2")
        return logDEBUG2;
    if (level == "DEBUG1")
        return logDEBUG1;
    if (level == "DEBUG")
        return logDEBUG;
    if (level == "INFO")
        return logINFO;
    if (level == "WARNING")
        return logWARNING;
    if (level == "ERROR")
        return logERROR;

    UbLogger<OutputPolicy>().get(logWARNING)
        << "UbLogger<OutputPolicy>::logLevelFromString(level) - unknown logging level '" << level
        << "'. Using INFO level as default.";
    return logINFO;
}

//////////////////////////////////////////////////////////////////////////
// logTimeString
//////////////////////////////////////////////////////////////////////////
#if defined(WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(_WIN64) || defined(__WIN64__)
template <typename OutputPolicy>
inline std::string UbLogger<OutputPolicy>::logTimeString()
{
    const int MAX_LEN = 200;
    char buffer[MAX_LEN];
    if (GetTimeFormatA(LOCALE_USER_DEFAULT, 0, 0, "HH':'mm':'ss", buffer, MAX_LEN) == 0) {
        return "Error in std::string UbLogger<OutputPolicy>::logTimeString()";
    }

    char result[100]   = { 0 };
    static DWORD first = (DWORD)GetTickCount64();
    std::sprintf(result, "%s.%03ld", buffer, (long)(GetTickCount64() - first) % 1000);
    return result;
}
#else
template <typename OutputPolicy>
inline std::string UbLogger<OutputPolicy>::logTimeString()
{
    char buffer[11];
    time_t t;
    time(&t);
    tm r; // = { 0 };
    strftime(buffer, sizeof(buffer), "%X", localtime_r(&t, &r));
    struct timeval tv;
    gettimeofday(&tv, 0);
    char result[100] = { 0 };
    std::sprintf(result, "%s.%03ld", buffer, (long)tv.tv_usec / 1000);
    return result;
}
#endif

//////////////////////////////////////////////////////////////////////////
//! Implementation of OutputPolicy)
//////////////////////////////////////////////////////////////////////////
class Output2Stream // implementation of OutputPolicy
{
public:
    static std::ostream *&getStream();
    static void output(const std::string &msg);

    //! creates output-file-stream (of file opening fails -> stream is set to std::cerr)
    static void setStream(const std::string &filename);

    //! direct set outputstream, gcControl = true -> object will be deleted by Output2Stream
    static void setStream(std::ostream *pStream, const bool &gcControl = false);

protected:
#if defined(CAB_BOOST) && !defined(NO_MT_LOGGING)
    static boost::mutex mtx;
#endif
};
/*==========================================================*/
inline std::ostream *&Output2Stream::getStream()
{
    static std::ostream *pStream = &std::clog;
    return pStream;
}
/*==========================================================*/
inline void Output2Stream::setStream(std::ostream *pFile, const bool &gcControl)
{
#if defined(CAB_BOOST) && !defined(NO_MT_LOGGING)
    boost::mutex::scoped_lock lock(mtx);
#endif
    static bool s_gcControl = false;

    if (s_gcControl && Output2Stream::getStream()) {
        delete Output2Stream::getStream();
    }

    s_gcControl = gcControl;

    Output2Stream::getStream() = pFile;
}
/*==========================================================*/
inline void Output2Stream::setStream(const std::string &filename)
{
    std::ofstream *file = new std::ofstream(filename.c_str());
    if (!(*file)) {
        delete file;
        Output2Stream::setStream(&std::cerr, false);
        UbLogger<Output2Stream>().get(logERROR)
            << " Output2Stream::setStream(const std::string& filename) could not open file " << filename
            << " -> std::cerr is used instead " << std::endl;
        return;
    }
    std::cout << "UbLog writes to " << filename << std::endl;
    Output2Stream::setStream(file, true);
}
/*==========================================================*/
inline void Output2Stream::output(const std::string &msg)
{
#if defined(CAB_BOOST) && !defined(NO_MT_LOGGING)
    boost::mutex::scoped_lock lock(mtx);
#endif
    std::ostream *pStream = getStream();
    if (!pStream)
        return;
    (*pStream) << msg << std::flush;
}

//////////////////////////////////////////////////////////////////////////
// UbLog
//////////////////////////////////////////////////////////////////////////
class UbLog : public UbLogger<Output2Stream>
{
};

// Macro to limit compiler-side maxLevel
#ifndef UBLOG_MAX_LEVEL
#define UBLOG_MAX_LEVEL logDEBUG5
#endif

//////////////////////////////////////////////////////////////////////////
// example UBLOG(logINFO) << "das ist ein log eintrag";
//////////////////////////////////////////////////////////////////////////
#define UBLOG(level, logtext)                                                                                          \
    if (level <= UBLOG_MAX_LEVEL && level <= UbLog::reportingLevel() && Output2Stream::getStream()) {                  \
        UbLog().get(level) << logtext;                                                                                 \
    }

//////////////////////////////////////////////////////////////////////////
// makro 2 fuer korrekten MultiLineOutput (teuer!!)
// example1: UBLOGML(logINFO, "line1"<<endl<<"line2"<<endl<<"line3" )
// example2: UBLOGML(logINFO, "line1\nline2\nendl\nline3" )
//////////////////////////////////////////////////////////////////////////
#define UBLOGML(level, multiline)                                                                                      \
    if (level > UBLOG_MAX_LEVEL || level > UbLog::reportingLevel() || !Output2Stream::getStream())                     \
        ;                                                                                                              \
    else {                                                                                                             \
        std::ostringstream output;                                                                                     \
        output << multiline;                                                                                           \
        std::istringstream input(output.str());                                                                        \
        while (!input.eof()) {                                                                                         \
            std::string dummy;                                                                                         \
            getline(input, dummy, '\n');                                                                               \
            UbLog().get(level) << dummy;                                                                               \
        }                                                                                                              \
    }
//////////////////////////////////////////////////////////////////////////
// makro3, falls auch bildschirmausgabe erwuenscht
//   -> es wird sowohl ins logfile als auch auf den "stream" geschrieben
//      wenn reporting level und level passen :D
// example1: UBLOG2ML(logINFO, std::cout,  "line1"<<endl<<"line2"<<endl<<"line3" )
// example2: UBLOG2ML(logINFO, std::cout,  "line1\nline2\nendl\nline3" )
//////////////////////////////////////////////////////////////////////////
#define UBLOG2(level, stream, text)                                                                                    \
    if (level > UBLOG_MAX_LEVEL || level > UbLog::reportingLevel() || !Output2Stream::getStream())                     \
        ;                                                                                                              \
    else {                                                                                                             \
        stream << text << std::endl;                                                                                   \
        UbLog().get(level) << text;                                                                                    \
    }

//////////////////////////////////////////////////////////////////////////
// makro4, wie 3 nur mit multiline
// example: UBLOG2(logINFO, std::cout,  "test" )
//////////////////////////////////////////////////////////////////////////
#define UBLOG2ML(level, stream, multiline)                                                                             \
    if (level > UBLOG_MAX_LEVEL || level > UbLog::reportingLevel() || !Output2Stream::getStream())                     \
        ;                                                                                                              \
    else {                                                                                                             \
        stream << multiline << std::endl;                                                                              \
        std::ostringstream output;                                                                                     \
        output << multiline;                                                                                           \
        std::istringstream input(output.str());                                                                        \
        while (!input.eof()) {                                                                                         \
            std::string dummy;                                                                                         \
            getline(input, dummy, '\n');                                                                               \
            UbLog().get(level) << dummy;                                                                               \
        }                                                                                                              \
    }

//////////////////////////////////////////////////////////////////////////
// example 2
//////////////////////////////////////////////////////////////////////////
// try
// {
//    UbLog::reportingLevel() = UbLog::logLevelFromString("DEBUG3");
//    //UbLog::output_policy::setStream(&std::cerr); //<- clog ist stdandard
//    UbLog::output_policy::setStream("c:/temp/out.txt");  //kann man diese nicht oeffnen -> fehlermeldung -> Log wird
//    in cerr ausgegben
//
//    int count = 3;
//    UBLOG(logINFO, "A loop with " << count << " iterations");
//    for (int i = 0; i != count; ++i)
//    {
//        UBLOG(logERROR , "error  - the counter i = " << i );
//        UBLOG(logDEBUG1, "debug1 - the counter i = " << i );
//        UBLOG(logDEBUG2, "debug2 - the counter i = " << i );
//        UBLOG(logDEBUG3, "debug3 - the counter i = " << i );
//        //fuer MultiLine Eintraege: --> koerrekte formatierung im logfile
//        UBLOGML(logDEBUG3, "debug3 - the counter i = "<<endl<<" 2 zeile "<< "3. Zeile" << i);
//        UBLOGML(logDEBUG3, "debug3 - the counter i = "<<endl<<" 2 zeile "<< "3. Zeile" << i);
//        UBLOG2ML(logDEBUG3,std:cout,"debug3 - the counter i = "<<endl<<" 2 zeile "<< "3. Zeile" << i);
//    }
//    return 0;
// }
// catch(const std::exception& e)
// {
//    UBLOG(logERROR) << e.what();
// }

#endif // UBLOGGER_H

//! \}

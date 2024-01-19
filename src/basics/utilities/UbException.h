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
#ifndef UBEXCEPTION_H
#define UBEXCEPTION_H

#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "./UbTuple.h"

//=========================================================================
//
//! \brief UbException
//! usage: UB_THROW( UbException("error message") );
//!        UB_THROW( UbException(__FILE__, __LINE__,"error message") );
//!        UB_THROW( UbException(__FILE__, __LINE__,UB_FUNCTION,"error message") );
//!        UB_THROW( UbException(UB_EXARGS,"error") ); //same as above
//
//=========================================================================

// Macro UB_FUNCTION: figures out the method/function name (platform dependant)
#if defined(__GNUC__) || (defined(__MWERKS__) && (__MWERKS__ >= 0x3000)) || (defined(__ICC) && (__ICC >= 600))
#define UB_FUNCTION __PRETTY_FUNCTION__
#elif defined(__DMC__) && (__DMC__ >= 0x810)
#define UB_FUNCTION __PRETTY_FUNCTION__
#elif defined(__FUNCSIG__)
#define UB_FUNCTION __FUNCSIG__
#elif (defined(__INTEL_COMPILER) && (__INTEL_COMPILER >= 600)) || (defined(__IBMCPP__) && (__IBMCPP__ >= 500))
#define UB_FUNCTION __FUNCTION__
#elif defined(__BORLANDC__) && (__BORLANDC__ >= 0x550)
#define UB_FUNCTION __FUNC__
#elif defined(__STDC_VERSION__) && (__STDC_VERSION__ >= 199901)
#define UB_FUNCTION __func__
#else
#define UB_FUNCTION "(unknown)"
#endif

// Helper Marco
#define UB_EXARGS __FILE__, __LINE__, UB_FUNCTION

#ifdef CAB_BOOST
#define UB_THROW(e) throw boost::enable_current_exception(e)
#else
#define UB_THROW(e) throw e
#endif

class UbException : public std::runtime_error
{
public:
    using ExceptionData = UbTuple<std::string, int, std::string, std::string>;

public:
    //////////////////////////////////////////////////////////////////////////
    // constructors
    UbException() : std::runtime_error("") {}
    /*==========================================================*/
    UbException(const std::string &str) : std::runtime_error("") { this->addInfo(str); }
    /*==========================================================*/
    UbException(const std::string &file, const int &line, const std::string &err_str) : std::runtime_error("")
    {
        this->addInfo(file, line, "unknown", err_str);
    }
    /*==========================================================*/
    // UbException(const char* file, const int& line, const char* function, const std::string& err_str)
    UbException(const std::string &file, const int &line, const std::string &function, const std::string &err_str)
        : std::runtime_error("")
    {
        this->addInfo(file, line, function, err_str);
    }
    //////////////////////////////////////////////////////////////////////////
    // destructor
    ~UbException() noexcept override = default;
    //////////////////////////////////////////////////////////////////////////
    // virtual public methods
    // returns  exception-string
    const char *what() const noexcept override
    {
        exceptionString = this->toString();
        return exceptionString.c_str(); // ansonsten ist das Verhalten anschliessend undefiniert!
    }
    /*==========================================================*/
    virtual void addInfo(const std::string &err_str)
    {
        exceptionData.push_back(makeUbTuple(std::string("-"), 0, std::string("unknown"), err_str));
    }
    /*==========================================================*/
    // add exception
    virtual void addInfo(const std::string &file, const int &line, const std::string &function,
                         const std::string &err_str)
    {
        exceptionData.push_back(makeUbTuple(file, line, function, err_str));
    }
    /*==========================================================*/
    // returns exception-string with all calles exceptions
    virtual const std::vector<std::string> getInfo() const
    {
        std::vector<std::string> tmp;
        for (std::size_t i = 0; i < exceptionData.size(); i++) {
            std::stringstream str;
            str << val<1>(exceptionData[i]) << ", " << val<2>(exceptionData[i]) << ", " << val<3>(exceptionData[i])
                << ", " << val<4>(exceptionData[i]);
            tmp.push_back(str.str());
        }
        return tmp;
    }
    /*==========================================================*/
    // returns exception-string with all calles exceptions and detailes informations
    virtual std::string toString() const
    {
        std::stringstream str("UbExeption");

        for (std::size_t i = 0; i < exceptionData.size(); i++)
            str << (std::string) "caller[" << i << "]\n"
                << "  - file:     " << val<1>(exceptionData[i]) << "\n"
                << "  - line:     " << val<2>(exceptionData[i]) << "\n"
                << "  - function: " << val<3>(exceptionData[i]) << "\n"
                << "  - what:     " << val<4>(exceptionData[i]) << std::endl;

        return str.str();
    }

protected:
    //////////////////////////////////////////////////////////////////////////
    // protected member
    std::vector<ExceptionData> exceptionData;
    mutable std::string exceptionString;
};

// overlading operator <<
inline std::ostream &operator<<(std::ostream &os, const UbException &e) { return os << e.toString(); }

#endif // UBEXCEPTION_H

//! \}

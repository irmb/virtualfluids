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
//  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
//  for more details.
//
//  You should have received a copy of the GNU General Public License along
//  with VirtualFluids (see COPYING.txt). If not, see <http://www.gnu.org/licenses/>.
//
//! \author Martin Schoenherr
//=======================================================================================
#ifndef RestartObject_H
#define RestartObject_H

#include <memory>
#include <string>
#include <vector>

#include <basics/DataTypes.h>

class Parameter;

class RestartObject
{
public:
    virtual ~RestartObject() = default;

    void deserialize(const std::string &filename, std::shared_ptr<Parameter>& para);
    void serialize(const std::string &filename, const std::shared_ptr<Parameter>& para);

    std::vector<std::vector<real>> fs;

    virtual void serialize_internal(const std::string &filename)   = 0;
    virtual void deserialize_internal(const std::string &filename) = 0;

private:
    void clear(const std::shared_ptr<Parameter>& para);
};


class ASCIIRestartObject : public RestartObject
{
private:
    virtual void serialize_internal(const std::string &filename);
    virtual void deserialize_internal(const std::string &filename);
};

class BinaryRestartObject : public RestartObject
{
private:
    virtual void serialize_internal(const std::string &filename);
    virtual void deserialize_internal(const std::string &filename);
};

#endif

#ifdef BUILD_JSONCPP

#include "JsonInput.h"

#include <fstream>
#include <iterator>
#include <sstream>
#include <string>
#include <vector>

namespace input
{
template <typename Out>
void split(const std::string &s, char delim, Out result)
{
    std::stringstream ss;
    ss.str(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        *(result++) = item;
    }
}

std::vector<std::string> split(const std::string &s, char delim)
{
    std::vector<std::string> elems;
    split(s, delim, std::back_inserter(elems));
    return elems;
}

JsonInput::JsonInput(std::istream &stream)
{
    Json::Reader reader;
    reader.parse(stream, jsonValue);
}

bool JsonInput::hasValue(const std::string &key) const
{
    auto keys = split(key, ' ');

    if (keys.size() == 1 && !jsonValue[keys[0]].isNull())
        return true;
    else if (keys.size() == 2 && !jsonValue[keys[0]][keys[1]].isNull())
        return true;
    else if (keys.size() == 3 && !jsonValue[keys[0]][keys[1]][keys[2]].isNull())
        return true;
    else
        return false;
}

std::string JsonInput::getValue(const std::string &key)
{
    auto keys = split(key, ' ');

    if (keys.size() == 1)
        return jsonValue[keys[0]].asString();
    else if (keys.size() == 2)
        return jsonValue[keys[0]][keys[1]].asString();
    else if (keys.size() == 3)
        return jsonValue[keys[0]][keys[1]][keys[2]].asString();
    else
        return "";
}

} // namespace input

#endif

#ifdef BUILD_JSONCPP

#ifndef JsonInput_H
#define JsonInput_H


#include <string>
#include <json/json.h>


#include "basics_export.h"


#include "../Input.h"

namespace input
{
	class JsonInput : public Input
	{
	public:
        BASICS_EXPORT JsonInput(std::istream &stream);

        BASICS_EXPORT virtual bool hasValue(const std::string &key) const override;
        BASICS_EXPORT virtual std::string getValue(const std::string &key) override;

    private:
        Json::Value jsonValue;

    };
}

#endif

#endif

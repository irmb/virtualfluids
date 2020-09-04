#ifndef Input_H
#define Input_H

#include "basics_export.h"


#include <string>
#include <memory>
#include <istream>

namespace input
{
	class Input
	{
	public:
        static BASICS_EXPORT std::unique_ptr<Input> makeInput(std::istream &stream, const std::string &inputType);

        virtual ~Input() = default;

        virtual bool hasValue(const std::string &key) const = 0;
        virtual std::string getValue(const std::string &key) = 0;
	};
}

#endif

#ifndef COLOR_OUTPUT_H
#define COLOR_OUTPUT_H

#include <gtest/gtest.h>

// https://stackoverflow.com/questions/16491675/how-to-send-custom-message-in-google-c-testing-framework/29155677#29155677
namespace testing
{
	namespace internal
	{
		enum GTestColor {
			COLOR_DEFAULT,
			COLOR_RED,
			COLOR_GREEN,
			COLOR_YELLOW
		};

		// in case of unresoved external while using shared libraries
		// add in gtest.h line 167 and 168:
		// enum GTestColor;
		//       void GTEST_API_ ColoredPrintf(GTestColor color, const char* fmt, ...);
		// see commit: 4c0ed885ceab18b9df7a2495c77a51e236aee6f1
		extern void ColoredPrintf(GTestColor color, const char* fmt, ...);
	}
}
#endif
#ifndef TEST_COUT_H
#define TEST_COUT_H


#include <sstream>
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
#define PRINTTESTINFOSTART()  do { testing::internal::ColoredPrintf(testing::internal::COLOR_GREEN, "[TestInfo  ]");std::cout << std::endl;} while(0)
#define PRINTF(h)  do { testing::internal::ColoredPrintf(testing::internal::COLOR_GREEN, "[          ] "); testing::internal::ColoredPrintf(testing::internal::COLOR_DEFAULT, h); std::cout << std::endl;} while(0)
#define PRINTTESTINFOEND()  do { testing::internal::ColoredPrintf(testing::internal::COLOR_GREEN, "[  TestInfo]");std::cout << std::endl;} while(0)

#define PRINTRAUTE()  do { testing::internal::ColoredPrintf(testing::internal::COLOR_GREEN, "#################################################");std::cout << std::endl;} while(0)
#define PRINTC(h)  do { testing::internal::ColoredPrintf(testing::internal::COLOR_GREEN, h); std::cout << std::endl;} while(0)

// C++ stream interface
class TestCout : public std::stringstream
{
public:
	TestCout(std::string testName1, int l1, int l2, std::string nameWerte1, std::string nameWerte2, std::string nameWerte3, double testWert1, double testWert2, double testWert3);
	TestCout(std::string testName, int l);
	~TestCout();

private:
};

#define TEST_COUT(testName1, l1, l2, nameWerte1, nameWerte2, nameWerte3, testWert1, testWert2, testWert3)  TestCout(testName1, l1, l2, nameWerte1, nameWerte2, nameWerte3,testWert1, testWert2, testWert3)
#define TEST_HEAD(testName, l) TestCout(testName, l)

#endif 
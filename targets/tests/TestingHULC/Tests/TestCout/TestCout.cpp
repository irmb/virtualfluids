#include "TestCout.h"

#include <iomanip>
#include <ctime>

TestCout::TestCout(std::string testName1, int l1, int l2, std::string nameWerte1, std::string nameWerte2, std::string nameWerte3, double testWert1, double testWert2, double testWert3)
{
	std::ostringstream oss1;
	oss1 << testName1;
	std::string output1 = oss1.str();

	std::ostringstream oss2;
	oss2 << "L: " << l1 << "\t" << "\t" << "\t" << "L: " << l2;
	std::string output2 = oss2.str();

	std::ostringstream oss3;
	oss3 << nameWerte1 << ": " << testWert1 << std::setw(5) << "\t" << nameWerte2 << ": " << testWert2;
	std::string output3 = oss3.str();

	std::ostringstream oss4;
	oss4 << nameWerte3 << ": " << testWert3;
	std::string output4 = oss4.str();

	PRINTTESTINFOSTART();
	PRINTF(output1.c_str());
	PRINTF(output2.c_str());
	PRINTF(output3.c_str());
	PRINTF(output4.c_str());
	PRINTTESTINFOEND();
}

TestCout::TestCout(std::string testName, int l)
{
	std::ostringstream oss1;
	oss1 << "# TEST: " << std::setfill(' ') << std::left << std::setw(40) << testName << "#";
	std::string output1 = oss1.str();

	std::ostringstream oss2;
	oss2 << "# L: " << std::setfill(' ') << std::left << std::setw(43) << l << "#";
	std::string output2 = oss2.str();

	std::ostringstream oss3;
	time_t now;
	struct tm nowLocal;
	now = time(NULL);
	nowLocal = *localtime(&now);
	oss3 << "# DATUM: " << std::setfill('0') << std::setw(2) << nowLocal.tm_mday << "." << std::setw(2) << nowLocal.tm_mon + 1 << "." << nowLocal.tm_year + 1900 << "   TIME: " << std::setw(2) << nowLocal.tm_hour << ":" << std::setw(2) << nowLocal.tm_min << ":" << std::setw(2) << nowLocal.tm_sec << "\t" << "\t#";
	std::string output3 = oss3.str();

	PRINTRAUTE();
	PRINTC(output1.c_str());
	PRINTC(output2.c_str());
	PRINTC(output3.c_str());
	PRINTRAUTE();
}

TestCout::~TestCout()
{
	
}

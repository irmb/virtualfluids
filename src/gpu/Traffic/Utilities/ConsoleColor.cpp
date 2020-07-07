#include "ConsoleColor.h"



#include <VirtualFluidsDefinitions.h>
#include "Core/DataTypes.h"


//// Windows //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef WIN32

#include <iostream>
#include <iomanip>	//formatting output streams
#include <windows.h> //for colourful console output


void ConsoleColor::setDefaultWhite()
{
	SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE), 7); // set output default white 7;
}


void ConsoleColor::setDarkGrey()
{
	SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE), 8); //set output dark grey 8, dark blue 1, black 0;
}


void ConsoleColor::setBrightRed()
{
	SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE), 12); //set output bright green 10, bright red 12;
}


void ConsoleColor::setBlack()
{
	SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE), 0); //set output dark grey 8, dark blue 1, black 0;
}


void ConsoleColor::setBrightGreen()
{
	SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE), 10); //set output bright green 10, bright red 12;
}


#endif WIN32





//// Linux //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef WIN32

void ConsoleColor::setDefaultWhite()
{}

void ConsoleColor::setDarkGrey()
{}

void ConsoleColor::setBrightRed()
{}

void ConsoleColor::setBlack()
{}

void ConsoleColor::setBrightGreen()
{}

#endif // !WIN32



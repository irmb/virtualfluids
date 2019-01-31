#include "Utilities\LogFileData\LogFileData.h"
#include "Utilities\LogFileReader\LogFileReader.h"

#include <mpi.h>
#include <memory>
#include <gmock/gmock.h>
#include <cstdio>
#include "wstp.h"

#include <iostream>
#include <fstream>

void waitForLinkActivity(WSLINK link)
{
	switch (WSWaitForLinkActivity(link)){
		case WSWAITERROR:
			fprintf(stderr, "Something went wrong when waiting for link");
			break;
		case WSWAITSUCCESS:
			fprintf(stderr, "\nwait succes");
	}
}

WSENV initWSTP() {
	WSENV env = WSInitialize((WSEnvironmentParameter)0);
	return env;
}

WSLinkServer initMathematicaServer(WSENV env, int &error) {
	unsigned short port = 8000;
	WSLinkServer server = WSNewLinkServerWithPort(env, port, NULL, &error);
	if (error != WSEOK) {
		fprintf(stderr, "Something went wrong when starting up the server");
	}

	return server;
}

WSLINK initMathematicaLink(WSLinkServer server, int &error){
	const char *interface = WSInterfaceFromLinkServer(server, &error);
	unsigned short port = WSPortFromLinkServer(server, &error);

	std::ofstream myMathematicaFile("C:\\Users\\Timon\\Desktop\\myFile1.txt");
	myMathematicaFile << "link = LinkConnect[\"" << port << "\", LinkProtocol -> \"TCPIP\", LinkHost -> \"" << interface << "\", LinkOptions -> 4];" << std::endl;
	myMathematicaFile << "If[LinkReadyQ[link], LinkRead[link]]";
	myMathematicaFile.close();

	std::rename("C:\\Users\\Timon\\Desktop\\myFile1.txt", "C:\\Users\\Timon\\Desktop\\myFile1.nb");
	system("C:\\Users\\Timon\\Desktop\\myFile1.nb");

	fprintf(stderr, "link = LinkConnect[\"%d\", LinkProtocol -> \"TCPIP\", LinkHost -> \"%s\", LinkOptions -> 4];\nIf[LinkReadyQ[link], LinkRead[link]]", port, interface);

	WSLINK link = WSWaitForNewLinkFromLinkServer(server, &error);
	if (link == (WSLINK)0 || error != WSEOK)
		std::cout << "unable to open link" << std::endl << std::flush;
	else
		std::cout << "\nConected to new link...";

	WSActivate(link);
	WSPutFunction(link, "Print", 1);
		WSPutString(link, "Hello client program.\nRead instructions in programm output.");
	WSEndPacket(link);
	WSFlush(link);

	std::cout << std::endl << std::endl << "Copy to Mathematica to read Data:" << std::endl << std::endl;
	std::cout << "Do[If[LinkReadyQ[link], Print[LinkRead[link]], LinkWrite[link, \"End\"]; Break[]], 20]" << std::endl;

	return link;
}

void DeinitializeMathematica(WSENV env, WSLinkServer server, WSLINK link) {
	WSClose(link);
	WSShutdownLinkServer(server);
	WSDeinitialize(env);
}

void variante1() {
	int error;
	WSENV env = initWSTP();
	WSLinkServer server = initMathematicaServer(env, error);
	WSLINK link = initMathematicaLink(server, error);

	WSPutFunction(link, "Plot", 3);
	WSPutSymbol(link, "x");
	WSPutFunction(link, "List", 3);
	WSPutSymbol(link, "x");
	WSPutInteger(link, 0);
	WSPutInteger(link, 1);
	WSPutFunction(link, "Rule", 2);
	WSPutSymbol(link, "AxesLabel");
	WSPutFunction(link, "List", 2);
	WSPutSymbol(link, "x");
	WSPutSymbol(link, "y");
	WSEndPacket(link);
	WSFlush(link);

	WSPutFunction(link, "Plot", 2);
	WSPutFunction(link, "Times", 2);
	WSPutSymbol(link, "x");
	WSPutSymbol(link, "x");
	WSPutFunction(link, "List", 3);
	WSPutSymbol(link, "x");
	WSPutInteger(link, 0);
	WSPutInteger(link, 2);
	WSEndPacket(link);
	WSFlush(link);

	int number = 2;

	WSPutFunction(link, "Plot", number);
	WSPutSymbol(link, "x");
	WSPutFunction(link, "List", 3);
	WSPutSymbol(link, "x");
	WSPutInteger(link, 0);
	WSPutInteger(link, 3);
	WSEndPacket(link);
	WSFlush(link);

	double list[20];

	for (int i = 0; i < 20; i++)
		list[i] = i + .9;

	WSPutFunction(link, "Set", 2);
	WSPutSymbol(link, "list");
	if (!WSPutReal64List(link, (double *)list, 20))
		std::cout << "unable to put the list to link" << std::endl << std::flush;


	std::cout << std::endl << "Waiting for Mathematica to read Data";
	waitForLinkActivity(link);
	DeinitializeMathematica(env, server, link);
	std::cout << "\nend" << std::endl;
}

void variante2() {
	std::ofstream myMathematicaFile("C:\\Users\\Timon\\Desktop\\myFile2.txt");

	myMathematicaFile << "Plot[x^2, {x, 0, 1}, AxesLabel -> {x, y}]" << std::endl;
	myMathematicaFile << "list = {0, 1, 2, 3}" << std::endl;
	myMathematicaFile << "Plot[list, {x, 0, 1}]" << std::endl;

	myMathematicaFile.close();

	std::rename("C:\\Users\\Timon\\Desktop\\myFile2.txt", "C:\\Users\\Timon\\Desktop\\myFile2.nb");
	system("C:\\Users\\Timon\\Desktop\\myFile2.nb");
}


int main(int argc, char **argv)
{
	std::shared_ptr<LogFileReader> logFileReader = LogFileReader::getInstance();

	std::shared_ptr<LogFileData> logFileData = logFileReader->readLogFileToLogFileData("C:\\Users\\Timon\\Documents\\studienarbeitIRMB\\logFiles\\NumericalTestLogFiles\\TaylorGreenVortexU0\\viscosity_0.0001\\u0_ 0.032_Amplitude_ 0.01\\CumulantAA2016CompSP27\\logfile_20181212_164220_CumulantAA2016CompSP27_vis_0.0001.txt");

	

	//variante1();
	variante2();
	
	return 0;
}

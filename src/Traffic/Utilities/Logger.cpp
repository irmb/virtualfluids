#include "Logger.h"

#include <stdexcept>
#include <iostream>
#include <time.h>

TrafficLogger TrafficLogger::instance;

void TrafficLogger::startLogger(std::string filename)
{
	instance.filename = filename;

	instance.file.open(filename.c_str(), std::ios::app);
	try { if (instance.file.fail()) throw std::runtime_error("couldn't open file for logger: " + filename); }
	catch (const std::exception& e) {
		std::cerr << e.what() << std::endl;
		std::cin.get();
		exit(EXIT_FAILURE);
	}
}


void TrafficLogger::writeSimulationStart(std::string info, bool useGPU)
{
	time_t now = time(0);

	instance.file << "Simulation started at: " << ctime(&now);
	instance.file << "Info: " << info << "\t \t" << "simulating on the ";
	if (useGPU) instance.file << "GPU \t\t";
	else instance.file << "CPU \t\t";
	#ifdef NDEBUG
		instance.file << "Release \n";
	#else
		instance.file << "Debug \n";
	#endif
}


void TrafficLogger::writeError(std::string error, uint currentTimestep)
{
	instance.file << "Error: " << error << "\t timestep: " << currentTimestep << "\n";
}


void TrafficLogger::writeSimulationEnd(uint numRoadCells, uint numTimesteps, double duration)
{
	uint hours = static_cast<uint>(std::floor(duration / 3600));
	uint minutes = static_cast<uint>(std::floor((duration - 3600 * hours) / 60));
	uint seconds = static_cast<uint>(duration - 3600 * hours - 60 * minutes);

	std::string durationString = std::to_string(hours) + " h \t" + std::to_string(minutes) + " m \t" + std::to_string(seconds) + " s";

	instance.file << "Simulation finished: Number of roadcells : \t" << numRoadCells << "\t total number of timesteps: " << numTimesteps;
	instance.file << "\t duration: " << duration << " s" << "\t\t => " << durationString << "\n\n";
	instance.file.close();
}
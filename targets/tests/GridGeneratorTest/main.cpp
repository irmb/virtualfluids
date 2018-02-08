#include "gmock/gmock.h"
#include "mpi.h"
#include <fstream>
#include "VirtualFluidsBasics/utilities/logger/Logger.h"

int main(int argc, char **argv)
{
	MPI_Init(&argc, &argv);

    std::ofstream file("log.txt", std::ofstream::out);
    logging::Logger::addStream(&file);

	::testing::InitGoogleTest(&argc, argv);
	int i = RUN_ALL_TESTS();

    file.close();
	MPI_Finalize();

    return i;
}

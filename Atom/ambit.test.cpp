#ifdef AMBIT_USE_MPI
#include <mpi.h>
#endif
#include "Include.h"
#include "OutStreams.h"
#include "gitInfo.h"
#include "Universal/Enums.h"
#include "Atom/GetPot"
#include <boost/math/special_functions/bessel.hpp>

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <boost/bind.hpp>

#include <gtest/gtest.h>

#ifdef AMBIT_USE_MPI
    #ifdef AMBIT_USE_SCALAPACK
    #if !(_FUS)
        #define blacs_exit_ blacs_exit
    #endif
    extern "C"{
    void blacs_exit_(const int*);
    }
    #endif
#endif

// MPI details (if not used, we can have NumProcessors == 1)
int NumProcessors;
int ProcessorRank;

// The debug options for the whole program.
Debug DebugOptions;

void PrintHelp(const std::string& ApplicationName);

int main(int argc, char* argv[])
{
    #ifdef AMBIT_USE_MPI
        MPI_Init(&argc, &argv);
        MPI_Comm_size(MPI_COMM_WORLD, &NumProcessors);
        MPI_Comm_rank(MPI_COMM_WORLD, &ProcessorRank);
    #else
        NumProcessors = 1;
        ProcessorRank = 0;
    #endif
    
    OutStreams::InitialiseStreams();
    outstream = logstream;

    testing::InitGoogleTest(&argc, argv);
    int ret = RUN_ALL_TESTS();

    OutStreams::FinaliseStreams();

    #ifdef AMBIT_USE_MPI
        MPI_Finalize();
    #endif

    return ret;
}

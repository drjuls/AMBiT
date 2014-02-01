#ifdef _MPI
#include <mpi.h>
#endif
#include "Include.h"
#include "OutStreams.h"
#include "gitInfo.h"
//#include "Atom.h"
//#include "RateCalculator.h"
#include "Universal/Enums.h"
//#include "Universal/Function.h"
//#include "Universal/Integrator.h"
#include "Atom/GetPot"
#include <boost/math/special_functions/bessel.hpp>

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <boost/bind.hpp>

#include <gtest/gtest.h>

#ifdef _MPI
    #ifdef _SCALAPACK
    #if !(_FUS)
        #define blacs_exit_ blacs_exit
    #endif
    extern "C"{
    void blacs_exit_(const int*);
    }
    #endif
#endif

// MPI details (if not used, we can have NumProcessors == 1)
unsigned int NumProcessors;
unsigned int ProcessorRank;

// The debug options for the whole program.
Debug DebugOptions;

void PrintHelp(const std::string& ApplicationName);

int main(int argc, char* argv[])
{
    #ifdef _MPI
        MPI::Init(argc, argv);
        MPI::Intracomm& comm_world = MPI::COMM_WORLD; // Alias
        NumProcessors = comm_world.Get_size();
        ProcessorRank = comm_world.Get_rank();
    #else
        NumProcessors = 1;
        ProcessorRank = 0;
    #endif
    
    OutStreams::InitialiseStreams();

    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

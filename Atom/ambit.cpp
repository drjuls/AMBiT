#ifdef _MPI
#include <mpi.h>
#endif
#include "Include.h"
#include "OutStreams.h"
#include "Atom.h"

#include <Atom/GetPot>

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

    try
    {
        GetPot lineInput(argc, argv, ",");

        // Check for help message
        if(lineInput.size() == 1 || lineInput.search(2, "--help", "-h"))
            PrintHelp(lineInput[0]);

        // Find .input file
        std::string inputFileName = "";
        const std::string extension = ".input";

        if(lineInput.search("-f"))
        {   inputFileName = lineInput.next("");
            if(inputFileName == "")
                PrintHelp(lineInput[0]);
        }
        else
        {   // Else look for .input arguments
            bool foundFile = false;
            bool endOfInput = false;
            while(!foundFile && !endOfInput)
            {
                inputFileName = lineInput.next_nominus();
                if(inputFileName == "")
                    endOfInput = true;
                else if(inputFileName.find(extension) != std::string::npos)
                    foundFile = true;
            }
            
            if(endOfInput)
            {   *errstream << "Input file not found" << std::endl;
                PrintHelp(lineInput[0]);
            }
        }

        GetPot fileInput(inputFileName.c_str(), "//", "\n", ",");
        fileInput.absorb(lineInput);

        // Must-have arguments. This also checks on the file's existence.
        int Z, N, Charge;
        Z = fileInput("Z", 0);
        N = fileInput("N", -1);  // Number of electrons
        if(Z != 0 && N == -1)
        {   Charge = fileInput("Charge", Z);
            N = Z - Charge;
        }
        if(Z == 0 || N == -1 || N > Z)
        {   *errstream << inputFileName << " must specify Z, and either N (number of electrons) or Charge.\n"
                       << "    Furthermore Z >= N." << std::endl;
            exit(1);
        }

        // Identifier
        std::string identifier = fileInput("ID", "");
        if(identifier == "")
        {
            if(inputFileName.find(extension) != std::string::npos)
                identifier = inputFileName.substr(0, inputFileName.find(extension));
            else
            {   *errstream << "No ID could be found." << std::endl;
                exit(1);
            }
        }

        // That's all we need to start the Atom class
        Atom A(fileInput, Z, N, identifier);
        A.Run();
/*
        InputParser userinput("Sm15+.input");
        *outstream << userinput.Parse() << std::endl;
        int val[20];
        unsigned int nval = userinput.GetMultipleValues("OddSymmetryTwoJ", val, 20);
        *outstream << nval;
        for(unsigned int i=0; i<nval; i++)
            *outstream << ":" << val[i];
        *outstream << std::endl;
*/
    }
    catch(std::bad_alloc& ba)
    {   *errstream << ba.what() << std::endl;
        exit(1);
    }

    #ifdef _MPI
        comm_world.Barrier();
        #ifdef _SCALAPACK
            // Continue should be non-zero, otherwise MPI_Finalise is called.
            int cont = 1;
            blacs_exit_(&cont);
        #endif

        MPI::Finalize();
    #endif

    *outstream << "\nFinished" << std::endl;
    OutStreams::FinaliseStreams();

    PAUSE
    return 0;
}

void PrintHelp(const std::string& ApplicationName)
{
    *outstream << std::endl;
    *outstream << ApplicationName << " Usage Message" << std::endl;
    exit(0);
}

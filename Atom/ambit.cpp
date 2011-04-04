#ifdef _MPI
#include <mpi.h>
#endif
#include "Include.h"
#include "OutStreams.h"
#include "Atom.h"
#include "RateCalculator.h"
#include "Universal/Enums.h"

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
        int Z, ZIterations, N, Charge;
        Z = fileInput("Z", 0);
        N = fileInput("HF/N", -1);  // Number of electrons
        ZIterations = fileInput("ZIterations", 1);
        if(Z != 0 && N == -1)
        {   Charge = fileInput("HF/Charge", -1);
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

        if(ZIterations > 1 && fileInput("NumValenceElectrons", 1) > 1)
        {
            *errstream << "Configuration interaction for multiple atoms with differing charge is not supported yet." << std::endl;
            exit(1);
        }

        // That's all we need to start the Atom class
        Atom A(fileInput, Z, N, identifier);
        A.Run();
        if(ZIterations > 1)
        {
            int ZCounter = 1;
            while(ZCounter < ZIterations) {
                Atom aA(fileInput, Z + ZCounter, N, identifier);
                aA.Run();
                ZCounter++;
            }
        }
        
        // At this stage, we are ready to run calculations with the
        // calculated energy levels

        // Print the solutions currently available after CI calculation
        DisplayOutputType::Enum OutputType;
        std::string TypeString = fileInput("CI/Output/Type", "");
        if(TypeString == "Standard")
        {
            OutputType = DisplayOutputType::Standard;
        }
        else if(TypeString == "CommaSeparated")
        {
            OutputType = DisplayOutputType::CommaSeparated;
        }
        else if(TypeString == "SpaceSeparated")
        {
            OutputType = DisplayOutputType::SpaceSeparated;
        }
        else if(TypeString == "TabSeparated")
        {
            OutputType = DisplayOutputType::TabSeparated;
        }

        if(TypeString != "" && fileInput("NumValenceElectrons", 1) > 1)
        {
            A.GetSolutionMap()->Print(OutputType);
        }
        
        A.GetSolutionMap()->PrintID();
        
        // Interactive mode
        if(fileInput.search("--interactive"))
        {
            std::string input;
            while(input != "quit")
            {
                std::cin >> input;
                if(input != "quit")
                    A.GetSolutionMap()->PrintSolution(A.GetSolutionMap()->FindByIdentifier(input));
            }
        }
        
        // Transition calculations should be performed here!

//        RateCalculator RC(A.GetBasis());
//        RC.CalculateAllDipoleStrengths(&A, Symmetry(6, odd), 0, false, true);
//        *outstream << std::endl;
//        RC.CalculateAllDipoleStrengths(&A, Symmetry(6, odd), 1, false, true);
//        *outstream << std::endl;
//        RC.CalculateAllDipoleStrengths(&A, Symmetry(8, odd), 0, false, true);
//        *outstream << std::endl;
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

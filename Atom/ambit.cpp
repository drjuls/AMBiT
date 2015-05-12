#ifdef _MPI
#include <mpi.h>
#endif
#include "Include.h"
#include "OutStreams.h"
#include "gitInfo.h"
#include "ambit.h"
#include "Atom.h"
#include "Atom/MultirunOptions.h"
#include <boost/math/special_functions/bessel.hpp>

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <boost/bind.hpp>

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
int NumProcessors;
int ProcessorRank;

// The debug options for the whole program.
Debug DebugOptions;

int main(int argc, char* argv[])
{
    #ifdef _MPI
        MPI_Init(&argc, &argv);
        MPI_Comm_size(MPI_COMM_WORLD, &NumProcessors);
        MPI_Comm_rank(MPI_COMM_WORLD, &ProcessorRank);
    #else
        NumProcessors = 1;
        ProcessorRank = 0;
    #endif

    OutStreams::InitialiseStreams();

    try
    {   MultirunOptions lineInput(argc, argv, ",");

        // Check for help message
        if(lineInput.size() == 1 || lineInput.search(2, "--help", "-h"))
        {   Ambit::PrintHelp(lineInput[0]);
            return 0;
        }

        *outstream << "AMBiT version:    " << GIT_LAST_TAG << std::endl;
        *outstream << "      git branch: " << GIT_SOURCE_DESCRIPTION << std::endl;
        *outstream << "      compiled:   " << COMPILE_DATE << std::endl;

        if(lineInput.search("--version"))
            return 0;

        // Find .input file
        std::string inputFileName = "";
        const std::string extension = ".input";

        if(lineInput.search("-f"))
        {   inputFileName = lineInput.next("");
            if(inputFileName == "")
                Ambit::PrintHelp(lineInput[0]);
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
            
            if(endOfInput && !lineInput.search("--recursive-build"))
            {   *errstream << "Input file not found" << std::endl;
                Ambit::PrintHelp(lineInput[0]);
            }
        }

        MultirunOptions fileInput(inputFileName.c_str(), "//", "\n", ",");
        fileInput.absorb(lineInput);

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

        Ambit ambit(fileInput, identifier);

        if(fileInput.search("--recursive-build"))
        {
            ambit.Recursive();
        }
        else
        {
            ambit.EnergyCalculations();

            if(!fileInput.search("--check-sizes"))
            {
                ambit.TransitionCalculations();
                ambit.Recombination();
            }
        }
    /*
        else
        {
            if(fileInput.search("--print-wf"))
            {                
                StateIterator it = A.GetCore()->GetStateIterator();
                it.First();
                while(!it.AtEnd())
                {
                    it.GetState()->Print(identifier + "_core" + it.GetOrbitalInfo().Name() + ".out", A.GetLattice());
                    it.Next();
                }
                
                it = A.GetBasis()->GetStateIterator();
                it.First();
                while(!it.AtEnd())
                {
                    it.GetState()->Print(identifier + "_excited" + it.GetOrbitalInfo().Name() + ".out", A.GetLattice());
                    it.Next();
                }
            }
            
            // At this stage, we are ready to run calculations with the
            // calculated energy levels
            if(A.NumberRunsSelected() <= 1)
            {
                A.WriteEigenstatesToSolutionMap();
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
        
                // Interactive mode
                if(fileInput.search("--interactive"))
                {
                    A.GetSolutionMap()->PrintID();
                    std::string input;
                    while(input != "quit" && input != "q")
                    {
                        std::cin >> input;
                        if(input != "quit" && input != "q")
                        {
                            if(strcspn(input.c_str(), "eo") > 0 && strlen(input.c_str()) > strcspn(input.c_str(), "eo") + 1 && strchr(input.c_str(), '>') == NULL)
                            {
                                if(A.GetSolutionMap()->FindByIdentifier(input) != A.GetSolutionMap()->end())
                                {
                                    A.GetSolutionMap()->PrintSolution(A.GetSolutionMap()->FindByIdentifier(input));
                                }
                                else
                                {
                                    *outstream << "Level " << input << " not found." << std::endl;
                                }
                            }
                            else if(strchr(input.c_str(), '>') != NULL)
                            {
                                char input_cstring[strlen(input.c_str()) + 1];
                                strcpy(input_cstring, input.c_str());
                                char* token1 = strtok(input_cstring, " ->");
                                char* token2 = strtok(NULL, " ->");
                                if(strchr(input.c_str(), '|') != NULL)
                                {
                                    strcpy(input_cstring, input.c_str());
                                    char* token3 = strtok(input_cstring, " |");
                                    char* token4 = strtok(NULL, " |");
                                    if(TransitionType::StringSpecifiesTransitionType(token4))
                                    {
                                        A.GetSolutionMap()->FindByIdentifier(token1)->second.GetTransitionSet()->insert(Transition(&A, TransitionType(token4), A.GetSolutionMap()->FindByIdentifier(token1)->first, A.GetSolutionMap()->FindByIdentifier(token2)->first));
                                    }
                                    else
                                    {
                                        A.GetSolutionMap()->FindByIdentifier(token1)->second.GetTransitionSet()->insert(Transition(&A, A.GetSolutionMap()->FindByIdentifier(token1)->first, A.GetSolutionMap()->FindByIdentifier(token2)->first));
                                    }
                                }
                                else
                                {
                                    A.GetSolutionMap()->FindByIdentifier(token1)->second.GetTransitionSet()->insert(Transition(&A, A.GetSolutionMap()->FindByIdentifier(token1)->first, A.GetSolutionMap()->FindByIdentifier(token2)->first));
                                }
                            }
                        }
                    }
                }
            }
            else
            {
                A.RunIndexBegin(true);
                SolutionMapMap SMapMap;
                while(!A.RunIndexAtEnd())
                {
                    if(!A.RestoreAllEigenstates())
                    {
                        *errstream << "int main: An error occurred while trying to restore all eigenstates." << std::endl;
                    }
                    
                    A.WriteEigenstatesToSolutionMap();
                    SMapMap.insert(std::pair<unsigned int, SolutionMap>(A.GetCurrentRunSelection()[A.GetCurrentRunIndex()],(*(A.GetSolutionMap()))));
                    
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
                    
                    if(TypeString != "")
                    {
                        A.GetSolutionMap()->Print(OutputType);
                    }
                    
                    A.RunIndexNext(true);
                }
                SMapMap.Print();
            }
        } // End not recursive build
    */

        // Transition calculations should be performed here!
        //A.GetSolutionMap()->FindByIdentifier("0e0")->second.GetTransitionSet()->insert(Transition(&A, TransitionType(MultipolarityType::E, 1), A.GetSolutionMap()->FindByIdentifier("0e0")->first, A.GetSolutionMap()->FindByIdentifier("2o0")->first));
        //A.GetSolutionMap()->FindByIdentifier("0e0")->second.GetTransitionSet()->insert(Transition(&A, TransitionType(MultipolarityType::E, 1), A.GetSolutionMap()->FindByIdentifier("0e0")->first, A.GetSolutionMap()->FindByIdentifier("2o1")->first));

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
        MPI_Barrier(MPI_COMM_WORLD);
        #ifdef _SCALAPACK
            // Continue should be non-zero, otherwise MPI_Finalise is called.
            int cont = 1;
            blacs_exit_(&cont);
        #endif

        MPI_Finalize();
    #endif

    *outstream << "\nFinished" << std::endl;
    OutStreams::FinaliseStreams();

    return 0;
}

void Ambit::EnergyCalculations()
{
    int Z = user_input("Z", 0);
    if(Z <= 0)
    {   *errstream << "Z not specified correctly." << std::endl;
        exit(1);
    }

    // Choose which of the multiple runs are being done in the current calculation
    unsigned int total_run_selections = user_input.vector_variable_size("-r");

    if((user_input.GetNumRuns() <= 1) || user_input.search("--check-sizes"))
    {   // Only do a single run if --check-sizes option is set (doesn't matter which one).
        run_indexes.push_back(0);
    }
    else if(total_run_selections)
    {
        for(int i = 0; i < total_run_selections; i++)
        {
            int selection = user_input("-r", 0, i);

            if(selection == 0)
            {   // "-r=0" is shorthand for "just do the zero variation variety"
                if(total_run_selections > 1)
                {   *outstream << "USAGE: \"-r=0\" can only be used by itself." << std::endl;
                    exit(1);
                }

                // Try to find a good candidate for -r=0
                auto pair = user_input.FindZeroParameterRun();
                if(pair.second < 0)
                {   *errstream << "ERROR: Multiple run option couldn't find zero-variation option." << std::endl;
                    exit(1);
                }
                else
                {   *outstream << "Running case " << pair.second+1 << ": " << pair.first << " = 0.\n" << std::endl;
                    run_indexes.push_back(pair.second);
                }
            }
            else if((selection < 0) || (selection > user_input.GetNumRuns()))
            {   *outstream << "USAGE: Option \"-r\" included index outside of multiple run range." << std::endl;
                exit(1);
            }
            else
            {   // Subtract one from the user's selection to make it start at zero
                run_indexes.push_back(selection-1);
            }
        }
    }
    else
    {   // If "-r" is not found, do all of the multiple run options
        for(int i = 0; i < user_input.GetNumRuns(); i++)
            run_indexes.push_back(i);
    }

    // That's all we need to start the Atom class
    for(int run: run_indexes)
    {
        user_input.SetRun(run);
        std::string id = identifier + "_" + itoa(run);
        atoms.emplace_back(user_input, Z, id);      // This copies the user_input, so each atom has its own.
    }

    // Start with -r=0 if possible, otherwise just first one
    int zero_index = user_input.FindZeroParameterRun().second;
    first_run_index = run_indexes.size();
    if(zero_index > 0)
    {   // Check whether it is in current runs
        for(first_run_index = 0; first_run_index < run_indexes.size(); first_run_index++)
            if(run_indexes[first_run_index] == zero_index)
                break;
    }
    if(first_run_index >= run_indexes.size())   // -r=0 not found
        first_run_index = 0;

    // Log first build
    DebugOptions.LogFirstBuild(true);
    DebugOptions.LogHFIterations(true);
    DebugOptions.OutputHFExcited(true);
    DebugOptions.HartreeEnergyUnits(true);

    user_input.SetRun(run_indexes[first_run_index]);        // Just for the printing (each atom has its own user_input)
    user_input.PrintCurrentRunCondition(*outstream, "\n");
    pCoreConst hf_open_core = atoms[first_run_index].MakeBasis();
    *outstream << std::endl;

    // Don't output the others
    DebugOptions.LogFirstBuild(false);
    DebugOptions.LogHFIterations(false);
    DebugOptions.OutputHFExcited(false);

    for(unsigned int index = 0; index < run_indexes.size(); index++)
        if(index != first_run_index)
            atoms[index].MakeBasis(hf_open_core);

    // Print basis only option
//    if(user_input.search(2, "--print-basis", "-p"))
//    {   // Check follower for option
//        std::string print_option = user_input.next("");
//        if(print_option == "Cowan")
//            atom.GenerateCowanInputFile();
//        else
//            atom.WriteGraspMCDF();
//    }

    // Generate Brueckner orbitals and MBPT integrals
    for(auto& atom: atoms)
        atom.MakeMBPTIntegrals();

    // CI
    if(user_input.search("--check-sizes"))
        atoms[0].CheckMatrixSizes();
    else
    {   pLevelMap levels = atoms[0].ChooseHamiltoniansAndRead();

        for(int i = 1; i < run_indexes.size(); i++)
        {   user_input.SetRun(run_indexes[i]);
            atoms[i].ChooseHamiltoniansAndRead();
        }

        for(auto& pair: *levels)
        {
            pHamiltonianID sym = pair.first;

            for(int i = 0; i < run_indexes.size(); i++)
            {
                // This if statement is just to switch off printing the run condition for only one run
                if(user_input.GetNumRuns() > 1)
                {   user_input.SetRun(run_indexes[i]);
                    user_input.PrintCurrentRunCondition(*outstream, "\n");
                }
                atoms[i].CalculateEnergies(sym);
            }
        }
    }
}

void Ambit::TransitionCalculations()
{
    TransitionMap transitions(atoms[first_run_index]);

    // Calculate default types
    std::vector<std::pair<std::string, TransitionType>> default_types =
           {{"E1", TransitionType(MultipolarityType::E, 1)},
            {"M1", TransitionType(MultipolarityType::M, 1)},
            {"E2", TransitionType(MultipolarityType::E, 2)},
            {"M2", TransitionType(MultipolarityType::M, 2)},
            {"E3", TransitionType(MultipolarityType::E, 3)},
           };

    for(const auto& type_pair: default_types)
    {
        std::string user_string = "Transitions/" + type_pair.first;
        int num_transitions = user_input.vector_variable_size(user_string.c_str());
        for(int i = 0; i < num_transitions; i++)
        {
            std::string transition = user_input(user_string.c_str(), "", i);
            int pos = transition.find("->");
            if(pos == std::string::npos)
                *errstream << "TransitionCalculations: " << user_string << " = " << transition << " not properly formed." << std::endl;
            else
            {   LevelID left(make_LevelID(transition.substr(0, pos)));
                LevelID right(make_LevelID(transition.substr(pos+2)));

                transitions.CalculateTransition(left, right, type_pair.second);
            }
        }
    }

    // Calculate non-identified transitions
    int num_transitions = user_input.vector_variable_size("Transitions/Transitions");
    for(int i = 0; i < num_transitions; i++)
    {
        std::string transition = user_input("Transitions/Transitions", "", i);
        int pos = transition.find("->");
        if(pos == std::string::npos)
            *errstream << "TransitionCalculations: " << transition << " not properly formed." << std::endl;
        else
        {   LevelID left(make_LevelID(transition.substr(0, pos)));
            LevelID right(make_LevelID(transition.substr(pos+2)));

            transitions.CalculateTransition(left, right);
        }
    }

    // Calculate all transitions of a certain type below a given energy
    pLevelMap levels = atoms[first_run_index].GetLevels();
    for(const auto& type_pair : default_types)
    {
        std::string user_string = "Transitions/All" + type_pair.first + "Below";
        double max_energy = user_input(user_string.c_str(), 0.0);
        if(max_energy)
        {
            auto left_it = levels->begin();
            while(left_it != levels->end())
            {
                auto right_it = left_it;
                right_it++;
                while(right_it != levels->end())
                {
                    if(transitions.TransitionExists(left_it->first->GetSymmetry(), right_it->first->GetSymmetry(), type_pair.second))
                    {
                        for(int i = 0; i < left_it->second.size(); i++)
                        {
                            if(left_it->second[i]->GetEnergy() > max_energy)
                                break;

                            for(int j = 0; j < right_it->second.size(); i++)
                            {
                                if(right_it->second[i]->GetEnergy() > max_energy)
                                    break;

                                transitions.CalculateTransition(std::make_pair(left_it->first, i),
                                                                std::make_pair(right_it->first, j), type_pair.second);

                            }
                        }
                    }
                    right_it++;
                }
                left_it++;
            }
        }
    }

    if(transitions.size())
        transitions.Print();
}

void Ambit::Recombination()
{
    std::string target_file = user_input("DR/Target/Filename", "");
    if(target_file.length() == 0)
        return;

    int target_two_j = user_input("DR/Target/TwoJ", -1);
    if(target_two_j < 0)
    {   *errstream << "DR/Target/TwoJ must be set in input file." << std::endl;
        exit(1);
    }

    std::string parity_string = user_input("DR/Target/Parity", "");
    Parity target_parity;
    if(parity_string == "even")
        target_parity = Parity::even;
    else if (parity_string == "odd")
        target_parity = Parity::odd;
    else
    {   *errstream << "DR/Target/Parity must be set in input file to either \"even\" or \"odd\"." << std::endl;
        exit(1);
    }

    // Get target atom
    MultirunOptions target_input(target_file.c_str(), "//", "\n", ",");
    std::string target_id = target_input("ID", "");
    if(target_id == "")
    {   if(target_file.find(".input") != std::string::npos)
           target_id = target_file.substr(0, target_file.find(".input"));
       else
       {   *errstream << "No ID could be found." << std::endl;
           exit(1);
       }
    }

    Ambit target_calculator(target_input, target_id);
    target_calculator.EnergyCalculations();

    Atom& target_atom = target_calculator.atoms[target_calculator.first_run_index];

    // Get target level
    pHamiltonianID target_key = std::make_shared<HamiltonianID>(target_two_j, target_parity);
    const auto& target_it = target_atom.GetLevels()->find(target_key);
    if(target_it == target_atom.GetLevels()->end() || target_it->second.size() == 0)
    {   *errstream << "Recombination: Cannot find target level." << std::endl;
        exit(1);
    }
    pLevelConst target_level = target_it->second[0];

    // Get widths of current levels
    atoms[first_run_index].Autoionization(target_level);
}

void Ambit::PrintHelp(const std::string& ApplicationName)
{
    *outstream << ApplicationName << std::endl;
    *outstream
        << "usage: ambit [--version] [-h|--help] [-f] <file.input>\n"
        << "             <options> <atomic data and commands>\n"
        << "\n"
        << "Input file must be specified either with \"-f anyfile\" or \"file.input\".\n"
        << "Options may be specified in the command line or input file.\n"
        << "Most commands, data, and options are structured.\n\n"
        << "Important data:\n"
        << "   ID =       identifier for all save files\n"
        << "   Z  =       nuclear charge\n"
        << "   NumValenceElectrons =\n"
        << "              determines the type of calculation to run (CI or MBPT)\n"
        << "   HF/N =     number of electrons included in HF potential\n"
        << "\n"
        << "Common data structures in the input file include:\n"
        << "   [Lattice]  numerical grid\n"
        << "   [HF]       (mandatory) details for how to perform Hartree-Fock\n"
        << "   [Basis]    basis set (e.g. B-splines)\n"
        << "   [CI]       configuration interaction\n"
        << "   [MBPT]     many-body perturbation theory\n"
        << "\n"
        << "Options include:\n"
        << "   -d|--dont-save   don't save any files\n"
        << "   -s12             include one-body and two-body MBPT diagrams\n"
        << "   --check-sizes    calculate CI matrix size and/or\n"
        << "                        number of MBPT integrals\n"
        << "   --generate-integrals-mbpt\n"
        << "                    calculate specified MBPT integrals\n\n"
        << "Sample input files are included with the documentation." << std::endl;
}

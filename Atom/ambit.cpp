#ifdef AMBIT_USE_MPI
#include <mpi.h>
#endif
#include "Include.h"
#include "OutStreams.h"
#include "gitInfo.h"
#include "ambit.h"
#include "Atom.h"
#include "ExternalField/EJOperator.h"
#include "ExternalField/Hyperfine.h"

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

    try
    {   MultirunOptions lineInput(argc, argv, ",");

        // Check for help message
        if(lineInput.size() == 1 || lineInput.search(2, "--help", "-h"))
        {   Ambit::PrintHelp(lineInput[0]);
            return 0;
        }

        if(ProcessorRank == 0)
        {
            *outstream << "AMBiT version:    " << GIT_LAST_TAG << std::endl;
            *outstream << "      git branch: " << GIT_SOURCE_DESCRIPTION << std::endl;
            *outstream << "      compiled:   " << COMPILE_DATE << std::endl;

            // Print git diff output
            if(strlen(GIT_DIFF_OUTPUT))
            {   *outstream << "\n----------------------------------------------------------"
                << "\ngit diff:"
                << "\n----------------------------------------------------------\n";
                *outstream << GIT_DIFF_OUTPUT;
            }
        }

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
    }
    catch(std::bad_alloc& ba)
    {   *errstream << ba.what() << std::endl;
        exit(1);
    }

    #ifdef AMBIT_USE_MPI
        MPI_Barrier(MPI_COMM_WORLD);
        #ifdef AMBIT_USE_SCALAPACK
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

Ambit::Ambit(MultirunOptions& user_input, const std::string& identifier):
    user_input(user_input), identifier(identifier)
{
    if(ProcessorRank == 0)
    {
        *outstream << "----------------------------------------------------------"
                   << "\ninput:"
                   << "\n----------------------------------------------------------\n";
        user_input.print();
        *outstream << "----------------------------------------------------------" << std::endl;
    }
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
    *outstream << "Core orbitals: \n";
    hf_open_core->Print();

    *outstream << "Excited orbitals: \n";
    atoms[first_run_index].GetBasis()->excited->Print();
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
    {   pLevelStore levels = atoms[0].ChooseHamiltoniansAndRead();
        pAngularDataLibrary angular_data_lib = atoms[0].GetAngularDataLibrary();

        // Don't calculate extra levels
        if(user_input.search(2, "--ci-complete", "--CI-complete"))
            return;

        for(int i = 1; i < run_indexes.size(); i++)
        {   user_input.SetRun(run_indexes[i]);
            atoms[i].ChooseHamiltoniansAndRead(angular_data_lib);
        }

        for(auto& key: *levels)
        {
            for(int i = 0; i < run_indexes.size(); i++)
            {
                // This if statement is just to switch off printing the run condition for only one run
                if(user_input.GetNumRuns() > 1)
                {   user_input.SetRun(run_indexes[i]);
                    user_input.PrintCurrentRunCondition(*outstream, "\n");
                }
                atoms[i].CalculateEnergies(key);
            }
        }
    }
}

// Used below to run all transition types. Can only be used if they have constructor
//      TransitionCalculatorType(user_input, atom)
#define RUN_AND_STORE_TRANSITION(id, TRANSITION_CALCULATOR_TYPE) \
do { \
  std::string prefix = std::string("Transitions/") + (#id); \
  if(user_input.SectionExists(prefix)) \
  { \
    user_input.set_prefix(prefix); \
    std::unique_ptr<TRANSITION_CALCULATOR_TYPE> calculon(new TRANSITION_CALCULATOR_TYPE(user_input, atom)); \
    calculon->CalculateAndPrint(); \
    calculators.push_back(std::move(calculon)); \
  } \
} while(0)

void Ambit::TransitionCalculations()
{
    Atom& atom = atoms[first_run_index];

    // EM types
    std::array<std::tuple<std::string, MultipolarityType, int>, 5> EM_transition_types =
           {std::make_tuple("E1", MultipolarityType::E, 1),
            std::make_tuple("M1", MultipolarityType::M, 1),
            std::make_tuple("E2", MultipolarityType::E, 2),
            std::make_tuple("M2", MultipolarityType::M, 2),
            std::make_tuple("E3", MultipolarityType::E, 3),
           };

    std::vector<std::unique_ptr<TransitionCalculator>> calculators;

    for(const auto& types: EM_transition_types)
    {
        std::string user_string = "Transitions/" + std::get<0>(types);

        if(user_input.SectionExists(user_string))
        {
            user_input.set_prefix(user_string);
            std::unique_ptr<EMCalculator> calculon(new EMCalculator(std::get<1>(types), std::get<2>(types), user_input, atom.GetBasis(), atom.GetLevels(), atom.GetHFOperator()->GetIntegrator()));

            calculon->CalculateAndPrint();
            calculators.push_back(std::move(calculon));
        }
    }

    // Other operators
    RUN_AND_STORE_TRANSITION(HFS1, HyperfineDipoleCalculator);

    user_input.set_prefix("");

    for(auto& calc: calculators)
    {
        calc->PrintAll();
    }
}

#undef RUN_AND_STORE_TRANSITION

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

    *outstream << "\nTarget:" << std::endl;
    Ambit target_calculator(target_input, target_id);
    target_calculator.EnergyCalculations();
    *outstream << "----------------------------------------------------------" << std::endl;

    Atom& target_atom = target_calculator.atoms[target_calculator.first_run_index];

    // Get target level
    pHamiltonianID target_key = std::make_shared<HamiltonianID>(target_two_j, target_parity);
    LevelVector target_level_vec = target_atom.GetLevels()->GetLevels(target_key);
    int target_level_index = user_input("DR/Target/Index", 0);
    if(target_level_vec.size() <= target_level_index)
    {   *errstream << "Recombination: Cannot find target level." << std::endl;
        exit(1);
    }
    pLevelConst target_level = target_level_vec[target_level_index];

    // Get widths of current levels
    if(user_input.search("--configuration-average"))
        atoms[first_run_index].AutoionizationConfigurationAveraged(target_level);
    else if(user_input.search("DR/--energy-grid"))
        atoms[first_run_index].AutoionizationEnergyGrid(target_level);
    else
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
        << "   -s12             include one-body and two-body MBPT diagrams\n"
        << "   --check-sizes    calculate CI matrix size and/or\n"
        << "                        number of MBPT integrals\n"
        << "   --generate-integrals-mbpt\n"
        << "                    calculate specified MBPT integrals\n\n"
        << "Sample input files are included with the documentation." << std::endl;
}

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
#include "ExternalField/FieldShift.h"
#include "ExternalField/RadiativePotential.h"
#include "ExternalField/YukawaPotential.h"
#include "ExternalField/KineticEnergy.h"
#include "ExternalField/LorentzInvarianceT2.h"
#include "ExternalField/NormalMassShiftDecorator.h"

// Headers to get stack-traces
#ifdef UNIX
#include <unistd.h>
#include <execinfo.h>
#include <signal.h>
#endif

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
int Ambit::NumProcessors;
int Ambit::ProcessorRank;

// The debug options for the whole program.
Ambit::Debug Ambit::DebugOptions;

int main(int argc, char* argv[])
{
    using namespace Ambit;

    #ifdef AMBIT_USE_MPI
        #ifdef AMBIT_USE_OPENMP
            int MPI_thread_safety;
            // Tell MPI to use thread-safe functions
            MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &MPI_thread_safety);
            MPI_Comm_size(MPI_COMM_WORLD, &NumProcessors);
            MPI_Comm_rank(MPI_COMM_WORLD, &ProcessorRank);

            /* Check if the MPI implementation provides the required level of thread safety (only the master thread needs to make calls to MPI subroutines so this should almost always be supported)
            */
            if(MPI_thread_safety != MPI_THREAD_FUNNELED){
                std::cerr << "Warning: this MPI implementation is incompatible OpenMP." << std::endl;
                exit(-1);
            }
        #else
            MPI_Init(&argc, &argv);
            MPI_Comm_size(MPI_COMM_WORLD, &NumProcessors);
            MPI_Comm_rank(MPI_COMM_WORLD, &ProcessorRank);
        #endif
    #else
        NumProcessors = 1;
        ProcessorRank = 0;
    #endif

    // Set the file permissions so generated files can be read and written to by group, as well as user
    // (so AngularData files can be reused by different users) and register our custom signal handler
    // for SIGTERM and SIGSEGV.
    // N.B. these currently only work on posix systems.
    #ifdef UNIX
        umask(0003);

        // NOTE: these print to stderr rather than logstream, since that is the place stack-traces would
        // be written to
        if(signal(SIGTERM, ::signal_handler) == SIG_ERR){
            fputs("Note: This system does not allow AMBiT to produce stack-traces on termination", stderr);
        }
        if(signal(SIGSEGV, ::signal_handler) == SIG_ERR){
            fputs("Note: This system does not allow AMBiT to produce stack-traces on a segmentation fault", stderr);
        }
    #endif

    // Check whether the environment variable OMP_NUM_THREADS has been explicitly set - OpenMP threading
    // should only be used if explicitly requested by the user
    #ifdef AMBIT_USE_OPENMP
        if(!std::getenv("OMP_NUM_THREADS"))
            omp_set_num_threads(1);
    #endif

    OutStreams::InitialiseStreams();

    // Write the number of threads to the log file. This is useful for debugging.
    #ifdef AMBIT_USE_OPENMP
    *logstream << "AMBiT running with " << omp_get_max_threads() << " OpenMP threads." << std::endl;
    #endif

    try
    {   MultirunOptions lineInput(argc, argv, ",");

        // Check for help message
        if(lineInput.size() == 1 || lineInput.search(2, "--help", "-h"))
        {   Ambit::AmbitInterface::PrintHelp(lineInput[0]);
            return 0;
        }

        if(ProcessorRank == 0)
        {
            *outstream << "AMBiT  Copyright (C) 2018  Julian Berengut\n"
                       << "This program comes with ABSOLUTELY NO WARRANTY; for details see LICENSE.\n\n";
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
                Ambit::AmbitInterface::PrintHelp(lineInput[0]);
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
                Ambit::AmbitInterface::PrintHelp(lineInput[0]);
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

        Ambit::AmbitInterface ambit(fileInput, identifier);
        ambit.EnergyCalculations();

        if(!fileInput.search("--check-sizes"))
        {
            ambit.TransitionCalculations();
            ambit.Recombination();
            ambit.InternalConversion();
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

namespace Ambit
{
AmbitInterface::AmbitInterface(MultirunOptions& user_input, const std::string& identifier):
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

void AmbitInterface::EnergyCalculations()
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

    // Print basis option
    if(user_input.search(2, "--print-basis", "-p"))
    {
        auto fp = file_err_handler->fopen("Orbitals.txt", "wt");

        auto lattice = atoms[first_run_index].GetLattice();
        for(auto& var: *atoms[first_run_index].GetBasis()->all)
        {
            // Print pqn, kappa, Energy, number of points
            fprintf(fp, "  %i %i %12.6f %i\n", var.first.PQN(), var.first.Kappa(), var.second->Energy(), var.second->size());
            var.second->Print(fp, lattice, true);
        }
        file_err_handler->fclose(fp);
    }
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

        for(auto& key: levels->keys)
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
    TRANSITION_CALCULATOR_TYPE calculon(user_input, atom); \
    calculon.CalculateAndPrint(); \
    calculon.PrintAll(); \
  } \
} while(0)

void AmbitInterface::TransitionCalculations()
{
    Atom& atom = atoms[first_run_index];

    // EM types
    std::array<std::tuple<std::string, MultipolarityType, int>, 8> EM_transition_types =
           {std::make_tuple("E1", MultipolarityType::E, 1),
            std::make_tuple("M1", MultipolarityType::M, 1),
            std::make_tuple("E2", MultipolarityType::E, 2),
            std::make_tuple("M2", MultipolarityType::M, 2),
            std::make_tuple("E3", MultipolarityType::E, 3),
            std::make_tuple("M3", MultipolarityType::M, 3),
            std::make_tuple("E4", MultipolarityType::E, 4),
            std::make_tuple("M4", MultipolarityType::M, 4),
           };

    std::vector<std::unique_ptr<EMCalculator>> calculators;

    for(const auto& types: EM_transition_types)
    {
        std::string user_string = "Transitions/" + std::get<0>(types);

        if(user_input.SectionExists(user_string))
        {
            user_input.set_prefix(user_string);
            std::unique_ptr<EMCalculator> calculon(new EMCalculator(std::get<1>(types), std::get<2>(types), user_input, atom));

            calculon->CalculateAndPrint();
            calculators.push_back(std::move(calculon));
        }
    }

    // Print all EM transitions
    for(auto& calc: calculators)
    {
        user_input.set_prefix(std::string("Transitions/") + calc->Name());
        calc->PrintAll();
    }

    // Other operators
    RUN_AND_STORE_TRANSITION(HFS1, HyperfineDipoleCalculator);
    RUN_AND_STORE_TRANSITION(HFS2, HyperfineQuadrupoleCalculator);
    RUN_AND_STORE_TRANSITION(HFI, GeneralisedHyperfineCalculator);
    RUN_AND_STORE_TRANSITION(FS, FieldShiftCalculator);
    RUN_AND_STORE_TRANSITION(QED, QEDCalculator);
    RUN_AND_STORE_TRANSITION(Yukawa, YukawaCalculator);
    RUN_AND_STORE_TRANSITION(KE, KineticEnergyCalculator);
    RUN_AND_STORE_TRANSITION(LLIT2, LorentzInvarianceT2Calculator);
    RUN_AND_STORE_TRANSITION(NMS, NormalMassShiftCalculator);

    user_input.set_prefix("");
}

#undef RUN_AND_STORE_TRANSITION

void AmbitInterface::Recombination()
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
    AmbitInterface target_calculator(target_input, target_id);
    target_calculator.EnergyCalculations();
    *outstream << "----------------------------------------------------------" << std::endl;

    Atom& target_atom = target_calculator.atoms[target_calculator.first_run_index];

    // Get target level
    pHamiltonianID target_key = std::make_shared<HamiltonianID>(target_two_j, target_parity);
    LevelVector target_level_vec = target_atom.GetLevels()->GetLevels(target_key);

    int target_level_index = user_input("DR/Target/Index", 0);
    if(target_level_vec.levels.size() <= target_level_index)
    {   *errstream << "Recombination: Cannot find target level." << std::endl;
        exit(1);
    }

    auto keep = std::next(target_level_vec.levels.begin(), target_level_index);
    target_level_vec.levels.erase(target_level_vec.levels.begin(), keep);
    target_level_vec.levels.resize(1);

    // Get widths of current levels
    if(user_input.search("--configuration-average"))
        atoms[first_run_index].AutoionizationConfigurationAveraged(target_level_vec);
    else if(user_input.search("DR/--energy-grid"))
        atoms[first_run_index].AutoionizationEnergyGrid(target_level_vec);
    else
        atoms[first_run_index].Autoionization(target_level_vec);
}

void AmbitInterface::InternalConversion()
{
    std::string source_file = user_input("IC/Source/Filename", "");
    if(source_file.length() == 0)
        return;

    int source_two_j = user_input("IC/Source/TwoJ", -1);
    if(source_two_j < 0)
    {   *errstream << "IC/Source/TwoJ must be set in input file." << std::endl;
        exit(1);
    }

    std::string parity_string = user_input("IC/Source/Parity", "");
    Parity source_parity;
    if(parity_string == "even")
        source_parity = Parity::even;
    else if (parity_string == "odd")
        source_parity = Parity::odd;
    else
    {   *errstream << "IC/Source/Parity must be set in input file to either \"even\" or \"odd\"." << std::endl;
        exit(1);
    }

    // Get source atom and level
    MultirunOptions source_input(source_file.c_str(), "//", "\n", ",");
    std::string source_id = source_input("ID", "");
    if(source_id == "")
    {   if(source_file.find(".input") != std::string::npos)
           source_id = source_file.substr(0, source_file.find(".input"));
       else
       {   *errstream << "No ID could be found." << std::endl;
           exit(1);
       }
    }

    *outstream << "\nSource:" << std::endl;
    AmbitInterface source_calculator(source_input, source_id);
    source_calculator.EnergyCalculations();
    *outstream << "----------------------------------------------------------" << std::endl;

    Atom& source_atom = source_calculator.atoms[source_calculator.first_run_index];

    // Get source level
    pHamiltonianID source_key = std::make_shared<HamiltonianID>(source_two_j, source_parity);
    LevelVector source_level_vec = source_atom.GetLevels()->GetLevels(source_key);

    int source_level_index = user_input("IC/Source/Index", 0);
    if(source_level_vec.levels.size() <= source_level_index)
    {   *errstream << "IC: Cannot find source level." << std::endl;
        exit(1);
    }

    auto keep = std::next(source_level_vec.levels.begin(), source_level_index);
    source_level_vec.levels.erase(source_level_vec.levels.begin(), keep);
    source_level_vec.levels.resize(1);

    if(user_input.search("--configuration-average"))
        atoms[first_run_index].InternalConversionConfigurationAveraged(source_level_vec);
    else
        atoms[first_run_index].InternalConversion(source_level_vec);
}

void AmbitInterface::PrintHelp(const std::string& ApplicationName)
{
    *outstream << ApplicationName << std::endl;
    *outstream
        << "AMBiT  Copyright (C) 2018  Julian Berengut\n"
        << "This program comes with ABSOLUTELY NO WARRANTY; for details see LICENSE.\n\n"
        << "usage: ambit [--version] [-f] <file.input>\n"
        << "             <options> <atomic data and commands>\n"
        << "\n"
        << "Input file must be specified either with \"-f anyfile\" or \"file.input\".\n"
        << "Options may be specified in the command line or input file.\n"
        << "Consult the user guide in ./Documentation for usage information.\n"
        << "Sample input files are included with the documentation." << std::endl;
}
}

#ifdef UNIX
static void signal_handler(int signo){
    /* Function to catch SIGTERM and print a stacktrace. This doesn't neatly handle mutlithreading, but
     * there isn't really any neat way to do this anyway. The stack-traces produced by this function are
     * very rudimentary and need to be passed through addr2line or a similar utility to obtain
     * line-numbers (although it will include function names in the stack-trace). Finally, this is
     * basically useless if AMBiT was not compiled with debugging symbols, as the stack-trace will just
     * be a bunch of unadorned hex addresses.
     */
    if(signo == SIGTERM){
        void *buff[256];
        size_t numFrames;

        const char* errmsg = "AMBiT terminated prematurely. Stacktrace:\n";
        // C/C++ stdio functions are not safe to call in a signal handler, so use a bare write() instead
        write(STDERR_FILENO, errmsg, strlen(errmsg));

        // Write a stack-trace. This is very rudimentary and needs to be passed through addr2line to be
        // properly readable
        numFrames = backtrace(buff, 256);
        backtrace_symbols_fd(buff, numFrames, STDERR_FILENO);

        _exit(EXIT_FAILURE);
    } else if(signo == SIGSEGV){
        void *buff[256];
        size_t numFrames;

        const char* errmsg = "AMBiT terminated on a segmentation fault. Stacktrace:\n";
        // C/C++ stdio functions are not safe to call in a signal handler, so use a bare write() instead
        write(STDERR_FILENO, errmsg, strlen(errmsg));

        // Write a stack-trace. This is very rudimentary and needs to be passed through addr2line to be
        // properly readable
        numFrames = backtrace(buff, 256);
        backtrace_symbols_fd(buff, numFrames, STDERR_FILENO);
        _exit(EXIT_FAILURE);
    }
}
#endif

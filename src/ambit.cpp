#ifdef AMBIT_USE_MPI
#include <mpi.h>
#endif
#include "Include.h"
#include "Atom/OutStreams.h"
#include "gitInfo.h"
#include "ambit.h"
#include "Atom/Atom.h"
#include "ExternalField/EJOperator.h"
#include "ExternalField/Hyperfine.h"
#include "ExternalField/FieldShift.h"
#include "ExternalField/RadiativePotential.h"
#include "ExternalField/YukawaPotential.h"
#include "ExternalField/KineticEnergy.h"
#include "ExternalField/LorentzInvarianceT2.h"
#include "ExternalField/NormalMassShiftDecorator.h"
#include "Specification/Specification.h"
#include "Specification/TransitionSpecifications.h"
#include "Specification/SpecificationMap.h"
#include "Universal/Profiler.h"

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

    // Start the global timer
    auto profiler = Profiler::Instance();
    profiler->total.start();

    // Write the number of threads to the log file. This is useful for debugging.
    #ifdef AMBIT_USE_OPENMP
    *logstream << "AMBiT running with " << omp_get_max_threads() << " OpenMP threads." << std::endl;
    #endif

    try
    {   if(ProcessorRank == 0)
        {
            *outstream << "AMBiT  Copyright (C) 2018-2026  Julian Berengut\n"
                       << "This program comes with ABSOLUTELY NO WARRANTY; for details see LICENSE.\n";
            *outstream << "\nAMBiT version:    " << GIT_LAST_TAG;
            *outstream << "\n      git branch: " << GIT_SOURCE_DESCRIPTION;
            *outstream << "\n      compiled:   " << COMPILE_DATE << std::endl;

            // Print git diff output
            if(strlen(GIT_DIFF_OUTPUT))
            {   *outstream << "\n----------------------------------------------------------"
                           << "\ngit diff:"
                           << "\n----------------------------------------------------------\n";
                *outstream << GIT_DIFF_OUTPUT << std::flush;
            }
        }

        // Find .input file
        std::string input_file_name;
        const std::string extension = ".input";
        std::vector<std::string> line_args;

        // Grab super-urgent argv options, wait for the rest
        for(char **argp = &argv[1]; *argp; ++argp)
        {
            std::string_view arg(*argp);

            if(arg == "-h" || arg == "--help")
            {   Ambit::AmbitInterface::PrintHelp(argv[0]);
                return 0;
            }
            else if(arg == "--version")
                return 0;
            else if(arg == "-f")
            {   ++argp;
                if (!*argp)
                    break;
                else
                    input_file_name = arg;
            }
            else if(arg.length() >= extension.length() && arg.substr(arg.length()-extension.length()) == extension)
                input_file_name = arg;
            else
                line_args.emplace_back(arg);
        }

        if(input_file_name.empty())
        {   Ambit::AmbitInterface::PrintHelp(argv[0]);
            return 0;
        }

        // We print the user specification as we read it into GlobalSpecification.
        *outstream << "\n----------------------------------------------------------"
                   << "\ninput (" << input_file_name << " and command line):"
                   << "\n----------------------------------------------------------" << std::endl;

        // Read input into parapara-based parameter handling.
        auto specification = GlobalSpecification::Instance();

        std::string perr = importSpecificationFile(specification, input_file_name);
        if (!perr.empty()) {
            *errstream << "importSpecificationFile:\n" << perr << std::endl;
            exit(1);
        }
        for(auto arg: line_args)
        {
            perr = importSpecificationKV(specification, arg);
            if(!perr.empty())
            {
                *errstream << "importSpecificationKV:\n" << perr << std::endl;
                exit(1);
            }
        }
        *outstream << "----------------------------------------------------------" << std::endl;

        // Validate and normalise the spec
        perr = specification->validateAndNormaliseSpecification();
        if (!perr.empty()) {
            *errstream << "validateAndNormaliseSpecification:\n" << perr << std::endl;
            exit(1);
        }

        Ambit::AmbitInterface ambit(specification->Z, specification->ID);
        ambit.EnergyCalculations();

        if(!specification->check_sizes)
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
    #endif

    // Print timing information now that all processes have finished
    profiler->total.stop();
    *logstream << "\n|=========================================================|" << std::endl;
    *logstream << absl::StreamFormat("|%-35s %21d|", "Profiling information for rank", ProcessorRank)
               << std::endl;

    *logstream << "|=========================================================|"
               << std::endl;
    *logstream << absl::StreamFormat("|%-35s", "Hartree-Fock and basis set")
               << absl::StreamFormat("|%20gs|", absl::ToDoubleSeconds(profiler->hf.elapsed()))
               << std::endl;

    *logstream << "|---------------------------------------------------------|"
               << std::endl;
    *logstream << absl::StreamFormat("|%-35s", "Two-electron Slater integrals")
               << absl::StreamFormat("|%20gs|", absl::ToDoubleSeconds(profiler->slater.elapsed()))
               << std::endl;

    *logstream << "|---------------------------------------------------------|"
               << std::endl;
    *logstream << absl::StreamFormat("|%-35s", "One-body MBPT integrals")
               << absl::StreamFormat("|%20gs|", absl::ToDoubleSeconds(profiler->one_body_mbpt.elapsed()))
               << std::endl;

    *logstream << "|---------------------------------------------------------|"
               << std::endl;
    *logstream << absl::StreamFormat("|%-35s", "Two-body MBPT integrals")
               << absl::StreamFormat("|%20gs|", absl::ToDoubleSeconds(profiler->two_body_mbpt.elapsed()))
               << std::endl;

    *logstream << "|---------------------------------------------------------|"
               << std::endl;
    *logstream << absl::StreamFormat("|%-35s", "CI")
               << absl::StreamFormat("|%20gs|", absl::ToDoubleSeconds(profiler->ci.elapsed()))
               << std::endl;

    *logstream << "|---------------------------------------------------------|"
               << std::endl;
    *logstream << absl::StreamFormat("|%-35s", "Transition matrix elements")
               << absl::StreamFormat("|%20gs|", absl::ToDoubleSeconds(profiler->transitions.elapsed()))
               << std::endl;

    *logstream << "|=========================================================|"
               << std::endl;
    *logstream << absl::StreamFormat("|%-35s", "Total wall time")
               << absl::StreamFormat("|%20gs|", absl::ToDoubleSeconds(profiler->total.elapsed()))
               << std::endl;

    *logstream << "|=========================================================|"
               << std::endl;

    #ifdef AMBIT_USE_MPI
        MPI_Finalize();
    #endif

    *outstream << "\nFinished" << std::endl;
    OutStreams::FinaliseStreams();

    return 0;
}

namespace Ambit
{
void AmbitInterface::EnergyCalculations()
{
    auto specification= GlobalSpecification::Instance();

    // Log first build
    DebugOptions.LogFirstBuild(true);
    DebugOptions.LogHFIterations(true);
    DebugOptions.OutputHFExcited(true);
    DebugOptions.HartreeEnergyUnits(true);

    pCoreConst hf_open_core = atom.MakeBasis();
    *outstream << "Core orbitals: \n";
    hf_open_core->Print();

    *outstream << "Excited orbitals: \n";
    atom.GetBasis()->excited->Print();
    *outstream << std::endl;

    // Print basis option
    if(specification->print_basis)
    {
        auto fp = file_err_handler->fopen("Orbitals.txt", "wt");

        auto lattice = atom.GetLattice();
        for(auto& var: *atom.GetBasis()->all)
        {
            // Print pqn, kappa, Energy, number of points
            fprintf(fp, "  %i %i %12.6f %i\n", var.first.PQN(), var.first.Kappa(), var.second->Energy(), var.second->size());
            var.second->Print(fp, lattice, true);
        }
        file_err_handler->fclose(fp);
    }

    // Generate MBPT integrals and Brueckner orbitals
    bool make_new_integrals = !specification->no_new_mbpt;
    bool make_brueckner = specification->mbpt_brueckner;
    bool mbpt_use_breuckner = specification->mbpt_use_brueckner;
    bool check_sizes = specification->check_sizes;

    if((!mbpt_use_breuckner && make_new_integrals) || check_sizes)
        atom.MakeMBPTIntegrals();

    if (make_brueckner && !check_sizes)
    {   // Show details of HF iterated Brueckner orbitals
        DebugOptions.OutputHFExcited(true);
        atom.GenerateBruecknerOrbitals(make_new_integrals);
    }

    if(mbpt_use_breuckner && make_new_integrals && !check_sizes)
        atom.MakeMBPTIntegrals();

    // CI
    if(check_sizes)
        atom.CheckMatrixSizes();
    else
    {   pLevelStore levels = atom.ChooseHamiltoniansAndRead();
        pAngularDataLibrary angular_data_lib = atom.GetAngularDataLibrary();

        // Don't calculate extra levels
        if(specification->ci_complete)
            return;

        // Start the CI timer
        auto profiler = Profiler::Instance();
        profiler->ci.start();

        atom.ChooseHamiltoniansAndRead(angular_data_lib);

        for(auto& key: levels->keys)
        {
            atom.CalculateEnergies(key);
        }

        // Stop the CI timer
        profiler->ci.stop();
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

template<typename SpecificationOptType, typename CalculatorType>
void run_and_store_transitions(SpecificationOptType specification_opt, Atom atom)
{
    if(specification_opt)
    {
        auto specification = specification_opt.value();
        CalculatorType calculon(atom, specification);
        calculon.CalculateAndPrint();
        calculon.PrintAll();
    }
}

void AmbitInterface::TransitionCalculations()
{
    // Start the timer
    auto profiler = Profiler::Instance();
    profiler->transitions.start();

    auto specification= GlobalSpecification::Instance();

    // EM types
    using EJTransitionsOpt = std::optional<TransitionsEJSpecification>;
    std::array<std::tuple<EJTransitionsOpt, MultipolarityType, int>, 8> EM_transition_types =
           {std::make_tuple(specification->transitions_do_e1 ? EJTransitionsOpt{specification->transitions_e1} : std::nullopt, MultipolarityType::E, 1),
            std::make_tuple(specification->transitions_do_e2 ? EJTransitionsOpt{specification->transitions_e2} : std::nullopt, MultipolarityType::M, 1),
            std::make_tuple(specification->transitions_do_e3 ? EJTransitionsOpt{specification->transitions_e3} : std::nullopt, MultipolarityType::E, 2),
            std::make_tuple(specification->transitions_do_e4 ? EJTransitionsOpt{specification->transitions_e4} : std::nullopt, MultipolarityType::M, 2),
            std::make_tuple(specification->transitions_do_m1 ? EJTransitionsOpt{specification->transitions_m1} : std::nullopt, MultipolarityType::E, 3),
            std::make_tuple(specification->transitions_do_m2 ? EJTransitionsOpt{specification->transitions_m2} : std::nullopt, MultipolarityType::M, 3),
            std::make_tuple(specification->transitions_do_m3 ? EJTransitionsOpt{specification->transitions_m3} : std::nullopt, MultipolarityType::E, 4),
            std::make_tuple(specification->transitions_do_m4 ? EJTransitionsOpt{specification->transitions_m4} : std::nullopt, MultipolarityType::M, 4),
           };

    std::vector<std::unique_ptr<EMCalculator>> calculators;

    for(const auto& types: EM_transition_types)
    {
        if(std::get<0>(types))
        {
            auto specification = std::get<0>(types).value();
            std::unique_ptr<EMCalculator> calculon(new EMCalculator(std::get<1>(types), std::get<2>(types), atom, specification));

            calculon->CalculateAndPrint();
            calculators.push_back(std::move(calculon));
        }
    }

    // Print all EM transitions
    for(auto& calc: calculators)
    {
        //user_input.set_prefix(std::string("Transitions/") + calc->Name());
        calc->PrintAll();
    }

    // Other operators
    //RUN_AND_STORE_TRANSITION(HFS1, HyperfineDipoleCalculator);
    //RUN_AND_STORE_TRANSITION(HFS2, HyperfineQuadrupoleCalculator);
    //RUN_AND_STORE_TRANSITION(HFI, GeneralisedHyperfineCalculator);
    //RUN_AND_STORE_TRANSITION(FS, FieldShiftCalculator);
    //RUN_AND_STORE_TRANSITION(QED, QEDCalculator);
    //RUN_AND_STORE_TRANSITION(Yukawa, YukawaCalculator);
    //RUN_AND_STORE_TRANSITION(KE, KineticEnergyCalculator);
    //RUN_AND_STORE_TRANSITION(LLIT2, LorentzInvarianceT2Calculator);
    //RUN_AND_STORE_TRANSITION(NMS, NormalMassShiftCalculator);

    // Convenient typedefs
    using HFS1Opt = std::optional<TransitionsHFS1Specification>;
    using HFS2Opt = std::optional<TransitionsHFS2Specification>;
    using HFIOpt = std::optional<TransitionsHFISpecification>;
    using FSOpt = std::optional<TransitionsFSSpecification>;
    using QEDOpt = std::optional<TransitionsQEDSpecification>;
    using YukawaOpt = std::optional<TransitionsYukawaSpecification>;
    using KEOpt = std::optional<TransitionsKESpecification>;
    using LLIT2Opt = std::optional<TransitionsLLIT2Specification>;
    using NMSOpt = std::optional<TransitionsNMSSpecification>;

    run_and_store_transitions<HFS1Opt, HyperfineDipoleCalculator>
        (specification->transitions_do_hfs1 ? HFS1Opt{specification->transitions_hfs1} : std::nullopt, atom);

    run_and_store_transitions<HFS2Opt, HyperfineQuadrupoleCalculator>
        (specification->transitions_do_hfs2 ? HFS2Opt{specification->transitions_hfs2} : std::nullopt, atom);

    run_and_store_transitions<HFIOpt, GeneralisedHyperfineCalculator>
        (specification->transitions_do_hfi ? HFIOpt{specification->transitions_hfi} : std::nullopt, atom);

    run_and_store_transitions<FSOpt, FieldShiftCalculator>
        (specification->transitions_do_fs ? FSOpt{specification->transitions_fs} : std::nullopt, atom);

    run_and_store_transitions<QEDOpt, QEDCalculator>
        (specification->transitions_do_qed ? QEDOpt{specification->transitions_qed} : std::nullopt, atom);

    run_and_store_transitions<YukawaOpt, YukawaCalculator>
        (specification->transitions_do_yukawa ? YukawaOpt{specification->transitions_yukawa} : std::nullopt, atom);

    run_and_store_transitions<KEOpt, KineticEnergyCalculator>
        (specification->transitions_do_ke ? KEOpt{specification->transitions_ke} : std::nullopt, atom);

    run_and_store_transitions<LLIT2Opt, LorentzInvarianceT2Calculator>
        (specification->transitions_do_llit2 ? LLIT2Opt{specification->transitions_llit2} : std::nullopt, atom);

    run_and_store_transitions<NMSOpt, NormalMassShiftCalculator>
        (specification->transitions_do_nms ? NMSOpt{specification->transitions_nms} : std::nullopt, atom);

    // Calculate additional Hamiltonian matrix elements
    atom.ResidualCoulomb();

    // Stop the timer
    profiler->transitions.stop();
}

#undef RUN_AND_STORE_TRANSITION
void AmbitInterface::Recombination()
{
    auto specification= GlobalSpecification::Instance();
    std::string target_file = specification->dr_target_file;
    if(target_file.empty())
        return;

    std::optional<int> target_two_j = specification->dr_target_twoj;
    if(!target_two_j)
    {   *errstream << "DR/Target/TwoJ must be set in input file." << std::endl;
        exit(1);
    }

    std::string parity_string = specification->dr_target_parity;
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
    GlobalSpecification target_specification;
    importSpecificationFile(&target_specification, target_file);
    target_specification.validateAndNormaliseSpecification();

    LevelVector target_level_vec;

    // Scope the second instance of GlobalSpecification
    {
        GlobalSpecification::ScopedOverride guard(&target_specification);

        *outstream << "\n----------------------------------------------------------";
        *outstream << "\nDR Target:" << std::endl;
        AmbitInterface target_calculator(target_specification.Z, target_specification.ID);
        target_calculator.EnergyCalculations();
        *outstream << "----------------------------------------------------------" << std::endl;

        Atom& target_atom = target_calculator.atom;

        // Get target level
        pHamiltonianID target_key = std::make_shared<HamiltonianID>(target_two_j.value(), target_parity);
        target_level_vec = target_atom.GetLevels()->GetLevels(target_key);
    }

    int target_level_index = specification->dr_target_index;
    if(target_level_index < target_level_vec.NumLevels())
    {   *errstream << "Recombination: Cannot find target level." << std::endl;
        exit(1);
    }

    target_level_vec = target_level_vec.getSubset(target_level_index);

    // Get widths of current levels
    if(specification->configuration_average)
        atom.AutoionizationConfigurationAveraged(target_level_vec);
    else if(specification->dr_use_energygrid)
        atom.AutoionizationEnergyGrid(target_level_vec);
    else
        atom.Autoionization(target_level_vec);
}

void AmbitInterface::InternalConversion()
{
    auto specification= GlobalSpecification::Instance();
    std::string source_file = specification->ic.source_filename;
    if(source_file.empty())
        return;

    std::optional<int> source_two_j = specification->ic.source_twoj;
    if(!source_two_j)
    {   *errstream << "IC/Source/TwoJ must be set in input file." << std::endl;
        exit(1);
    }

    std::string parity_string = specification->ic.source_parity;
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
    GlobalSpecification source_specification;
    importSpecificationFile(&source_specification, source_file);
    source_specification.validateAndNormaliseSpecification();

    LevelVector source_level_vec;

    // Scope the second instance of GlobalSpecification
    {
        GlobalSpecification::ScopedOverride guard(&source_specification);

        *outstream << "\n----------------------------------------------------------";
        *outstream << "\nIC Source:" << std::endl;
        AmbitInterface source_calculator(source_specification.Z, source_specification.ID);
        source_calculator.EnergyCalculations();
        *outstream << "----------------------------------------------------------" << std::endl;

        Atom& source_atom = source_calculator.atom;

        // Get source level
        pHamiltonianID source_key = std::make_shared<HamiltonianID>(source_two_j.value(), source_parity);
        source_level_vec = source_atom.GetLevels()->GetLevels(source_key);
    }

    int source_level_index = specification->ic.source_index;
    if(source_level_index < source_level_vec.NumLevels())
    {   *errstream << "IC: Cannot find source level." << std::endl;
        exit(1);
    }
    source_level_vec = source_level_vec.getSubset(source_level_index);

    if(specification->configuration_average)
        atom.InternalConversionConfigurationAveraged(source_level_vec);
    else
        atom.InternalConversion(source_level_vec);
}

void AmbitInterface::PrintHelp(const std::string& ApplicationName)
{
    *outstream << ApplicationName << std::endl;
    *outstream
        << "usage: ambit [--version] [-f] <file.input>\n"
        << "             <options> <atomic data and commands>\n"
        << "\n"
        << "Input file must be specified either with \"-f anyfile\" or \"file.input\".\n"
        << "Options may be specified in the command line or input file.\n"
        << "Consult the user guide in doc/UserGuide for usage information.\n"
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

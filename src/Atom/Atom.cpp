#include "Specification/Specification.h"
#ifdef AMBIT_USE_MPI
#include <mpi.h>
#endif
#include "Include.h"
#include "Atom.h"
#include "Basis/BasisGenerator.h"
#include "Universal/ExpLattice.h"
#include "MBPT/BruecknerDecorator.h"
#include "HartreeFock/ConfigurationParser.h"
#include "Universal/Profiler.h"

namespace Ambit
{
Atom::Atom(unsigned int atomic_number, const std::string& atom_identifier):
    Z(atomic_number), identifier(atom_identifier)
{}

Atom::~Atom(void)
{}

pCore Atom::MakeBasis(pCoreConst hf_open_core_start)
{
    // Start the timer
    auto profiler = Profiler::Instance();
    profiler->hf.start();
    auto specification = GlobalSpecification::Instance();
    bool use_read = !specification->clean_run;

    if((ProcessorRank == 0) && (!use_read || !ReadBasis()))
    {
        if(specification->lattice_exponential)
        {
            unsigned num_points = specification->lattice_num_points;
            double start_point = specification->lattice_start_point;
            double h = specification->lattice_H;
            lattice = pLattice(new ExpLattice(num_points, start_point, h));
        }
        else
        {
            unsigned num_points = specification->lattice_num_points;
            double start_point = specification->lattice_start_point;
            double end_point = specification->lattice_end_point;
            lattice = pLattice(new Lattice(num_points, start_point, end_point));
        }

        // Relativistic Hartree-Fock
        // Basis options from input
        basis_generator = std::make_shared<BasisGenerator>(lattice);
        open_core = basis_generator->GenerateHFCore(hf_open_core_start);
        hf_open = basis_generator->GetOpenHFOperator();

        // Create Basis
        orbitals = basis_generator->GenerateBasis();

        // Make closed HF operator
        hf = basis_generator->GetClosedHFOperator();

        // HartreeY operator
        hartreeY = basis_generator->GetHartreeY();

        // Nucleus
        nucleus = basis_generator->GetNucleusDecorator();

        std::string filename = identifier + ".basis";
        orbitals->Write(filename);
    }

#ifdef AMBIT_USE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
    ReadBasis();
#endif
    // Stop the timer
    profiler->hf.stop();

    return open_core;
}

bool Atom::ReadBasis()
{
    auto specification = GlobalSpecification::Instance();
    // Import lattice and all orbitals
    std::string filename = identifier + ".basis";
    FILE* fp = file_err_handler->fopen(filename.c_str(), "rb");
    if(!fp)
        return false;
    else
        file_err_handler->fclose(fp);

    pOrbitalManager modifiable_orbitals(new OrbitalManager(filename));
    lattice = modifiable_orbitals->GetLattice();

    // Generate HF operator
    // Basis and HF options from input
    basis_generator = std::make_shared<BasisGenerator>(lattice);
    hf_open = basis_generator->RecreateBasis(modifiable_orbitals);

    orbitals = modifiable_orbitals;
    hf = basis_generator->GetClosedHFOperator();

    open_core = basis_generator->GetHFCore();

    // HartreeY operator
    hartreeY = basis_generator->GetHartreeY();

    // Nucleus
    nucleus = basis_generator->GetNucleusDecorator();

    // Finally, go over the orbitals we've read and make sure they're consistent with the user input
    // First, parse the largest basis string in the user input
    std::optional<std::string> basis_size = specification->basis_size;
    std::optional<std::string> mbpt_basis = specification->mbpt_basis;
    std::optional<std::string> valence_basis = specification->basis_valence;
    // If we haven't got an explicit basis size, then get this from either the MBPT or Valence
    // Basis input options, or an empty string, in that order of preference
    std::string all_states;
    if(basis_size)
    {
        all_states = basis_size.value();
    }
    else if (mbpt_basis)
    {
        all_states = mbpt_basis.value();
    }
    else
    {
        all_states = valence_basis.value_or("");
    }

    std::vector<int> max_pqn_per_l = ConfigurationParser::ParseBasisSize(all_states);

    // Make a set to hold all of the PQN and L combinations in our basis so we can compare it against
    // what the user has specified in the input file
    auto all_orbs = orbitals->all;
    auto orb_it = all_orbs->begin();
    std::set<std::pair<int, int> > pqns;
    while(orb_it != all_orbs->end())
    {
        // Store all the PQN and l
        int pqn = orb_it->first.PQN();
        int l = orb_it->first.L();
        pqns.insert(std::make_pair(pqn, l));
        orb_it++;
    }

    // Run through all the user specified orbitals and complain if any are missing
    for(int l = 0; l < max_pqn_per_l.size(); l++)
    {
        auto key = std::make_pair(max_pqn_per_l[l], l);
        if(pqns.find(key) == pqns.end()){
            *outstream << "\nWarning: couldn't find all requested orbitals  " << filename << " basis file." << std::endl;
            *outstream << "Consider re-running AMBiT with the -c flag to recalculate the basis.\n" << std::endl;
            break;
        }
    }

    return true;
}

void Atom::GenerateBruecknerOrbitals(bool generate_sigmas)
{
    auto specification = GlobalSpecification::Instance();
    pBruecknerDecorator brueckner(new BruecknerDecorator(hf_open));
    bool use_fg = specification->mbpt_brueckner_use_lower;
    bool use_gg = specification->mbpt_brueckner_use_lower_lower;
    brueckner->IncludeLower(use_fg, use_gg);

    double sigma_start_r = specification->mbpt_brueckner_startpoint;
    double sigma_end_r   = specification->mbpt_brueckner_endpoint;
    int stride = specification->mbpt_brueckner_stride;
    brueckner->SetMatrixParameters(stride, sigma_start_r, sigma_end_r);

    pOrbitalMap orbitals_to_update = orbitals->valence;
    if(specification->mbpt_brueckner_excited)
        orbitals_to_update = orbitals->excited;

    // Attempt to read all requested kappas
    // Get max PQN for all kappas in valence orbitals
    std::map<int, int> valence_bounds;
    for(auto& pair: *orbitals_to_update)
    {
        if(valence_bounds[pair.first.Kappa()] < pair.first.PQN())
            valence_bounds[pair.first.Kappa()] = pair.first.PQN();
    }

    for(auto& kappa_pqn: valence_bounds)
        brueckner->Read(identifier, kappa_pqn.first);

    // Make new sigma potentials if they haven't been read (slowly)
    if(generate_sigmas)
    {
        std::optional<std::string> fermi_orbitals = specification->mbpt_energy_denom_orbitals;
        for(auto& kappa_maxpqn: valence_bounds)
        {   brueckner->CalculateSigma(kappa_maxpqn.first, orbitals, hartreeY, fermi_orbitals);
            brueckner->Write(identifier, kappa_maxpqn.first);
        }
    }

    // Get scalings or energies
    unsigned int scaling_length = specification->mbpt_brueckner_scaling.size();

    // Run this one unless MBPT/Brueckner/EnergyScaling option is used
    if(!specification->mbpt_brueckner_scaling.empty())
    {
        for(int i = 0; i < scaling_length-1; i+=2)
        {
            int kappa = specification->mbpt_brueckner_energy_scaling[i];
            double scale = specification->mbpt_brueckner_scaling[i+1];
            brueckner->SetSigmaScaling(kappa, scale);
        }
    }

    basis_generator->CreateBruecknerOrbitals(brueckner);

    // Replace hf operator for rest of calculation
    hf_open = basis_generator->GetOpenHFOperator();
    hf = basis_generator->GetClosedHFOperator();

    *outstream << "Brueckner orbitals:\n";
    orbitals_to_update->Print();
    *outstream << std::endl;
}
}

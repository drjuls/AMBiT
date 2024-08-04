#ifdef AMBIT_USE_MPI
#include <mpi.h>
#endif
#include "Include.h"
#include "Atom.h"
#include "Basis/BasisGenerator.h"
#include "Basis/BSplineBasis.h"
#include "Universal/ExpLattice.h"
#include "MBPT/BruecknerDecorator.h"
#include "HartreeFock/ConfigurationParser.h"

namespace Ambit
{
Atom::Atom(const MultirunOptions userInput, unsigned int atomic_number, const std::string& atom_identifier):
    user_input(userInput), Z(atomic_number), identifier(atom_identifier)
{}

Atom::~Atom(void)
{}

pCore Atom::MakeBasis(pCoreConst hf_open_core_start)
{
    bool use_read = true;
    if(user_input.search(2, "--clean", "-c"))
        use_read = false;

    if((ProcessorRank == 0) && (!use_read || !ReadBasis()))
    {
        int N = user_input("HF/N", -1);  // Number of electrons
        if(N == -1 || N > Z)
        {   *errstream << "Must specify N (number of electrons) with Z >= N" << std::endl;
            exit(1);
        }

        if(user_input.search("HF/--read-grasp0"))
        {   // Read lattice and core and basis orbitals
    //        ReadGraspMCDF("MCDF.DAT");
        }
        else
        {   // Lattice parameters
            if(user_input.search("Lattice/--exp-lattice"))
            {   int num_points = user_input("Lattice/NumPoints", 300);
                double first_point = user_input("Lattice/StartPoint", 1.e-5);
                double h = user_input("Lattice/H", 0.05);
                lattice = pLattice(new ExpLattice(num_points, first_point, h));
            }
            else
            {   int num_points = user_input("Lattice/NumPoints", 1000);
                double first_point = user_input("Lattice/StartPoint", 1.e-6);
                double lattice_size = user_input("Lattice/EndPoint", 50.);
                lattice = pLattice(new Lattice(num_points, first_point, lattice_size));
            }
        }

        // Relativistic Hartree-Fock
        basis_generator = std::make_shared<BasisGenerator>(lattice, user_input);
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

    return open_core;
}

bool Atom::ReadBasis()
{
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
    basis_generator = std::make_shared<BasisGenerator>(lattice, user_input);
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
    std::string all_states = user_input("Basis/BasisSize", "");
    if(all_states.empty())
        all_states = user_input("MBPT/Basis", "");
    if(all_states.empty())
        all_states = user_input("Basis/ValenceBasis", "");

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
    pBruecknerDecorator brueckner(new BruecknerDecorator(hf_open));
    bool use_fg = user_input.search("MBPT/Brueckner/--use-lower");
    bool use_gg = user_input.search("MBPT/Brueckner/--use-lower-lower");
    brueckner->IncludeLower(use_fg, use_gg);

    double sigma_start_r = user_input("MBPT/Brueckner/StartPoint", 4.35e-5);
    double sigma_end_r   = user_input("MBPT/Brueckner/EndPoint", 8.0);
    int stride = user_input("MBPT/Brueckner/Stride", 4);
    brueckner->SetMatrixParameters(stride, sigma_start_r, sigma_end_r);

    pOrbitalMap orbitals_to_update = orbitals->valence;
    if(user_input.search("MBPT/Brueckner/--excited"))
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
        std::string fermi_orbitals = user_input("MBPT/EnergyDenomOrbitals", "");
        for(auto& kappa_maxpqn: valence_bounds)
        {   brueckner->CalculateSigma(kappa_maxpqn.first, orbitals, hartreeY, fermi_orbitals);
            brueckner->Write(identifier, kappa_maxpqn.first);
        }
    }

    // Get scalings or energies
    unsigned int scaling_length = user_input.vector_variable_size("MBPT/Brueckner/Scaling");

    // Run this one unless MBPT/Brueckner/EnergyScaling option is used
    if(scaling_length)
    {
        for(int i = 0; i < scaling_length-1; i+=2)
        {
            int kappa = user_input("MBPT/Brueckner/Scaling", 0, i);
            double scale = user_input("MBPT/Brueckner/Scaling", 0.0, i+1);

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

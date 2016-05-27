#ifdef AMBIT_USE_MPI
#include <mpi.h>
#endif
#include "Include.h"
#include "Atom.h"
#include "Basis/BasisGenerator.h"
#include "Universal/ExpLattice.h"
#include "MBPT/BruecknerDecorator.h"
#include "HartreeFock/HartreeFocker.h"
#include "HartreeFock/ConfigurationParser.h"

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

    if(!use_read || !ReadBasis())
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
        BasisGenerator basis_generator(lattice, user_input);
        open_core = basis_generator.GenerateHFCore(hf_open_core_start);
        hf_open = basis_generator.GetOpenHFOperator();

        // Create Basis
        orbitals = basis_generator.GenerateBasis();

        // Make closed HF operator
        hf = basis_generator.GetClosedHFOperator();

        // HartreeY operator
        hartreeY = basis_generator.GetHartreeY();

        // Nucleus
        nucleus = basis_generator.GetNucleusDecorator();

        std::string filename = identifier + ".basis";
        orbitals->Write(filename);
    }

    return open_core;
}

bool Atom::ReadBasis()
{
    // Import lattice and all orbitals
    std::string filename = identifier + ".basis";
    FILE* fp = fopen(filename.c_str(), "rb");
    if(!fp)
        return false;
    else
        fclose(fp);

    pOrbitalManager modifiable_orbitals(new OrbitalManager(filename));
    lattice = modifiable_orbitals->GetLattice();

    // Generate HF operator
    BasisGenerator basis_generator(lattice, user_input);
    hf_open = basis_generator.RecreateBasis(modifiable_orbitals);

    orbitals = modifiable_orbitals;
    hf = basis_generator.GetClosedHFOperator();

    open_core = basis_generator.GetHFCore();

    // HartreeY operator
    hartreeY = basis_generator.GetHartreeY();

    // Nucleus
    nucleus = basis_generator.GetNucleusDecorator();

    return true;
}

void Atom::GenerateBruecknerOrbitals(bool generate_sigmas)
{
    pBruecknerDecorator brueckner(new BruecknerDecorator(hf));
    bool use_fg = user_input.search("MBPT/Brueckner/--use-lower");
    bool use_gg = user_input.search("MBPT/Brueckner/--use-lower-lower");
    brueckner->IncludeLower(use_fg, use_gg);

    double sigma_start_r = user_input("MBPT/Brueckner/StartPoint", 4.35e-5);
    double sigma_end_r   = user_input("MBPT/Brueckner/EndPoint", 8.0);
    int stride = user_input("MBPT/Brueckner/Stride", 4);
    brueckner->SetMatrixParameters(stride, sigma_start_r, sigma_end_r);

    // Attempt to read all requested kappas
    std::set<int> valence_kappas;
    for(auto& pair: *orbitals->valence)
        valence_kappas.insert(pair.first.Kappa());

    for(int kappa: valence_kappas)
        brueckner->Read(identifier, kappa);

    // Replace hf operator for rest of calculation
    hf = brueckner;

    // Make new sigma potentials if they haven't been read (slowly)
    if(generate_sigmas)
    {
        std::string fermi_orbitals = user_input("MBPT/EnergyDenomOrbitals", "");
        for(int kappa: valence_kappas)
        {   brueckner->CalculateSigma(kappa, orbitals, hartreeY, fermi_orbitals);
            brueckner->Write(identifier, kappa);
        }
    }

    // And finally change all our valence orbitals to Brueckner orbitals
    pIntegrator integrator(new SimpsonsIntegrator(lattice));
    pODESolver ode_solver(new AdamsSolver(integrator));
    HartreeFocker hartree_focker(ode_solver);

    // Get scalings or energies
    unsigned int scaling_length = user_input.vector_variable_size("MBPT/Brueckner/Scaling");
    unsigned int energy_scaling_length = user_input.vector_variable_size("MBPT/Brueckner/EnergyScaling");

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
    else if(energy_scaling_length)
    {   pOrbital brueckner_orbital;

        *outstream << "Brueckner scaling:" << std::endl;

        // Get orbitals to scale and scale them
        for(int i = 0; i < energy_scaling_length-1; i+=2)
        {
            std::string orbital_info = user_input("MBPT/Brueckner/EnergyScaling", "", i);
            OrbitalInfo info(ConfigurationParser::ParseOrbital<OrbitalInfo>(orbital_info));
            double requested_energy = user_input("MBPT/Brueckner/EnergyScaling", 0.0, i+1);

            pOrbital orbital = orbitals->valence->GetState(info);
            if(orbital)
            {
                double initial_energy = orbital->Energy();
                double requested_energy_change = requested_energy - initial_energy;

                // Make a copy of the old orbital
                brueckner_orbital.reset(new Orbital(*orbital));

                double deltaE;
                double scale = 1.;
                do
                {   brueckner->SetSigmaScaling(info.Kappa(), scale);

                    // Iterate
                    hartree_focker.ConvergeOrbitalAndExchange(brueckner_orbital, brueckner, &HartreeFocker::IterateOrbital, 1.e-7 * fabs(requested_energy));

                    // Rescale
                    deltaE = brueckner_orbital->Energy() - initial_energy;
                    scale = requested_energy_change/deltaE * scale;
                } while(!std::isnan(scale) && fabs(scale) < 2. && fabs((deltaE - requested_energy_change)/initial_energy) > 1.e-6);

                *outstream << "  kappa = " << info.Kappa() << ": scaling = " << std::setprecision(6) << scale << std::endl;
            }
            else
                *errstream << "Brueckner: Orbital " << info.Name() << " not found in valence set.";
        }

        *outstream << std::endl;
    }

    pOrbital brueckner_orbital;
    for(auto& pair: *orbitals->valence)
    {
        // Make a copy of the old orbital
        brueckner_orbital.reset(new Orbital(*pair.second));

        // Iterate
        hartree_focker.ConvergeOrbitalAndExchange(brueckner_orbital, brueckner, &HartreeFocker::IterateOrbital, hartree_focker.EnergyTolerance);

        // Copy back to orbital manager
        *pair.second = *brueckner_orbital;
    }

    *outstream << "Brueckner orbitals:\n";
    orbitals->valence->Print();
    *outstream << std::endl;
}

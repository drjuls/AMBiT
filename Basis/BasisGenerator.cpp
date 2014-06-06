#include "BasisGenerator.h"
#include "Include.h"
#include "HartreeFock/ConfigurationParser.h"
#include "HartreeFock/OpIntegrator.h"
#include "HartreeFock/HartreeFocker.h"
#include "HartreeFock/NucleusDecorator.h"
#include "HartreeFock/MassShiftDecorator.h"
#include "HartreeFock/NonRelativisticSMSOperator.h"

BasisGenerator::BasisGenerator(pLattice lat, MultirunOptions& userInput, pPhysicalConstant physical_constant):
    hf(pHFOperator()), lattice(lat), user_input(userInput), physical_constant(physical_constant), open_core(nullptr)
{
    orbitals = pOrbitalManager(new OrbitalManager(lattice));
}

BasisGenerator::~BasisGenerator()
{}

void BasisGenerator::InitialiseHF(pHFOperator& undressed_hf)
{
    unsigned int Z = user_input("Z", 0);

    int Charge = user_input("HF/Charge", -1);
    if(Charge < 0)
    {   int N = user_input("HF/N", -1);
        if(Z >= N && N >= 0)
            Charge = Z - N;
    }

    //TODO: Error message if Charge or N is missing or incorrect.
    std::string config(user_input("HF/Configuration", ""));

    // Get orbitals and occupancies
    std::string open_shell_string;
    size_t colon_pos = config.find(':');
    if(colon_pos == std::string::npos)
        open_shell_string = config;
    else
    {   open_shell_string = config;
        open_shell_string.erase(colon_pos, 1);
    }

    OccupationMap open_shell_occupations = ConfigurationParser::ParseFractionalConfiguration(open_shell_string);

    // Set open_core occupancies
    open_core->SetOccupancies(open_shell_occupations);
    if(open_core->NumElectrons() != Z - Charge)
    {   *errstream << "Core::BuildFirstApproximation: Incorrect electron count in configuration." << std::endl;
        exit(1);
    }

    pOPIntegrator integrator(new SimpsonsIntegrator(lattice));
    pODESolver ode_solver(new AdamsSolver(integrator));
    pCoulombOperator coulomb(new CoulombOperator(lattice, ode_solver));

    if(physical_constant == nullptr)
    {
        physical_constant = pPhysicalConstant(new PhysicalConstant());
        double alpha_variation = user_input("AlphaSquaredVariation", 0.0);
        if(alpha_variation)
            physical_constant->SetAlphaSquaredIncreaseRatio(alpha_variation);
    }

    undressed_hf = pHFOperator(new HFOperator(Z, open_core, physical_constant, integrator, coulomb));
    hf = undressed_hf;

    // Add nuclear potential
    double nuclear_radius = user_input("NuclearRadius", 0.0);
    if(nuclear_radius)
    {
        pNucleusDecorator nucleus(new NucleusDecorator(hf));
        double nuclear_thickness = user_input("NuclearThickness", 2.3);
        nucleus->SetFermiParameters(nuclear_radius, nuclear_thickness);
        nucleus->SetCore(open_core);
        hf = nucleus;
        undressed_hf = hf;
    }

    // Hartree operator
    hartreeY = pHartreeY(new HartreeY(integrator, coulomb));

    // Add additional operators
    double NuclearInverseMass = user_input("NuclearInverseMass", 0.0);
    if(NuclearInverseMass)
    {   pMassShiftDecorator sms_op(new MassShiftDecorator(hf));
        sms_op->SetInverseMass(NuclearInverseMass);
        sms_op->SetCore(open_core);
        hf = sms_op;

        pHartreeY dressed(new NonRelativisticSMSOperator(hartreeY));
        hartreeY = dressed;
    }
}

void BasisGenerator::SetOrbitalMaps()
{
    // Get closed-shell core orbitals and occupancies
    std::string config(user_input("HF/Configuration", ""));
    std::string closed_shell_string;
    size_t colon_pos = config.find(':');
    if(colon_pos == std::string::npos)
        closed_shell_string = config;
    else
        closed_shell_string = config.substr(0, colon_pos);
    OccupationMap closed_shell_occupations = ConfigurationParser::ParseFractionalConfiguration(closed_shell_string);

    // Make closed shell core. Ensure that all shells are completely filled.
    for(OccupationMap::iterator it = closed_shell_occupations.begin(); it != closed_shell_occupations.end(); it++)
        it->second = 2. * abs(it->first.Kappa());

    // Create closed core with empty pointers for all occupied orbitals
    pCore closed_core = pCore(new Core(lattice));
    closed_core->SetOccupancies(closed_shell_occupations);

    // Transfer from all to closed core
    OrbitalMap& all = *orbitals->all;
    for(auto core_occupation: closed_shell_occupations)
    {
        closed_core->AddState(all.GetState(core_occupation.first));
    }
    orbitals->core = closed_core;

    // TODO: implement hole states.
    orbitals->hole = pOrbitalMap(new OrbitalMap(lattice));
    orbitals->deep = orbitals->core;

    // Transfer from all to excited states
    std::string valence_states = user_input("Basis/ValenceBasis", "");
    std::vector<int> max_pqn_per_l = ConfigurationParser::ParseBasisSize(valence_states);

    pOrbitalMap particle = pOrbitalMap(new OrbitalMap(lattice));
    for(auto orbital: all)
    {
        if(orbital.first.L() < max_pqn_per_l.size()
           && orbital.first.PQN() <= max_pqn_per_l[orbital.first.L()]
           && closed_core->GetState(orbital.first) == nullptr)
        {
            particle->AddState(orbital.second);
        }
    }
    orbitals->valence = particle;
    orbitals->particle = particle;

    // high (virtual) states.
    std::string virtual_states = user_input("MBPT/Basis", "");
    if(virtual_states.size())
    {
        max_pqn_per_l = ConfigurationParser::ParseBasisSize(virtual_states);

        pOrbitalMap excited = pOrbitalMap(new OrbitalMap(lattice));
        pOrbitalMap high = pOrbitalMap(new OrbitalMap(lattice));
        for(auto orbital: all)
        {
            if(orbital.first.L() < max_pqn_per_l.size()
               && orbital.first.PQN() <= max_pqn_per_l[orbital.first.L()]
               && closed_core->GetState(orbital.first) == nullptr)
            {
                excited->AddState(orbital.second);
                if(particle->GetState(orbital.first) == nullptr)
                    high->AddState(orbital.second);
            }
        }
        orbitals->excited = excited;
        orbitals->high = high;
    }
    else
    {   orbitals->high = pOrbitalMap(new OrbitalMap(lattice));
        orbitals->excited = orbitals->particle;
    }
}

pCore BasisGenerator::GenerateHFCore(pCoreConst open_shell_core)
{
    open_core = pCore(new Core(lattice));
    hf = nullptr;
    hartreeY = nullptr;

    if(open_shell_core)
    {   // Copy, use same lattice
        open_core->Copy(*open_shell_core);
        lattice = open_core->GetLattice();
        orbitals = pOrbitalManager(new OrbitalManager(lattice));
    }

    pHFOperator undressed_hf;
    InitialiseHF(undressed_hf);

    // Create Hartree-Fock solver; define integrators.
    pOPIntegrator integrator(new SimpsonsIntegrator(lattice));
    pODESolver ode_solver(new AdamsSolver(integrator));
    HartreeFocker HF_Solver(ode_solver);

    // TODO: Check occupancies match
    if(!open_shell_core)
    {   HF_Solver.StartCore(open_core, undressed_hf);
        HF_Solver.SolveCore(open_core, undressed_hf);
    }

    HF_Solver.SolveCore(open_core, hf);

    return open_core;
}

pHFOperator BasisGenerator::RecreateBasis(pOrbitalManager orbital_manager)
{
    lattice = orbital_manager->GetLattice();
    open_core = pCore(new Core(lattice));

    pHFOperator undressed_hf;
    InitialiseHF(undressed_hf);

    // Copy orbitals from orbital_manager to open_core
    for(auto pair: *open_core)
    {
        pOrbital state = orbital_manager->all->GetState(pair.first);
        if(state == nullptr)
        {   *errstream << "BasisGenerator::CreateHFOperator(): orbital " << pair.first.Name() << " not found." << std::endl;
            exit(1);
        }
        *pair.second = *state;
    }

    hf->SetCore(open_core);

    // Modify orbital manager maps according to input file
    orbitals = orbital_manager;
    SetOrbitalMaps();

    return hf;
}

pOrbitalManagerConst BasisGenerator::GenerateBasis()
{
    // Generate excited states
    std::string all_states = user_input("Basis/BasisSize", "");
    if(all_states.empty())
        all_states = user_input("MBPT/Basis", "");
    if(all_states.empty())
        all_states = user_input("Basis/ValenceBasis", "");

    std::vector<int> max_pqn_per_l = ConfigurationParser::ParseBasisSize(all_states);

    pOrbitalMap excited;

    if(user_input.search("Basis/--bspline-basis"))
    {   excited = GenerateBSplines(max_pqn_per_l);
    }
    else
        excited = pOrbitalMap(new OrbitalMap(lattice));   // Just to stop seg-faults

    // Place all orbitals in orbitals->all.
    // Finally create orbitals->all and the state index
    orbitals->all = pOrbitalMap(new OrbitalMap(lattice));
    orbitals->all->AddStates(*open_core);
    orbitals->all->AddStates(*excited);

    orbitals->MakeStateIndexes();

    // Organise orbitals
    SetOrbitalMaps();

    if(DebugOptions.OutputHFExcited())
    {   OrbitalInfo max_i(-1, 1), max_j(-1, 1);
        double orth = TestOrthogonality(max_i, max_j);
        *outstream << "<" << max_i.Name() << " | " << max_j.Name() << "> = " << orth << std::endl;
    }

    return orbitals;
}

void BasisGenerator::Orthogonalise(pOrbital current) const
{
    pOPIntegrator integrator(hf->GetOPIntegrator());
    current->ReNormalise(integrator);

    // Orthogonalise to core
    if(orbitals->core)
    {
        auto it = orbitals->core->begin();
        while(it != orbitals->core->end())
        {
            pOrbitalConst other = it->second;
            if((other->Kappa() == current->Kappa()) && (other->PQN() < current->PQN()))
            {
                double S = integrator->GetInnerProduct(*other, *current);
                (*current) -= (*other) * S;

                current->ReNormalise(integrator);
            }
            it++;
        }
    }

    // Orthogonalise to other excited states.
    if(orbitals->excited)
    {
        auto it = orbitals->excited->begin();
        while(it != orbitals->excited->end())
        {
            pOrbitalConst other = it->second;
            if((other->Kappa() == current->Kappa()) && (other->PQN() < current->PQN()))
            {
                double S = integrator->GetInnerProduct(*other, *current);
                (*current) -= (*other) * S;

                current->ReNormalise(integrator);
            }
            it++;
        }
    }

    current->SetEnergy(hf->GetMatrixElement(*current, *current));
}

double BasisGenerator::TestOrthogonality(OrbitalInfo& max_i, OrbitalInfo& max_j) const
{
    double max_orth = 0.;
    pOPIntegrator integrator = hf->GetOPIntegrator();

    pOrbitalMap all_states;
    if(orbitals->all)
        all_states = orbitals->all;
    else if(open_core)
        all_states = open_core;
    else
        return max_orth;

    auto it = all_states->begin();
    while(it != all_states->end())
    {
        auto jt = all_states->begin();
        while(jt != all_states->end() && (it->first != jt->first))
        {
            if(it->first.Kappa() == jt->first.Kappa())
            {
                double orth = fabs(integrator->GetInnerProduct(*it->second, *jt->second));
                if(orth > max_orth)
                {   max_orth = orth;
                    max_i = it->first;
                    max_j = jt->first;
                }
            }
            jt++;
        }

        it++;
    }

    return max_orth;
}

#include "BasisGenerator.h"
#include "Include.h"
#include "HartreeFock/ConfigurationParser.h"
#include "HartreeFock/OpIntegrator.h"
#include "HartreeFock/HartreeFocker.h"
#include "HartreeFock/NucleusDecorator.h"
#include "HartreeFock/MassShiftDecorator.h"

BasisGenerator::BasisGenerator(pLattice lat, MultirunOptions& userInput):
    hf(pHFOperator()), lattice(lat), user_input(userInput), open_core(nullptr)
{
    orbitals = pOrbitalManager(new OrbitalManager());
}

BasisGenerator::~BasisGenerator()
{}

pCore BasisGenerator::GenerateHFCore(pCoreConst open_shell_core)
{
    open_core = pCore(new Core(lattice));
    if(open_shell_core)
        open_core->Copy(*open_shell_core);

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

    open_core->SetOccupancies(open_shell_occupations);
    if(open_core->NumElectrons() != Z - Charge)
    {   *errstream << "Core::BuildFirstApproximation: Incorrect electron count in configuration." << std::endl;
        exit(1);
    }

    pOPIntegrator integrator(new SimpsonsIntegrator(lattice));
    pODESolver ode_solver(new AdamsSolver(integrator));
    pCoulombOperator coulomb(new CoulombOperator(lattice, ode_solver));

    hf.reset(new HFOperator(Z, open_core, integrator, coulomb));

    // Add nuclear potential
    double nuclear_radius = user_input("NuclearRadius", 0.0);
    if(nuclear_radius)
    {
        pNucleusDecorator nucleus(new NucleusDecorator(hf));
        double nuclear_thickness = user_input("NuclearThickness", 2.3);
        nucleus->SetFermiParameters(nuclear_radius, nuclear_thickness);
        nucleus->SetCore(open_core);
        hf = nucleus;
    }

    HartreeFocker HF_Solver(ode_solver);

    // TODO: Check occupancies match
    if(!open_shell_core)
    {   HF_Solver.StartCore(open_core, hf);
        HF_Solver.SolveCore(open_core, hf);
    }

    // Add additional operators
    double NuclearInverseMass = user_input("NuclearInverseMass", 0.0);
    if(NuclearInverseMass)
    {   pMassShiftDecorator sms_op(new MassShiftDecorator(hf));
        sms_op->SetInverseMass(NuclearInverseMass);
        sms_op->SetCore(open_core);
        hf = sms_op;
    }

    HF_Solver.SolveCore(open_core, hf);

    if(DebugOptions.LogHFIterations())
    {   OrbitalInfo max_i(-1, 1), max_j(-1, 1);
        double orth = TestOrthogonality(max_i, max_j);
        *logstream << "Core Orthogonality test: "
                   << "<" << max_i.Name() << " | " << max_j.Name() << "> = " << orth << std::endl;
    }

    return open_core;
}

pOrbitalManager BasisGenerator::GenerateBasis()
{
    // Core states first
    // Get closed shell orbitals and occupancies
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

    pCore closed_core = pCore(new Core(open_core->Copy()));
    closed_core->SetOccupancies(closed_shell_occupations);
    orbitals->core = closed_core;

    // TODO: implement hole states.
    orbitals->hole = pOrbitalMap(new OrbitalMap(lattice));
    orbitals->deep = orbitals->core;

    // Excited states
    std::string valence_states = user_input("Basis/ValenceBasis", "");
    std::vector<int> max_pqn_per_l = ConfigurationParser::ParseBasisSize(valence_states);

    if(user_input.search("Basis/--bspline-basis"))
    {   orbitals->excited = GenerateBSplines(max_pqn_per_l);
    }
    else
        orbitals->excited = pOrbitalMap(new OrbitalMap(lattice));   // Just to stop seg-faults

    if(DebugOptions.OutputHFExcited())
    {   OrbitalInfo max_i(-1, 1), max_j(-1, 1);
        double orth = TestOrthogonality(max_i, max_j);
        *outstream << "<" << max_i.Name() << " | " << max_j.Name() << "> = " << orth << std::endl;
    }

    // TODO: implement high (virtual) states.
    orbitals->high = pOrbitalMap(new OrbitalMap(lattice));
    orbitals->particle = orbitals->excited;
    orbitals->valence = orbitals->excited;

    // Finally create orbitals->all and the state index
    orbitals->all = pOrbitalMap(new OrbitalMap(lattice));
    orbitals->all->AddStates(*orbitals->core);
    orbitals->all->AddStates(*orbitals->excited);

    orbitals->MakeStateIndexes();

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

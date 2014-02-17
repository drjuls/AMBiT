#include "BasisGenerator.h"
#include "Include.h"
#include "HartreeFock/OpIntegrator.h"
#include "HartreeFock/HartreeFocker.h"
#include "HartreeFock/NucleusDecorator.h"
#include "HartreeFock/MassShiftDecorator.h"

BasisGenerator::BasisGenerator(pLattice lat, MultirunOptions& userInput):
    hf(pHFOperator()), lattice(lat), user_input(userInput)
{}

BasisGenerator::~BasisGenerator()
{}

pCore BasisGenerator::GenerateCore(pCoreConst open_shell_core)
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
    std::string closed_shell_string;

    size_t colon_pos = config.find(':');
    if(colon_pos == std::string::npos)
        closed_shell_string = open_shell_string = config;
    else
    {   closed_shell_string = config.substr(0, colon_pos);
        open_shell_string = config;
        open_shell_string.erase(colon_pos, 1);
    }

    OccupationMap open_shell_occupations = ConfigurationParser::ParseFractionalConfiguration(open_shell_string);
    OccupationMap closed_shell_occupations = ConfigurationParser::ParseFractionalConfiguration(closed_shell_string);

    open_core->SetOccupancies(open_shell_occupations);
    if(open_core->NumElectrons() != Z - Charge)
    {   *errstream << "Core::BuildFirstApproximation: Incorrect electron count in configuration." << std::endl;
        exit(1);
    }

    pOPIntegrator integrator(new SimpsonsIntegrator(lattice));
    pODESolver ode_solver(new AdamsSolver(lattice));
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

    closed_core = pCore(new Core(open_core->Copy()));
    closed_core->SetOccupancies(closed_shell_occupations);

    return open_core;
}

pStateManager BasisGenerator::GenerateBasis()
{
    excited = pStateManager(new StateManager(lattice));

    std::string valence_states = user_input("Basis/ValenceBasis", "");
    std::vector<int> max_pqn_per_l = ConfigurationParser::ParseBasisSize(valence_states);

    if(user_input.search("Basis/--bspline-basis"))
    {   GenerateBSplines(max_pqn_per_l);
    }

    if(DebugOptions.OutputHFExcited())
    {   OrbitalInfo max_i(-1, 1), max_j(-1, 1);
        double orth = TestOrthogonality(max_i, max_j);
        *outstream << "<" << max_i.Name() << " | " << max_j.Name() << "> = " << orth << std::endl;
    }

    return excited;
}

void BasisGenerator::Orthogonalise(pOrbital current) const
{
    current->ReNormalise(lattice);
    pOPIntegrator integrator(hf->GetOPIntegrator());

    // Orthogonalise to core
    if(closed_core)
    {
        ConstStateIterator it = closed_core->GetConstStateIterator();
        while(!it.AtEnd())
        {
            pOrbitalConst other = it.GetState();
            if((other->Kappa() == current->Kappa()) && (other->PQN() < current->PQN()))
            {
                double S = integrator->GetInnerProduct(*other, *current);
                (*current) -= (*other) * S;

                current->ReNormalise(lattice);
            }
            it.Next();
        }
    }

    // Orthogonalise to other excited states.
    if(excited)
    {
        ConstStateIterator it = excited->GetConstStateIterator();
        while(!it.AtEnd())
        {
            pOrbitalConst other = it.GetState();
            if((other->Kappa() == current->Kappa()) && (other->PQN() < current->PQN()))
            {
                double S = integrator->GetInnerProduct(*other, *current);
                (*current) -= (*other) * S;

                current->ReNormalise(lattice);
            }
            it.Next();
        }
    }

    current->SetEnergy(hf->GetMatrixElement(*current, *current));
}

double BasisGenerator::TestOrthogonality(OrbitalInfo& max_i, OrbitalInfo& max_j) const
{
    double max_orth = 0.;
    pOPIntegrator integrator = hf->GetOPIntegrator();

    StateManager all_states(lattice);
    if(excited)
    {   all_states = *excited;
        if(closed_core)
            all_states.AddStates(*closed_core);
        else if(open_core)
            all_states.AddStates(*open_core);
    }
    else if(open_core)
        all_states = *open_core;
    else if(closed_core)
        all_states = *closed_core;

    ConstStateIterator it = all_states.GetConstStateIterator();
    ConstStateIterator jt = all_states.GetConstStateIterator();

    it.First();
    while(!it.AtEnd())
    {
        jt.First();
        while(!jt.AtEnd() && (it.GetOrbitalInfo() != jt.GetOrbitalInfo()))
        {
            if(it.GetOrbitalInfo() != jt.GetOrbitalInfo())
            {
                double orth = fabs(integrator->GetInnerProduct(*it.GetState(), *jt.GetState()));
                if(orth > max_orth)
                {   max_orth = orth;
                    max_i = it.GetOrbitalInfo();
                    max_j = jt.GetOrbitalInfo();
                }
            }
            jt.Next();
        }

        it.Next();
    }

    return max_orth;
}

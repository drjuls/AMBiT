#include "CoreGenerator.h"
#include "Include.h"
#include "Core.h"
#include "OpIntegrator.h"
#include "HartreeFocker.h"
#include "NucleusDecorator.h"
#include "MassShiftDecorator.h"

CoreGenerator::CoreGenerator(pLattice lat):
    hf(pHFOperator()), lattice(lat)
{}

CoreGenerator::~CoreGenerator()
{}

void CoreGenerator::GenerateCore(MultirunOptions& userInput, Core* open_shell_core, Core* closed_shell_core)
{
    unsigned int Z = userInput("Z", 0);

    int Charge = userInput("HF/Charge", -1);
    if(Charge < 0)
    {   int N = userInput("HF/N", -1);
        if(Z >= N && N >= 0)
            Charge = Z - N;
    }

    //TODO: Error message if Charge or N is missing or incorrect.
    std::string config(userInput("HF/Configuration", ""));

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

    open_shell_core->SetOccupancies(open_shell_occupations);
    if(open_shell_core->NumElectrons() != Z - Charge)
    {   *errstream << "Core::BuildFirstApproximation: Incorrect electron count in configuration." << std::endl;
        exit(1);
    }

    pOPIntegrator integrator(new SimpsonsIntegrator(lattice));
    pODESolver ode_solver(new AdamsSolver(lattice));
    pCoulombOperator coulomb(new CoulombOperator(lattice, ode_solver));

    hf.reset(new HFOperator(Z, open_shell_core, integrator, coulomb));

    // Add nuclear potential
    double nuclear_radius = userInput("NuclearRadius", 0.0);
    if(nuclear_radius)
    {
        pNucleusDecorator nucleus(new NucleusDecorator(hf));
        double nuclear_thickness = userInput("NuclearThickness", 2.3);
        nucleus->SetFermiParameters(nuclear_radius, nuclear_thickness);
        nucleus->SetCore(open_shell_core);
        hf = nucleus;
    }

    HartreeFocker HF_Solver(ode_solver);
    HF_Solver.StartCore(open_shell_core, hf);
    HF_Solver.SolveCore(open_shell_core, hf);

    // Add additional operators
    double NuclearInverseMass = userInput("NuclearInverseMass", 0.0);
    if(NuclearInverseMass)
    {   pMassShiftDecorator sms_op(new MassShiftDecorator(hf));
        sms_op->SetInverseMass(NuclearInverseMass);
        sms_op->SetCore(open_shell_core);
        hf = sms_op;
    }

    HF_Solver.SolveCore(open_shell_core, hf);

    *closed_shell_core = *open_shell_core;
    closed_shell_core->SetOccupancies(closed_shell_occupations);
}

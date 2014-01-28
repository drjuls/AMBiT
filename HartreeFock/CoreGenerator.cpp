#include "CoreGenerator.h"
#include "Include.h"
#include "OpIntegrator.h"
#include "HartreeFocker.h"
#include "NucleusDecorator.h"
#include "MassShiftDecorator.h"

CoreGenerator::CoreGenerator(pLattice lat):
    hf(pHFOperator()), lattice(lat)
{}

CoreGenerator::~CoreGenerator()
{}

Core* CoreGenerator::GenerateCore(MultirunOptions& userInput)
{
    unsigned int Z = userInput("Z", 0);

    int Charge = userInput("HF/Charge", -1);
    if(Charge < 0)
    {   int N = userInput("HF/N", -1);
        if(Z >= N && N >= 0)
            Charge = Z - N;
    }

    //TODO: Error message if Charge or N is missing or incorrect.
    Core* core = new Core(lattice, Z, Charge);

    std::string config(userInput("HF/Configuration", ""));
    core->Initialise(config);

    pOPIntegrator integrator(new SimpsonsIntegrator(lattice));
    pODESolver ode_solver(new AdamsSolver(lattice));
    pCoulombOperator coulomb(new CoulombOperator(lattice, ode_solver));

    hf.reset(new HFOperator(Z, core, integrator, coulomb));

    // Add nuclear potential
    double nuclear_radius = userInput("NuclearRadius", 0.0);
    if(nuclear_radius)
    {
        pNucleusDecorator nucleus(new NucleusDecorator(hf));
        double nuclear_thickness = userInput("NuclearThickness", 2.3);
        nucleus->SetFermiParameters(nuclear_radius, nuclear_thickness);
        nucleus->SetCore(core);
        hf = nucleus;
    }

    HartreeFocker HF_Solver(ode_solver);
    HF_Solver.SolveCore(core, hf);

    // Add additional operators
    double NuclearInverseMass = userInput("NuclearInverseMass", 0.0);
    if(NuclearInverseMass)
    {   pMassShiftDecorator sms_op(new MassShiftDecorator(hf));
        sms_op->SetInverseMass(NuclearInverseMass);
        sms_op->SetCore(core);
        hf = sms_op;
    }

    HF_Solver.SolveCore(core, hf);

    return core;
}

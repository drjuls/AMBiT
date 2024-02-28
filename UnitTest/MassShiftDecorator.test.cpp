#include "ExternalField/SpecificMassShiftDecorator.h"
#include "gtest/gtest.h"
#include "HartreeFock/Core.h"
#include "Include.h"
#include "HartreeFock/ODESolver.h"
#include "Universal/PhysicalConstant.h"
#include "HartreeFock/HartreeFocker.h"
#include "HartreeFock/ConfigurationParser.h"

using namespace Ambit;

TEST(MassShiftDecoratorTester, CaII)
{
    pLattice lattice(new Lattice(1000, 1.e-6, 50.));

    // CaII
    unsigned int Z = 20;
    int Charge = 2;
    OccupationMap occ = ConfigurationParser::ParseFractionalConfiguration("1s2 2s2 2p6 3s2 3p6");

    pCore core(new Core(lattice));
    core->SetOccupancies(occ);

    pIntegrator integrator(new SimpsonsIntegrator(lattice));
    pODESolver ode_solver(new AdamsSolver(integrator));
    pCoulombOperator coulomb(new CoulombOperator(lattice, ode_solver));
    pPhysicalConstant physical_constant(new PhysicalConstant());
    pHFOperator hf(new HFOperator(Z, core, physical_constant, integrator, coulomb));

    HartreeFocker HF_Solver(ode_solver);
    HF_Solver.StartCore(core, hf);
    HF_Solver.SolveCore(core, hf);

    unsigned int pqn = 4;
    double trialE = -0.5 * Charge/(pqn*pqn);
    pOrbital start_4s(new Orbital(-1, pqn, trialE));
    HF_Solver.CalculateExcitedState(start_4s, hf);

    pOrbital new_4s(new Orbital(*start_4s));
    new_4s->SetEnergy(-0.4);

    // Construct operator, check nothing is broken with lambda = 0
    pSpecificMassShiftDecorator t(new SpecificMassShiftDecorator(hf));
    HF_Solver.CalculateExcitedState(new_4s, t);

    EXPECT_NEAR(new_4s->Energy(), -0.41663154, 1.e-6 * 0.41663154);
    EXPECT_NEAR(t->GetMatrixElement(*new_4s, *new_4s), -0.41663154, 1.e-6  * 0.41663154);

    // Check that < 4s | t | 4s > = t.Energy()
    t->SetInverseMass(0.001);
    HF_Solver.CalculateExcitedState(new_4s, t);
    EXPECT_NEAR(new_4s->Norm(integrator), 1.0, 1.e-8);
    EXPECT_NEAR(t->GetMatrixElement(*new_4s, *new_4s), new_4s->Energy(), 1.e-6 * fabs(new_4s->Energy()));
}

TEST(MassShiftDecoratorTester, SrIISlow)
{
    pLattice lattice(new Lattice(1000, 1.e-6, 50.));
    DebugOptions.LogHFIterations(true);

    // SrII
    unsigned int Z = 38;
    int Charge = 2;
    std::string filling = "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6";

    pCore core(new Core(lattice, filling));

    pIntegrator integrator(new SimpsonsIntegrator(lattice));
    pODESolver ode_solver(new AdamsSolver(integrator));
    pCoulombOperator coulomb(new CoulombOperator(lattice, ode_solver));
    pPhysicalConstant physical_constant(new PhysicalConstant());
    HartreeFocker HF_Solver(ode_solver);

    pHFOperator hf(new HFOperator(Z, core, physical_constant, integrator, coulomb));
    pSpecificMassShiftDecorator t(new SpecificMassShiftDecorator(hf, true, false));

    HF_Solver.StartCore(core, hf);
    HF_Solver.SolveCore(core, t);

    unsigned int pqn = 5;
    double trialE = -0.5 * Charge/(pqn*pqn);
    pOrbital new_5s(new Orbital(-1, pqn, trialE));

    HF_Solver.CalculateExcitedState(new_5s, t);
    new_5s->CheckSize(lattice, 1.e-10);

    // Test that gradient of SMS operator is correct (compared to old values)
    t->SetInverseMass(0.001);
    HF_Solver.SolveCore(core, t);
    HF_Solver.CalculateExcitedState(new_5s, t);
    double Eplus = new_5s->Energy();

    t->SetInverseMass(-0.001);
    HF_Solver.SolveCore(core, t);
    HF_Solver.CalculateExcitedState(new_5s, t);
    double Eminus = new_5s->Energy();

    double sms = -660./3609.;
    EXPECT_NEAR((Eplus-Eminus)/0.002, sms, 1.e-3 * fabs(sms));
}

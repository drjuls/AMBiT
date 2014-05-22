#include "HFOperator.h"
#include "gtest/gtest.h"
#include "Core.h"
#include "Include.h"
#include "ODESolver.h"
#include "Universal/MathConstant.h"
#include "HartreeFocker.h"
#include "NucleusDecorator.h"

TEST(HartreeFockerTester, CaIIOrbital)
{
    pLattice lattice(new Lattice(1500, 1.e-6, 100.));

    // Ca
    unsigned int Z = 20;
    int Charge = 2;
    std::string filling = "1s2 2s2 2p6 3s2 3p6";

    DebugOptions.LogFirstBuild(true);
    DebugOptions.LogHFIterations(true);
    DebugOptions.HartreeEnergyUnits(true);

    pCore core(new Core(lattice, filling));

    // Set up HF ODE and HartreeFocker
    pOPIntegrator integrator(new SimpsonsIntegrator(lattice));
    pODESolver ode_solver(new AdamsSolver(integrator));
    pCoulombOperator coulomb(new CoulombOperator(lattice, ode_solver));
    pPhysicalConstant physical_constant(new PhysicalConstant());
    pHFOperator t(new HFOperator(Z, core, physical_constant, integrator, coulomb));
    HartreeFocker HF_Solver(ode_solver);

    HF_Solver.StartCore(core, t);
    HF_Solver.SolveCore(core, t);

    EXPECT_NEAR(core->GetState(OrbitalInfo(1, -1))->Energy(), -150.717923115, 1.e-6 * 150.717923115);
    EXPECT_NEAR(core->GetState(OrbitalInfo(2, -1))->Energy(), -17.5158164976, 1.e-6 * 17.5158164976);
    EXPECT_NEAR(core->GetState(OrbitalInfo(3, -1))->Energy(), -2.79675319452, 1.e-6 * 2.79675319452);
    EXPECT_NEAR(core->GetState(OrbitalInfo(2,  1))->Energy(), -14.2827896685, 1.e-6 * 14.2827896685);
    EXPECT_NEAR(core->GetState(OrbitalInfo(3, -2))->Energy(), -1.87184451641, 1.e-6 * 1.87184451641);

    unsigned int pqn = 4;
    double trialE = -0.5 * Charge/(pqn*pqn);
    pOrbital start_4s(new Orbital(-1, pqn, trialE));

    HF_Solver.CalculateExcitedState(start_4s, t);

    EXPECT_NEAR(start_4s->Energy(), -0.41663154, 0.000001 * 0.41663154);
    EXPECT_NEAR(start_4s->Norm(integrator), 1.0, 1.e-8);

    // Another 4s orbital
    pOrbital new_4s(new Orbital(-1, 4, -0.4));
    HF_Solver.CalculateExcitedState(new_4s, t);
    EXPECT_NEAR(new_4s->Energy(), -0.41663154, 1.e-6 * 0.41663154);
    EXPECT_NEAR(new_4s->Norm(integrator), 1.0, 1.e-8);
    EXPECT_NEAR(t->GetMatrixElement(*new_4s, *new_4s), -0.41663154, 1.e-6 * 0.41663154);

    // Create 5d orbital
//    pOrbital new_5d(new Orbital(2, 5, -0.1));
//    HF_Solver.CalculateExcitedState(new_5d, t);
//    EXPECT_NEAR(new_5d->Energy(), -0.10135136, 1.e-6 * 0.10135136);

    // Check core orbital
    pOrbital new_2p(new Orbital(*core->GetState(OrbitalInfo(2, 1))));
    new_2p->SetEnergy(-12.0);
    HF_Solver.SolveOrbital(new_2p, t);
    EXPECT_NEAR(new_2p->Energy(), -14.282789, 1.e-6 * 14.282789);
}

TEST(HartreeFockerTester, ContinuumOrbital)
{
    pLattice lattice(new Lattice(1500, 1.e-6, 100.));

    // Ca
    unsigned int Z = 12;
    std::string filling = "1s2 2s2 2p6";

    DebugOptions.LogFirstBuild(true);
    DebugOptions.LogHFIterations(true);
    DebugOptions.HartreeEnergyUnits(true);

    pCore core(new Core(lattice, filling));

    // Set up HF ODE and HartreeFocker
    pOPIntegrator integrator(new SimpsonsIntegrator(lattice));
    pODESolver ode_solver(new AdamsSolver(integrator));
    pCoulombOperator coulomb(new CoulombOperator(lattice, ode_solver));
    pPhysicalConstant physical_constant(new PhysicalConstant());
    pHFOperator t(new HFOperator(Z, core, physical_constant, integrator, coulomb));
    HartreeFocker HF_Solver(ode_solver);

    HF_Solver.StartCore(core, t);
    HF_Solver.SolveCore(core, t);

    pContinuumWave epsilon(new ContinuumWave(-1, 100, 0.2));    // kappa, pqn, energy

    unsigned int loop = HF_Solver.CalculateContinuumWave(epsilon, t);
    EXPECT_NE(0, loop);
    epsilon->Print();
}

TEST(HartreeFockerTester, AuContinuum)
{
    pLattice lattice(new Lattice(1500, 1.e-6, 50.));
    lattice->resize(600.);

    // Au25+
    unsigned int Z = 79;
    std::string filling = "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d10 4f8";

    DebugOptions.LogFirstBuild(true);
    DebugOptions.LogHFIterations(true);
    DebugOptions.HartreeEnergyUnits(true);

    pCore core(new Core(lattice, filling));

    // Set up HF ODE and HartreeFocker
    pOPIntegrator integrator(new SimpsonsIntegrator(lattice));
    pODESolver ode_solver(new AdamsSolver(integrator));
    pCoulombOperator coulomb(new CoulombOperator(lattice, ode_solver));
    pPhysicalConstant physical_constant(new PhysicalConstant());
    pHFOperator t(new HFOperator(Z, core, physical_constant, integrator, coulomb));
    HartreeFocker HF_Solver(ode_solver);

    pNucleusDecorator nuc(new NucleusDecorator(t));
    nuc->SetFermiParameters(6.5544);
    t = nuc;

    HF_Solver.StartCore(core, t);
    HF_Solver.SolveCore(core, t);

    double eV = 1./MathConstant::Instance()->HartreeEnergyIneV();
    pContinuumWave epsilon(new ContinuumWave(-1, 100, 1. * eV));    // kappa, pqn, energy

    unsigned int loop = HF_Solver.CalculateContinuumWave(epsilon, t);
    EXPECT_NE(0, loop);
    epsilon->Print(lattice);
}

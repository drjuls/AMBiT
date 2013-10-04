#include "HFOperator.h"
#include "gtest/gtest.h"
#include "Core.h"
#include "Include.h"
#include "ODESolver.h"
#include "Universal/PhysicalConstant.h"
#include "HartreeFocker.h"

TEST(HFOperatorTester, ODESolver)
{
    Lattice lattice(1000, 1.e-6, 50.);

    // Ca
    unsigned int Z = 20;
    int Charge = 2;
    std::string filling = "1s2 2s2 2p6 3s2 3p6";

    DebugOptions.LogFirstBuild(true);
    DebugOptions.LogHFIterations(true);

    Core core(&lattice, Z, Charge);
    core.Initialise(filling);
    core.Update();

    unsigned int pqn = 4;
    double trialE = -0.5 * Charge/(pqn*pqn);
    Orbital start_4s(-1, trialE, pqn);
    core.CalculateExcitedState(&start_4s);
    start_4s.CheckSize(&lattice, 1.e-10);

    EXPECT_NEAR(start_4s.GetEnergy(), -0.41663154, 0.000001 * 0.41663154);
    EXPECT_NEAR(start_4s.Norm(&lattice), 1.0, 1.e-8);

    Orbital new_4s(start_4s);
    new_4s.SetEnergy(-0.4);

    SimpsonsIntegrator integrator(&lattice);
    AdamsSolver ode_solver(&lattice);
    CoulombOperator coulomb(&lattice, &ode_solver);
    HFOperator t(Z, &core, &integrator, &coulomb);

    HartreeFocker HF_Solver;
    HF_Solver.SolveOrbital(&new_4s, &t);

    EXPECT_NEAR(new_4s.GetEnergy(), -0.41663154, 1.e-6 * 0.41663154);
    EXPECT_NEAR(new_4s.Norm(&lattice), 1.0, 1.e-8);
    EXPECT_NEAR(t.GetMatrixElement(new_4s, new_4s), -0.41663154, 1.e-6  * 0.41663154);

    // Check core orbital
    Orbital new_2p(*core.GetState(OrbitalInfo(2, 1)));
    new_2p.SetEnergy(-12.0);
    HF_Solver.SolveOrbital(&new_2p, &t);

    EXPECT_NEAR(new_2p.GetEnergy(), -14.282789, 1.e-6 * 14.282789);    
}
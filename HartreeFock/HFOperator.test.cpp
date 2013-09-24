#include "HFOperator.h"
#include "gtest/gtest.h"
#include "Core.h"
#include "Include.h"
#include "ODESolver.h"
#include "GreensMethodOperator.h"
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

    EXPECT_NEAR(start_4s.GetEnergy(), -0.41663154257325596, 0.000001 * 0.41663154257325596);
    EXPECT_NEAR(start_4s.Norm(&lattice), 1.0, 1.e-4);

    Orbital new_4s(start_4s);
    new_4s.SetEnergy(-0.4);

    OPIntegrator* integrator = NULL;
    HFOperator t(Z, &core, integrator);

    HartreeFocker HF_Solver;
    HF_Solver.IterateOrbital(&new_4s, &t);

    EXPECT_NEAR(new_4s.GetEnergy(), -0.41663154257325596, 0.000001 * 0.41663154257325596);
    EXPECT_NEAR(new_4s.Norm(&lattice), 1.0, 1.e-8);
}

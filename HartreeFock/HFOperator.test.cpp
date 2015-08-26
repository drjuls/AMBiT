#include "HFOperator.h"
#include "gtest/gtest.h"
#include "Core.h"
#include "Include.h"
#include "ODESolver.h"
#include "Universal/PhysicalConstant.h"
#include "HartreeFocker.h"
#include "ConfigurationParser.h"

TEST(HFOperatorTester, ODESolver)
{
    pLattice lattice(new Lattice(1000, 1.e-6, 50.));

    // Ca
    unsigned int Z = 20;
    int Charge = 2;
    OccupationMap filling = ConfigurationParser::ParseFractionalConfiguration("1s2 2s2 2p6 3s2 3p6");

    DebugOptions.LogFirstBuild(true);
    DebugOptions.LogHFIterations(true);

    pCore core(new Core(lattice));
    core->SetOccupancies(filling);

    pOPIntegrator integrator(new SimpsonsIntegrator(lattice));
    pODESolver ode_solver(new AdamsSolver(integrator));
    pCoulombOperator coulomb(new CoulombOperator(lattice, ode_solver));
    pPhysicalConstant physical_constant(new PhysicalConstant());
    pHFOperator t(new HFOperator(Z, core, physical_constant, integrator, coulomb));

    HartreeFocker HF_Solver(ode_solver);
    HF_Solver.StartCore(core, t);
    HF_Solver.SolveCore(core, t);

    unsigned int pqn = 4;
    double trialE = -0.5 * Charge/(pqn*pqn);
    pOrbital new_4s(new Orbital(-1, pqn, trialE));
    new_4s->SetEnergy(-0.4);

    HF_Solver.CalculateExcitedState(new_4s, t);

    EXPECT_NEAR(new_4s->Energy(), -0.41663154, 1.e-6 * 0.41663154);
    EXPECT_NEAR(new_4s->Norm(integrator), 1.0, 1.e-8);
    EXPECT_NEAR(t->GetMatrixElement(*new_4s, *new_4s), -0.41663154, 1.e-6  * 0.41663154);

    // Check core orbital
    pOrbital new_2p(new Orbital(*core->GetState(OrbitalInfo(2, 1))));
    new_2p->SetEnergy(-12.0);
    HF_Solver.ConvergeOrbitalAndExchange(new_2p, t, &HartreeFocker::IterateOrbital, 1.e-12);

    EXPECT_NEAR(new_2p->Energy(), -14.282789, 1.e-6 * 14.282789);
}

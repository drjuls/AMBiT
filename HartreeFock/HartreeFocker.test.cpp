#include "HFOperator.h"
#include "gtest/gtest.h"
#include "Core.h"
#include "Include.h"
#include "ODESolver.h"
#include "Universal/PhysicalConstant.h"
#include "HartreeFocker.h"

TEST(HartreeFockerTester, CaIIOrbital)
{
    pLattice lattice(new Lattice(1500, 1.e-6, 100.));

    // Ca
    unsigned int Z = 20;
    int Charge = 2;
    std::string filling = "1s2 2s2 2p6 3s2 3p6";

    // Bootstrap CaII core using old code
    DebugOptions.LogFirstBuild(true);
    DebugOptions.LogHFIterations(true);

    Core core(lattice, Z, Charge);
    core.Initialise(filling);
    core.Update();

    unsigned int pqn = 4;
    double trialE = -0.5 * Charge/(pqn*pqn);
    pOrbital start_4s(new Orbital(-1, trialE, pqn));
    core.CalculateExcitedState(start_4s);
    start_4s->CheckSize(lattice, 1.e-10);

    EXPECT_NEAR(start_4s->GetEnergy(), -0.41663154, 0.000001 * 0.41663154);
    EXPECT_NEAR(start_4s->Norm(lattice), 1.0, 1.e-8);

    // Set up HF ODE and HartreeFocker
    pOPIntegrator integrator(new SimpsonsIntegrator(lattice));
    AdamsSolver ode_solver(lattice);
    pCoulombOperator coulomb(new CoulombOperator(lattice, &ode_solver));
    pHFOperator t(new HFOperator(Z, &core, integrator, coulomb));
    HartreeFocker HF_Solver(&ode_solver);

    // Create 4s orbital
    pOrbital new_4s(new Orbital(-1, -0.4, 4));
    HF_Solver.CalculateExcitedState(new_4s, t);
    EXPECT_NEAR(new_4s->GetEnergy(), -0.41663154, 1.e-6 * 0.41663154);
    EXPECT_NEAR(new_4s->Norm(lattice), 1.0, 1.e-8);
    EXPECT_NEAR(t->GetMatrixElement(*new_4s, *new_4s), -0.41663154, 1.e-6 * 0.41663154);

    // Create 5d orbital
    /* TODO: HartreeFocker::CalculateExcitedState() can't cope with this condition.
    Orbital new_5d(2, -0.1, 5);
    HF_Solver.CalculateExcitedState(&new_5d, &t);
    EXPECT_NEAR(new_5d.GetEnergy(), -0.10135136, 1.e-6 * 0.10135136);
    */

    // Check core orbital
    pOrbital new_2p(new Orbital(*core.GetState(OrbitalInfo(2, 1))));
    new_2p->SetEnergy(-12.0);
    HF_Solver.SolveOrbital(new_2p, t);
    EXPECT_NEAR(new_2p->GetEnergy(), -14.282789, 1.e-6 * 14.282789);
}

TEST(HartreeFockerTester, CaIICore)
{
    pLattice lattice(new Lattice(1000, 1.e-6, 50.));
    
    // Ca
    unsigned int Z = 20;
    int Charge = 2;
    std::string filling = "1s2 2s2 2p6 3s2 3p6";
    
    DebugOptions.LogFirstBuild(true);
    DebugOptions.LogHFIterations(true);
    
    Core core(lattice, Z, Charge);
    core.Initialise(filling);
    core.Update();
    
    pOPIntegrator integrator(new SimpsonsIntegrator(lattice));
    AdamsSolver ode_solver(lattice);
    pCoulombOperator coulomb(new CoulombOperator(lattice, &ode_solver));
    pHFOperator t(new HFOperator(Z, &core, integrator, coulomb));

    HartreeFocker HF_Solver(&ode_solver);
    HF_Solver.SolveCore(&core, t);

    EXPECT_NEAR(core.GetState(OrbitalInfo(1, -1))->GetEnergy(), -150.717923115, 1.e-6 * 150.717923115);
    EXPECT_NEAR(core.GetState(OrbitalInfo(2, -1))->GetEnergy(), -17.5158164976, 1.e-6 * 17.5158164976);
    EXPECT_NEAR(core.GetState(OrbitalInfo(3, -1))->GetEnergy(), -2.79675319452, 1.e-6 * 2.79675319452);
    EXPECT_NEAR(core.GetState(OrbitalInfo(2,  1))->GetEnergy(), -14.2827896685, 1.e-6 * 14.2827896685);
    EXPECT_NEAR(core.GetState(OrbitalInfo(3, -2))->GetEnergy(), -1.87184451641, 1.e-6 * 1.87184451641);
}

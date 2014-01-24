#include "MassShiftDecorator.h"
#include "gtest/gtest.h"
#include "Core.h"
#include "Include.h"
#include "ODESolver.h"
#include "Universal/PhysicalConstant.h"
#include "HartreeFocker.h"

TEST(MassShiftDecoratorTester, CaII)
{
    pLattice lattice(new Lattice(1000, 1.e-6, 50.));

    // CaII
    unsigned int Z = 20;
    int Charge = 2;
    std::string filling = "1s2 2s2 2p6 3s2 3p6";

    Core core(lattice, Z, Charge);
    core.Initialise(filling);
    core.Update();

    unsigned int pqn = 4;
    double trialE = -0.5 * Charge/(pqn*pqn);
    pOrbital start_4s(new Orbital(-1, trialE, pqn));
    core.CalculateExcitedState(start_4s);
    start_4s->CheckSize(lattice, 1.e-10);

    pOrbital new_4s(new Orbital(*start_4s));
    new_4s->SetEnergy(-0.4);

    pOPIntegrator integrator(new SimpsonsIntegrator(lattice));
    AdamsSolver ode_solver(lattice);
    pCoulombOperator coulomb(new CoulombOperator(lattice, &ode_solver));

    // Construct operator, check nothing is broken with lambda = 0
    pHFOperator hf(new HFOperator(Z, &core, integrator, coulomb));
    pMassShiftDecorator t(new MassShiftDecorator(hf));

    HartreeFocker HF_Solver(&ode_solver);
    HF_Solver.SolveOrbital(new_4s, t);

    EXPECT_NEAR(new_4s->GetEnergy(), -0.41663154, 1.e-6 * 0.41663154);
    EXPECT_NEAR(t->GetMatrixElement(*new_4s, *new_4s), -0.41663154, 1.e-6  * 0.41663154);

    // Check that < 4s | t | 4s > = t.GetEnergy()
    t->SetInverseMass(0.001);
    HF_Solver.SolveOrbital(new_4s, t);
    EXPECT_NEAR(new_4s->Norm(lattice), 1.0, 1.e-8);
    EXPECT_NEAR(t->GetMatrixElement(*new_4s, *new_4s), new_4s->GetEnergy(), 1.e-6 * fabs(new_4s->GetEnergy()));
}

TEST(MassShiftDecoratorTester, SrII)
{
    pLattice lattice(new Lattice(1000, 1.e-6, 50.));
    DebugOptions.LogHFIterations(true);
    
    // SrII
    unsigned int Z = 38;
    int Charge = 2;
    std::string filling = "1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6";

    Core core(lattice, Z, Charge);
    core.Initialise(filling);
    core.Update();

    unsigned int pqn = 5;
    double trialE = -0.5 * Charge/(pqn*pqn);
    pOrbital start_5s(new Orbital(-1, trialE, pqn));
    core.CalculateExcitedState(start_5s);
    start_5s->CheckSize(lattice, 1.e-10);
    
    pOrbital new_5s(new Orbital(*start_5s));
    
    pOPIntegrator integrator(new SimpsonsIntegrator(lattice));
    AdamsSolver ode_solver(lattice);
    pCoulombOperator coulomb(new CoulombOperator(lattice, &ode_solver));

    pHFOperator hf(new HFOperator(Z, &core, integrator, coulomb));
    pMassShiftDecorator t(new MassShiftDecorator(hf));
    
    // Test that gradient of SMS operator is correct (compared to old values)
    t->SetInverseMass(0.001);
    HartreeFocker HF_Solver(&ode_solver);
    HF_Solver.SolveCore(&core, t);
    HF_Solver.SolveOrbital(new_5s, t);
    double Eplus = new_5s->GetEnergy();

    t->SetInverseMass(-0.001);
    HF_Solver.SolveCore(&core, t);
    HF_Solver.SolveOrbital(new_5s, t);
    double Eminus = new_5s->GetEnergy();

    double sms = -660./3609.;
    EXPECT_NEAR((Eplus-Eminus)/0.002, sms, 1.e-3 * fabs(sms));
}

#include "HFOperator.h"
#include "gtest/gtest.h"
#include "Core.h"
#include "Include.h"
#include "ODESolver.h"
#include "GreensMethodOperator.h"
#include "Universal/PhysicalConstant.h"

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

    Orbital origin_4s(start_4s);
    Orbital inf_4s(start_4s);
    SpinorFunction exchange(start_4s.Kappa());
    core.CalculateExchange(start_4s, exchange);

    OPIntegrator* integrator = NULL;
    HFOperator t(Z, &core, integrator);

    ODESolver* odesolver = new AdamsSolver(&lattice);
    t.SetODEParameters(start_4s.Kappa(), start_4s.GetEnergy());
    odesolver->IntegrateBackwards(&t, &inf_4s);
    odesolver->IntegrateForwards(&t, &origin_4s);

    GreensMethodOperator greens(&lattice);
    greens.SetHomogenousSolutions(&origin_4s, &inf_4s);

    Orbital new_4s(start_4s);
    std::vector<double> greensfunction0(start_4s.Size());
    std::vector<double> greensfunctionInf(start_4s.Size());
    exchange *= PhysicalConstant::Instance()->GetAlpha();
    greens.SetSourceTerm(&exchange, true);
    odesolver->IntegrateForwards(&greens, &greensfunction0);
    greens.SetSourceTerm(&exchange, false);
    odesolver->IntegrateBackwards(&greens, &greensfunctionInf);

    for(unsigned int i = 0; i < start_4s.Size(); i++)
    {
        new_4s.f[i] = origin_4s.f[i] * greensfunctionInf[i] - inf_4s.f[i] * greensfunction0[i];
        new_4s.g[i] = origin_4s.g[i] * greensfunctionInf[i] - inf_4s.g[i] * greensfunction0[i];
    }

    EXPECT_NEAR(new_4s.Norm(&lattice), 1.0, 1.e-6);

    delete odesolver;
}

#include "HartreeFocker.h"
#include "Include.h"
#include "Universal/PhysicalConstant.h"
#include "ODESolver.h"
#include "GreensMethodOperator.h"

void HartreeFocker::IterateCore(Core* core, HFOperator* hf, SpinorODE* hf_decorated)
{
}

void HartreeFocker::IterateOrbital(Orbital* orbital, SpinorODE* hf)
{
    Lattice* lattice = hf->GetLattice();
    ODESolver* odesolver = new AdamsSolver(lattice);

    double delta_E = 0.0;

    do
    {   // Use greens method to iterate orbital
        orbital->ReNormalise(lattice);
        orbital->CheckSize(lattice, WavefunctionTolerance);

        hf->SetODEParameters(orbital);
        SpinorFunction exchange(hf->GetExchange(orbital));

        // Get solutions to homogenous equation (no exchange)
        Orbital originregular(*orbital);
        Orbital infinityregular(*orbital);

        hf->IncludeExchangeInODE(false);
        odesolver->IntegrateBackwards(hf, &infinityregular);
        odesolver->IntegrateForwards(hf, &originregular);

        GreensMethodOperator greens(lattice);
        greens.SetHomogenousSolutions(originregular, infinityregular);
        
        std::vector<double> G0(orbital->Size());
        std::vector<double> dG0dr(orbital->Size());
        std::vector<double> GInf(orbital->Size());
        std::vector<double> dGInfdr(orbital->Size());

        exchange *= PhysicalConstant::Instance()->GetAlpha();
        greens.SetSourceTerm(exchange, true);
        odesolver->IntegrateForwards(&greens, &G0, &dG0dr);
        greens.SetSourceTerm(exchange, false);
        odesolver->IntegrateBackwards(&greens, &GInf, &dGInfdr);

        *orbital = originregular.TimesVector(GInf, dGInfdr) - infinityregular.TimesVector(G0, dG0dr);

        // Now modify energy if required
        double norm = orbital->Norm(lattice);

        // Get delta_psi using Greens operator with psi as the source term
        greens.SetSourceTerm(*orbital, true);
        odesolver->IntegrateForwards(&greens, &G0, &dG0dr);
        greens.SetSourceTerm(*orbital, false);
        odesolver->IntegrateBackwards(&greens, &GInf, &dGInfdr);

        Orbital delta_psi = originregular.TimesVector(GInf, dGInfdr) - infinityregular.TimesVector(G0, dG0dr);
        delta_psi *= PhysicalConstant::Instance()->GetAlpha();

        double var = orbital->Overlap(delta_psi, lattice);

        delta_E = (1. - norm)/(2. * var);
        double energy = orbital->GetEnergy();
        orbital->SetEnergy(energy + delta_E);

    } while (fabs(delta_E) > WavefunctionEnergyTolerance);

    delete odesolver;
}

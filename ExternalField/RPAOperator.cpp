#include "RPAOperator.h"
#include "Include.h"

RPAOperator::RPAOperator(pSpinorOperator external, pHFOperatorConst hf, pHartreeY hartreeY, pRPASolver rpa_solver):
    TimeDependentSpinorOperator(external->GetK(), external->GetParity(), external->GetIntegrator()),
    external(external), hf(hf), hartreeY(hartreeY), core(nullptr), solver(rpa_solver), scale(1.0)
{
    // Check if static RPA
    pTimeDependentSpinorOperator tdop = std::dynamic_pointer_cast<TimeDependentSpinorOperator>(external);
    static_rpa = (tdop == nullptr);

    // Get rpa core: replace orbitals with RPAOrbitals
    MathConstant* math = MathConstant::Instance();
    int twoK = 2 * external->GetK();

    pCoreConst hf_core = hf->GetCore();
    core = std::make_shared<Core>(hf_core->GetLattice());

    for(const auto& orb: *hf_core)
    {
        pRPAOrbital rpa_orb = std::make_shared<RPAOrbital>(*orb.second);

        Parity Pdelta = rpa_orb->GetParity() * external->GetParity();
        rpa_orb->deltapsi.clear();

        for(int twoj = mmax(1, rpa_orb->TwoJ() - twoK); twoj <= rpa_orb->TwoJ() + twoK; twoj+=2)
        {
            int kappa = math->convert_to_kappa(twoj, Pdelta);
            if(rpa_orb->GetDeltaPsi(kappa) == rpa_orb->deltapsi.end())
            {
                if(static_rpa)
                    rpa_orb->deltapsi.push_back(std::make_pair(std::make_shared<DeltaOrbital>(kappa, rpa_orb), nullptr));
                else
                    rpa_orb->deltapsi.push_back(std::make_pair(std::make_shared<DeltaOrbital>(kappa, rpa_orb), std::make_shared<DeltaOrbital>(kappa, rpa_orb)));
            }
        }

        core->AddState(rpa_orb);
    }

    core->SetOccupancies(hf_core->GetOccupancies());
}

void RPAOperator::SolveRPA()
{
    // Solve RPA core.
    solver->SolveRPACore(hf, std::static_pointer_cast<RPAOperator>(shared_from_this()));
}

void RPAOperator::SetFrequency(double frequency)
{
    if(!static_rpa)
    {
        omega = fabs(frequency);

        pTimeDependentSpinorOperator tdop = std::dynamic_pointer_cast<TimeDependentSpinorOperator>(external);
        if(tdop)
            tdop->SetFrequency(omega);
    }
}

void RPAOperator::ClearRPACore()
{
    for(auto& orb: *core)
    {
        pRPAOrbital rpa_orb = std::make_shared<RPAOrbital>(*orb.second);
        for(auto& pair: rpa_orb->deltapsi)
        {
            if(pair.first)
            {   pair.first->Clear();
                pair.first->SetDeltaEnergy(0.0);
            }
            if(pair.second)
            {   pair.second->Clear();
                pair.second->SetDeltaEnergy(0.0);
            }
        }
    }
}

SpinorFunction RPAOperator::ReducedApplyTo(const SpinorFunction& a, int kappa_b) const
{
    return ReducedApplyTo(a, kappa_b, false);
}

SpinorFunction RPAOperator::ConjugateReducedApplyTo(const SpinorFunction& a, int kappa_b) const
{
    return ReducedApplyTo(a, kappa_b, !static_rpa);
}

SpinorFunction RPAOperator::ReducedApplyTo(const SpinorFunction& a, int kappa_b, bool conjugate) const
{
    SpinorFunction ret(kappa_b);
    MathConstant* math = MathConstant::Instance();

    // To save confusion: kappa_b -> alph
    int kappa_alph = kappa_b;
    int alph_L = (kappa_b > 0? kappa_b: -kappa_b-1);
    int alph_TwoJ = 2 * abs(kappa_alph) - 1;
    double alph_J = double(alph_TwoJ)/2;

    pSpinorFunctionConst pa = a.shared_from_this();

    // First term: <kappa || f || parent>
    if(conjugate)
        ret = std::static_pointer_cast<TimeDependentSpinorOperator>(external)->ConjugateReducedApplyTo(a, kappa_alph) * scale;
    else
        ret = external->ReducedApplyTo(a, kappa_alph) * scale;

    // Second term: <kappa || deltaV || parent>
    // Sum deltaV for all core states. Core states are `b', delta is `beta'.
    for(const auto& cs: *core)
    {
        pRPAOrbitalConst b(std::dynamic_pointer_cast<const RPAOrbital>(cs.second));
        if(b == nullptr)
            continue;

        // Sum over beta
        for(const auto& cs_delta: b->deltapsi)
        {
            pDeltaOrbitalConst beta, betaplus;

            if(static_rpa || conjugate == false)
            {
                beta = cs_delta.first;
                betaplus = cs_delta.second;
            }
            else
            {
                beta = cs_delta.second;
                betaplus = cs_delta.first;
            }

            double occupancy_factor = (a.TwoJ()+1) * (b->TwoJ()+1) * (alph_TwoJ+1) * (beta->TwoJ()+1);
            occupancy_factor = sqrt(occupancy_factor);
            // TODO: Open shells need to be scaled

            // Direct term: a -> alph, b -> beta
            if((alph_L + a.L() + K)%2 == 0 && (b->L() + beta->L() + K)%2 == 0)
            {
                double coeff = 1./(2. * K + 1.) * math->minus_one_to_the_power((a.TwoJ() + b->TwoJ())/2 + 1);
                coeff *= math->Electron3j(alph_TwoJ, a.TwoJ(), K) * math->Electron3j(b->TwoJ(), beta->TwoJ(), K);

                if(coeff)
                {
                    coeff *= occupancy_factor;

                    if(static_rpa)
                    {
                        hartreeY->SetParameters(K, b, beta);
                        ret += hartreeY->ApplyTo(a, ret.Kappa()) * coeff * 2.;
                    }
                    else
                    {
                        hartreeY->SetParameters(K, b, beta);
                        ret += hartreeY->ApplyTo(a, ret.Kappa()) * coeff;

                        hartreeY->SetParameters(K, betaplus, b);
                        ret += hartreeY->ApplyTo(a, ret.Kappa()) * coeff;
                    }
                }
            }

            // C: alph -> b; beta -> a
            // Sum over all k
            if((alph_L + b->L())%2 == (a.L() + beta->L())%2)
            {
                int mink = mmax(abs(a.TwoJ() - beta->TwoJ()), abs(alph_TwoJ - b->TwoJ()))/2;
                if((alph_L + b->L() + mink)%2)
                    mink++;
                int maxk = mmin(a.TwoJ() + beta->TwoJ(), alph_TwoJ + b->TwoJ());
                for(int k = mink; k <= maxk; k+=2)
                {
                    double coeff = math->Wigner6j(K, alph_J, a.J(), k, beta->J(), b->J()) *
                        math->minus_one_to_the_power((a.TwoJ() + b->TwoJ())/2 + k + K + 1);
                    coeff *= math->Electron3j(a.TwoJ(), beta->TwoJ(), k) * math->Electron3j(alph_TwoJ, b->TwoJ(), k);

                    if(coeff)
                    {
                        coeff *= occupancy_factor;

                        if(static_rpa)
                        {
                            hartreeY->SetParameters(k, beta, pa);
                            ret += hartreeY->ApplyTo(*b, ret.Kappa()) * coeff;
                        }
                        else
                        {
                            hartreeY->SetParameters(k, betaplus, pa);
                            ret += hartreeY->ApplyTo(*b, ret.Kappa()) * coeff;
                        }
                    }
                }
            }

            // D: alph->beta; b -> a
            // Sum over all k
            if((a.L() + b->L())%2 == (alph_L + beta->L())%2)
            {
                int mink = mmax(abs(a.TwoJ() - b->TwoJ()), abs(alph_TwoJ - beta->TwoJ()))/2;
                if((alph_L + beta->L() + mink)%2)
                    mink++;
                int maxk = mmin(a.L() + b->L(), alph_L + beta->L());
                for(int k = mink; k <= maxk; k+=2)
                {
                    double coeff = math->Wigner6j(K, a.J(), alph_J, k, beta->J(), b->J()) *
                        math->minus_one_to_the_power((a.TwoJ() + b->TwoJ())/2 + k + K + 1);
                    coeff *= math->Electron3j(a.TwoJ(), b->TwoJ(), k) * math->Electron3j(alph_TwoJ, beta->TwoJ(), k);

                    if(coeff)
                    {
                        coeff *= occupancy_factor;

                        hartreeY->SetParameters(k, b, pa);
                        ret += hartreeY->ApplyTo(*beta, ret.Kappa()) * coeff;
                    }
                }
            }
        }
    }
    
    return ret;
}

RadialFunction RPAOperator::GetRPAField() const
{
    RadialFunction ret;
    RadialFunction density, pot;
    MathConstant* math = MathConstant::Instance();

    pIntegrator integrator(new SimpsonsIntegrator(hf->GetLattice()));
    pODESolver ode_solver(new AdamsSolver(integrator));
    pCoulombOperator coulomb(new CoulombOperator(hf->GetLattice(), ode_solver));

    auto R = hf->GetLattice()->R();
    RadialFunction Rp1, Rp2, Rm1;
    Rp1.f = std::vector<double>(R, R+hf->GetLattice()->size());
    Rp1.dfdr = std::vector<double>(hf->GetLattice()->size(), 1.);

    Rp2 = Rp1 * Rp1;

    Rm1.resize(Rp1.size());
    for(int i = 0; i < Rp1.size(); i++)
    {
        Rm1.f[i] = 1./R[i];
        Rm1.dfdr[i] = -Rm1.f[i]/R[i];
    }

    // Sum deltaV for all core states. Core states are `b', delta is `beta'.
    for(const auto& cs: *core)
    {
        pRPAOrbitalConst b(std::dynamic_pointer_cast<const RPAOrbital>(cs.second));
        if(b == nullptr)
            continue;

        // Sum over beta
        for(const auto& cs_delta: b->deltapsi)
        {
            pDeltaOrbitalConst beta, betaplus;
            beta = cs_delta.first;
            betaplus = cs_delta.second;

            double occupancy_factor = (b->TwoJ()+1) * (beta->TwoJ()+1);
//            occupancy_factor = sqrt(occupancy_factor);
            // TODO: Open shells need to be scaled

            // Direct term
            if((b->L() + beta->L() + K)%2 == 0)
            {
                double coeff = 1./(2. * K + 1.) /** math->minus_one_to_the_power((a.TwoJ() + b->TwoJ())/2 + 1)*/;
//                coeff *= math->Electron3j(b->TwoJ(), beta->TwoJ(), K);
                coeff *= math->SphericalTensorReducedMatrixElement(b->Kappa(), beta->Kappa(), K);

                if(coeff)
                {
//                    coeff *= occupancy_factor;

                    if(static_rpa)
                    {
                        density = b->GetDensity(*beta);
                        density.resize(integrator->GetLattice()->size());
//                        coulomb->GetForwardPotential(K, density, pot);
//                        ret += pot * Rp2 * coeff * 2.;

//                        coulomb->GetBackwardPotential(K, density, pot);
//                        ret += pot * Rm1 * coeff * 2.;

                        coulomb->GetPotential(K, density, pot);
                        ret += pot * coeff * 2.;
                    }
                    else
                    {
                        density = b->GetDensity(*beta);
                        density.resize(integrator->GetLattice()->size());
//                        coulomb->GetForwardPotential(K, density, pot);
//                        ret += pot * Rp2 * coeff;
//                        coulomb->GetBackwardPotential(K, density, pot);
//                        ret += pot * Rm1 * coeff;
                        coulomb->GetPotential(K, density, pot);
                        ret += pot * coeff;

                        density = b->GetDensity(*betaplus);
                        density.resize(integrator->GetLattice()->size());
//                        coulomb->GetForwardPotential(K, density, pot);
//                        ret += pot * Rp2 * coeff;
//                        coulomb->GetBackwardPotential(K, density, pot);
//                        ret += pot * Rm1 * coeff;
                        coulomb->GetPotential(K, density, pot);
                        ret += pot * coeff;
                    }
                }
            }
        }
    }

    return ret * (1./scale);
}

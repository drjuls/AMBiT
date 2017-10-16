#include "CoreMBPTCalculator.h"
#include "Universal/PhysicalConstant.h"
#include "Universal/MathConstant.h"

#ifdef AMBIT_USE_OPENMP
#include <omp.h>
#endif

CoreMBPTCalculator::CoreMBPTCalculator(pOrbitalManagerConst orbitals, pHFIntegrals one_body, pSlaterIntegrals two_body, const std::string& fermi_orbitals):
    MBPTCalculator(orbitals, fermi_orbitals, two_body->OffParityExists()), one_body(one_body), two_body(two_body), core(orbitals->core), excited(orbitals->excited)
{}

CoreMBPTCalculator::~CoreMBPTCalculator()
{
    one_body->clear();
    two_body->clear();
}

unsigned int CoreMBPTCalculator::GetStorageSize()
{
    pOrbitalMap excited_valence = std::make_shared<OrbitalMap>(*excited);
    excited_valence->AddStates(*valence);

    unsigned int total = one_body->CalculateOneElectronIntegrals(core, excited_valence, true);
    total += two_body->CalculateTwoElectronIntegrals(core, valence, excited_valence, excited_valence, true);
    total += two_body->CalculateTwoElectronIntegrals(core, valence, excited_valence, core, true);

    if(!two_body->ReverseSymmetryExists())
    {
        total += two_body->CalculateTwoElectronIntegrals(core, valence, core, excited_valence, true);
        total += two_body->CalculateTwoElectronIntegrals(core, core, excited_valence, excited_valence, true);
    }

    return total;
}

void CoreMBPTCalculator::UpdateIntegrals()
{
    SetValenceEnergies();

    one_body->CalculateOneElectronIntegrals(core, excited);
    one_body->CalculateOneElectronIntegrals(core, valence);
    two_body->CalculateTwoElectronIntegrals(core, valence, excited, excited);
    two_body->CalculateTwoElectronIntegrals(core, valence, excited, core);
    two_body->CalculateTwoElectronIntegrals(core, valence, core, excited);
    two_body->CalculateTwoElectronIntegrals(core, core, excited, valence);
    two_body->CalculateTwoElectronIntegrals(core, core, valence, valence);
}

double CoreMBPTCalculator::GetOneElectronDiagrams(const OrbitalInfo& s1, const OrbitalInfo& s2) const
{
    if(s1.Kappa() != s2.Kappa())
        return 0.;

    if(DebugOptions.LogMBPT())
        *logstream << "\n < " << s1.Name() << " | " << s2.Name() << " > :" << std::endl;

    double term1 = CalculateCorrelation1and3(s1, s2);
    double term2 = CalculateCorrelation2(s1, s2);
    double term4 = CalculateCorrelation4(s1, s2);

    return (term1 + term2 + term4);
}

double CoreMBPTCalculator::GetOneElectronSubtraction(const OrbitalInfo& s1, const OrbitalInfo& s2) const
{
    if(s1.Kappa() != s2.Kappa())
        return 0.;

    if(DebugOptions.LogMBPT())
        *logstream << "\n < " << s1.Name() << " | " << s2.Name() << " > :" << std::endl;

    double term1 = CalculateSubtraction1(s1, s2);
    double term2 = CalculateSubtraction2(s1, s2);
    double term3 = CalculateSubtraction3(s1, s2);

    return (term1 + term2 + term3);
}

double CoreMBPTCalculator::GetTwoElectronDiagrams(unsigned int k, const OrbitalInfo& s1, const OrbitalInfo& s2, const OrbitalInfo& s3, const OrbitalInfo& s4) const
{
    if(DebugOptions.LogMBPT())
        *logstream << "\n R^" << k << " ( " << s1.Name() << " " << s2.Name()
                   << ", " << s3.Name() << " " << s4.Name() << ") :" << std::endl;

    double term = 0.;

    term += CalculateTwoElectron1(k, s1, s2, s3, s4);
    term += CalculateTwoElectron2(k, s1, s2, s3, s4);
    term += CalculateTwoElectron3(k, s1, s2, s3, s4);
    term += CalculateTwoElectron4(k, s1, s2, s3, s4);
    term += CalculateTwoElectron5(k, s1, s2, s3, s4);
    term += CalculateTwoElectron6(k, s1, s2, s3, s4);
    
    return term;
}

double CoreMBPTCalculator::GetTwoElectronBoxDiagrams(unsigned int k, const OrbitalInfo& s1, const OrbitalInfo& s2, const OrbitalInfo& s3, const OrbitalInfo& s4) const
{
    if(DebugOptions.LogMBPT())
        *logstream << "\n R^" << k << " ( " << s1.Name() << " " << s2.Name()
                   << ", " << s3.Name() << " " << s4.Name() << ") :" << std::endl;

    double term = 0.;
    term += CalculateTwoElectron4(k, s1, s2, s3, s4);
    term += CalculateTwoElectron5(k, s1, s2, s3, s4);
    term += CalculateTwoElectron6(k, s1, s2, s3, s4);

    return term;
}

double CoreMBPTCalculator::GetTwoElectronSubtraction(unsigned int k, const OrbitalInfo& s1, const OrbitalInfo& s2, const OrbitalInfo& s3, const OrbitalInfo& s4) const
{
    if(DebugOptions.LogMBPT())
        *logstream << "\n R^" << k << " ( " << s1.Name() << " " << s2.Name()
                   << ", " << s3.Name() << " " << s4.Name() << ") :" << std::endl;

    return CalculateTwoElectronSub(k, s1, s2, s3, s4);
}

double CoreMBPTCalculator::CalculateCorrelation1and3(const OrbitalInfo& sa, const OrbitalInfo& sb) const
{
    const bool debug = DebugOptions.LogMBPT();
    if(debug)
        *logstream << "Cor 1+3:  ";

    const double ValenceEnergy = ValenceEnergies.find(sa.Kappa())->second;
    MathConstant* constants = MathConstant::Instance();

    unsigned int k1, k1max;
    double energy1 = 0., energy3 = 0.;

    // Firstly, get the loop 24
    auto it_n = core->begin();
    while(it_n != core->end())
    {
        const OrbitalInfo& sn = it_n->first;
        const double En = it_n->second->Energy();

        auto it_alpha = excited->begin();
        while(it_alpha != excited->end())
        {
            const OrbitalInfo& salpha = it_alpha->first;
            const double Ealpha = it_alpha->second->Energy();

            k1 = kmin(sn, salpha);
            k1max = kmax(sn, salpha);

            while(k1 <= k1max)
            {
                double C_nalpha = constants->Electron3j(sn.TwoJ(), salpha.TwoJ(), k1);

                if(C_nalpha)
                {
                    C_nalpha = C_nalpha * C_nalpha * sn.MaxNumElectrons() * salpha.MaxNumElectrons()
                                                / (2. * k1 + 1.);

                    // Correlation 1 has excited state beta
                    auto it_beta = excited->begin();
                    while(it_beta != excited->end())
                    {
                        const OrbitalInfo& sbeta = it_beta->first;
                        const double Ebeta = it_beta->second->Energy();

                        double coeff;
                        if(InQSpace(sn, salpha, sbeta) && ParityCheck(sa, sbeta, k1))
                            coeff = constants->Electron3j(sa.TwoJ(), sbeta.TwoJ(), k1);
                        else
                            coeff = 0.;

                        if(coeff)
                        {
                            coeff = coeff * coeff * C_nalpha * sbeta.MaxNumElectrons();
                            double energy_denominator = ValenceEnergy + En - Ebeta - Ealpha + delta;

                            // R1 = R_k1 (a n, beta alpha)
                            // R2 = R_k1 (b n, beta alpha)
                            double R1 = two_body->GetTwoElectronIntegral(k1, sa, sn, sbeta, salpha);
                            double R2 = two_body->GetTwoElectronIntegral(k1, sb, sn, sbeta, salpha);

                            energy1 += TermRatio(R1 * R2 * coeff, energy_denominator, sa, sn, sbeta, salpha);
                        }
                        it_beta++;
                    }

                    // Correlation 3 has core state m
                    auto it_m = core->begin();
                    while(it_m != core->end())
                    {
                        const OrbitalInfo& sm = it_m->first;
                        const double Em = it_m->second->Energy();

                        double coeff;
                        if(InQSpace(sn, salpha, sm) && ParityCheck(sa, sm, k1))
                            coeff =  constants->Electron3j(sa.TwoJ(), sm.TwoJ(), k1);
                        else
                            coeff = 0.;

                        if(coeff)
                        {
                            coeff = coeff * coeff * C_nalpha * sm.MaxNumElectrons();
                            double energy_denominator = ValenceEnergy + Ealpha - En - Em - delta;

                            // R1 = R_k1 (a alpha, m n)
                            // R2 = R_k1 (b alpha, m n)
                            double R1 = two_body->GetTwoElectronIntegral(k1, sa, salpha, sm, sn);
                            double R2 = two_body->GetTwoElectronIntegral(k1, sb, salpha, sm, sn);

                            energy3 += TermRatio(R1 * R2 * coeff, energy_denominator, sa, salpha, sn, sm);
                        }
                        it_m++;
                    }
                }
                k1 += kstep;
            }
            it_alpha++;
        }
        it_n++;
    }

    if(debug)
        *logstream << "  " << std::setprecision(6) << energy1 * constants->HartreeEnergyInInvCm()
                   << "  " << std::setprecision(6) << energy3 * constants->HartreeEnergyInInvCm() << std::endl;
    return energy1 + energy3;
}

double CoreMBPTCalculator::CalculateCorrelation2(const OrbitalInfo& sa, const OrbitalInfo& sb) const
{
    const bool debug = DebugOptions.LogMBPT();

    if(debug)
        *logstream << "Cor 2:    ";

    const double ValenceEnergy = ValenceEnergies.find(sa.Kappa())->second;

    unsigned int k1, k1max;
    unsigned int k2, k2max;

    double energy = 0.;

    auto it_n = core->begin();
    while(it_n != core->end())
    {
        const OrbitalInfo& sn = it_n->first;
        const double En = it_n->second->Energy();

        auto it_alpha = excited->begin();
        while(it_alpha != excited->end())
        {
            const OrbitalInfo& salpha = it_alpha->first;
            const double Ealpha = it_alpha->second->Energy();

            k1 = kmin(sn, salpha);
            k1max = kmax(sn, salpha);

            while(k1 <= k1max)
            {
                double C_nalpha = MathConstant::Instance()->Electron3j(sn.TwoJ(), salpha.TwoJ(), k1);

                if(C_nalpha)
                {
                    C_nalpha = C_nalpha * sn.MaxNumElectrons() * salpha.MaxNumElectrons();

                    auto it_beta = excited->begin();
                    while(it_beta != excited->end())
                    {
                        const OrbitalInfo& sbeta = it_beta->first;
                        const double Ebeta = it_beta->second->Energy();

                        double C_abeta;
                        if(InQSpace(sn, salpha, sbeta) && ParityCheck(sa, sbeta, k1))
                            C_abeta = MathConstant::Instance()->Electron3j(sa.TwoJ(), sbeta.TwoJ(), k1);
                        else
                            C_abeta = 0.;

                        if(C_abeta && ((sa.L() + salpha.L() + sn.L() + sbeta.L())%2 == 0))
                        {
                            C_abeta = C_abeta * sbeta.MaxNumElectrons();
                            double energy_denominator = ValenceEnergy + En - Ebeta - Ealpha + delta;

                            k2 = kmin(sn, sbeta, sa, salpha);
                            k2max = kmax(sn, sbeta, sa, salpha);

                            while(k2 <= k2max)
                            {
                                double coeff
                                    = C_abeta * C_nalpha * MathConstant::Instance()->Electron3j(sa.TwoJ(), salpha.TwoJ(), k2)
                                    * MathConstant::Instance()->Electron3j(sbeta.TwoJ(), sn.TwoJ(), k2)
                                    * MathConstant::Instance()->Wigner6j(sa.J(), sbeta.J(), k1, sn.J(), salpha.J(), k2);
                                    // Note: The 6j symbol is given incorrectly in Berengut et al. PRA 73, 012504 (2006)

                                if(coeff)
                                {
                                    // Sign
                                    if((k1 + k2)%2)
                                        coeff = -coeff;

                                    // R1 = R_k1 (a n, beta alpha)
                                    double R1 = two_body->GetTwoElectronIntegral(k1, sa, sn, sbeta, salpha);

                                    // R2 = R_k2 (beta alpha, n b) = R_k2 (n b, beta alpha)
                                    double R2 = two_body->GetTwoElectronIntegral(k2, sn, sb, sbeta, salpha);

                                    energy += TermRatio(R1 * R2 * coeff, energy_denominator, sa, sn, sbeta, salpha);
                                }
                                k2 += kstep;
                            }
                        }
                        it_beta++;
                    }
                } // C_nalpha
                k1 += kstep;
            }
            it_alpha++;
        }
        it_n++;
    }

    if(debug)
        *logstream << "  " << std::setprecision(6) << energy * MathConstant::Instance()->HartreeEnergyInInvCm() << std::endl;
    return energy;
}

double CoreMBPTCalculator::CalculateCorrelation4(const OrbitalInfo& sa, const OrbitalInfo& sb) const
{
    const bool debug = DebugOptions.LogMBPT();

    if(debug)
        *logstream << "Cor 4:    ";

    const double ValenceEnergy = ValenceEnergies.find(sa.Kappa())->second;

    unsigned int k1, k1max;
    unsigned int k2, k2max;

    double energy = 0.;

    auto it_n = core->begin();
    while(it_n != core->end())
    {
        const OrbitalInfo& sn = it_n->first;
        const double En = it_n->second->Energy();

        auto it_alpha = excited->begin();
        while(it_alpha != excited->end())
        {
            const OrbitalInfo& salpha = it_alpha->first;
            const double Ealpha = it_alpha->second->Energy();

            k1 = kmin(sn, salpha);
            k1max = kmax(sn, salpha);

            while(k1 <= k1max)
            {
                double C_nalpha = MathConstant::Instance()->Electron3j(sn.TwoJ(), salpha.TwoJ(), k1);

                if(C_nalpha)
                {
                    C_nalpha = C_nalpha * sn.MaxNumElectrons() * salpha.MaxNumElectrons();

                    auto it_m = core->begin();
                    while(it_m != core->end())
                    {
                        const OrbitalInfo sm = it_m->first;
                        const double Em = it_m->second->Energy();

                        double C_am;
                        if(InQSpace(sn, salpha, sm) && ParityCheck(sa, sm, k1))
                            C_am = MathConstant::Instance()->Electron3j(sa.TwoJ(), sm.TwoJ(), k1);
                        else
                            C_am = 0.;

                        if(C_am && ((sa.L() + sn.L() + sm.L() + salpha.L())%2 == 0))
                        {
                            C_am = C_am * sm.MaxNumElectrons();
                            double energy_denominator = ValenceEnergy + Ealpha - En - Em - delta;

                            k2 = kmin(sm, salpha, sa, sn);
                            k2max = kmax(sm, salpha, sa, sn);

                            while(k2 <= k2max)
                            {
                                double coeff
                                    = C_am * C_nalpha * MathConstant::Instance()->Electron3j(sa.TwoJ(), sn.TwoJ(), k2)
                                        * MathConstant::Instance()->Electron3j(sm.TwoJ(), salpha.TwoJ(), k2)
                                        * MathConstant::Instance()->Wigner6j(sa.J(), sm.J(), k1, salpha.J(), sn.J(), k2);

                                if(coeff)
                                {   // Sign
                                    if((k1 + k2)%2)
                                        coeff = -coeff;

                                    // R1 = R_k1 (a alpha, m n)
                                    double R1 = two_body->GetTwoElectronIntegral(k1, sa, salpha, sm, sn);

                                    // R2 = R_k2 (m n, alpha b)
                                    double R2 = two_body->GetTwoElectronIntegral(k2, sm, sn, salpha, sb);

                                    energy += TermRatio(R1 * R2 * coeff, energy_denominator, sa, salpha, sn, sm);
                                }
                                k2 += kstep;
                            }
                        }
                        it_m++;
                    }
                } // C_nalpha
                k1 += kstep;
            }
            it_alpha++;
        }
        it_n++;
    }

    if(debug)
        *logstream << "  " << std::setprecision(6) << energy * MathConstant::Instance()->HartreeEnergyInInvCm() << std::endl;
    return energy;
}

double CoreMBPTCalculator::CalculateSubtraction1(const OrbitalInfo& sa, const OrbitalInfo& sb) const
{
    const bool debug = DebugOptions.LogMBPT();

    if(debug)
        *logstream << "Sub 1:    ";

    double energy = 0.;

    auto it_n = core->begin();
    while(it_n != core->end())
    {
        const OrbitalInfo& sn = it_n->first;
        const double En = it_n->second->Energy();

        auto it_alpha = excited->begin();
        while(it_alpha != excited->end())
        {
            const OrbitalInfo& salpha = it_alpha->first;
            const double Ealpha = it_alpha->second->Energy();

            if(sn.Kappa() == salpha.Kappa() && InQSpace(sn, salpha))
            {
                double coeff = one_body->GetMatrixElement(sn, salpha);
                coeff = coeff * sn.MaxNumElectrons();
                double energy_denominator = En - Ealpha + delta;

                // R1 = R_0 (a n, b alpha)
                double R1 = two_body->GetTwoElectronIntegral(0, sa, sn, sb, salpha);

                // Factor of 2 from identical mirror diagram:
                //   R_0 (a n, b alpha) = R_0 (a alpha, b n)
                // and no SMS since k == 0.
                energy += 2. * TermRatio(R1 * coeff, energy_denominator, sn, salpha);
            }
            it_alpha++;
        }
        it_n++;
    }

    if(debug)
        *logstream << "  " << std::setprecision(6) << energy * MathConstant::Instance()->HartreeEnergyInInvCm() << std::endl;
    return energy;
}

double CoreMBPTCalculator::CalculateSubtraction2(const OrbitalInfo& sa, const OrbitalInfo& sb) const
{
    const bool debug = DebugOptions.LogMBPT();
    if(debug)
        *logstream << "Sub 2:    ";

    double energy = 0.;

    unsigned int k1, k1max;

    auto it_n = core->begin();
    while(it_n != core->end())
    {
        const OrbitalInfo& sn = it_n->first;
        const double En = it_n->second->Energy();

        auto it_alpha = excited->begin();
        while(it_alpha != excited->end())
        {
            const OrbitalInfo& salpha = it_alpha->first;
            const double Ealpha = it_alpha->second->Energy();

            if(sn.Kappa() == salpha.Kappa() && InQSpace(sn, salpha))
            {
                double C_nalpha = one_body->GetMatrixElement(sn, salpha);
                C_nalpha = C_nalpha * sn.MaxNumElectrons();
                double energy_denominator = En - Ealpha + delta;

                k1 = kmin(sa, sn);
                k1max = kmax(sa, sn);

                while(k1 <= k1max)
                {
                    double coeff = MathConstant::Instance()->Electron3j(sa.TwoJ(), sn.TwoJ(), k1);
                    coeff = - coeff * coeff * C_nalpha;

                    if(coeff)
                    {
                        // R1 = R_k1 (a alpha, n b), R2 = R_k1 (a n, alpha b)
                        double R1 = two_body->GetTwoElectronIntegral(k1, sa, salpha, sn, sb);
                        double R2 = two_body->GetTwoElectronIntegral(k1, sa, sn, salpha, sb);

                        energy += TermRatio(R1 * coeff, energy_denominator, sn, salpha);
                        energy += TermRatio(R2 * coeff, energy_denominator, sn, salpha);
                    }
                    k1 += kstep;
                }
            }
            it_alpha++;
        }
        it_n++;
    }

    if(debug)
        *logstream << "  " << std::setprecision(6) << energy * MathConstant::Instance()->HartreeEnergyInInvCm() << std::endl;
    return energy;
}

double CoreMBPTCalculator::CalculateSubtraction3(const OrbitalInfo& sa, const OrbitalInfo& sb) const
{
    const bool debug = DebugOptions.LogMBPT();

    if(debug)
        *logstream << "Sub 3:    ";

    double energy = 0.;
    double ValenceEnergy = ValenceEnergies.find(sa.Kappa())->second;

    auto it_n = core->begin();
    while(it_n != core->end())
    {
        const OrbitalInfo& sn = it_n->first;
        const double En = it_n->second->Energy();

        if(sn.Kappa() == sa.Kappa() && InQSpace(sn))
        {
            double term = one_body->GetMatrixElement(sa, sn) * one_body->GetMatrixElement(sn, sb);
            double energy_denominator = En - ValenceEnergy + delta;

            energy -= TermRatio(term, energy_denominator, sn, sa);
        }
        it_n++;
    }

    if(debug)
        *logstream << "  " << std::setprecision(6) << energy * MathConstant::Instance()->HartreeEnergyInInvCm() << std::endl;
    return energy;
}

double CoreMBPTCalculator::CalculateTwoElectron1(unsigned int k, const OrbitalInfo& sa, const OrbitalInfo& sb, const OrbitalInfo& sc, const OrbitalInfo& sd) const
{
    const bool debug = DebugOptions.LogMBPT();

    if(debug)
        *logstream << "TwoE 1:   ";

    double energy = 0.;
    int nn;
#ifdef AMBIT_USE_OPENMP
    #pragma omp parallel for private(nn) reduction(+:energy)
#endif
    /* Note: this needs to be a for loop since OpenMP doesn't gracefully handle reductions over loops 
    with non-random-access iterators */
    for(int nn = 0; nn < core->size(); ++nn)
    {
        auto it_n = core->begin();
        std::advance(it_n, nn);
        const OrbitalInfo& sn = it_n->first;
        const double En = it_n->second->Energy();

        auto it_alpha = excited->begin();
        while(it_alpha != excited->end())
        {
            const OrbitalInfo& salpha = it_alpha->first;
            const double Ealpha = it_alpha->second->Energy();

            double coeff;
            if(InQSpace(sn, salpha) && ParityCheck(sn, salpha, k))
                coeff = MathConstant::Instance()->Electron3j(sn.TwoJ(), salpha.TwoJ(), k);
            else
                coeff = 0.;

            if(coeff)
            {
                coeff = coeff * coeff * sn.MaxNumElectrons() * salpha.MaxNumElectrons()
                                        / (2. * k + 1.);
                double energy_denominator = (En - Ealpha + delta);

                // There are two diagrams:
                //  1. R_k(a n, c alpha) * R_k(alpha b, n d)
                //  2. R_k(a alpha, c n) * R_k(n b, alpha d)
                double R1 = two_body->GetTwoElectronIntegral(k, sa, sn, sc, salpha)
                            * two_body->GetTwoElectronIntegral(k, salpha, sb, sn, sd);
                double R2 = two_body->GetTwoElectronIntegral(k, sa, salpha, sc, sn)
                            * two_body->GetTwoElectronIntegral(k, sn, sb, salpha, sd);

                energy += TermRatio(R1 * coeff, energy_denominator, sn, salpha);
                energy += TermRatio(R2 * coeff, energy_denominator, sn, salpha);
            }
            it_alpha++;
        }
    }

    if(debug)
        *logstream << "  " << std::setprecision(6) << energy * MathConstant::Instance()->HartreeEnergyInInvCm() << std::endl;
    return energy;
}

double CoreMBPTCalculator::CalculateTwoElectron2(unsigned int k, const OrbitalInfo& sa, const OrbitalInfo& sb, const OrbitalInfo& sc, const OrbitalInfo& sd) const
{
    const bool debug = DebugOptions.LogMBPT();

    if(debug)
        *logstream << "TwoE 2/3: ";

    double energy = 0.;
    const double coeff_ac = MathConstant::Instance()->Electron3j(sa.TwoJ(), sc.TwoJ(), k);
    const double coeff_bd = MathConstant::Instance()->Electron3j(sb.TwoJ(), sd.TwoJ(), k);
    if(!coeff_ac || !coeff_bd)
        return energy;

    unsigned int k1, k1max;

    int nn;
#ifdef AMBIT_USE_OPENMP
    #pragma omp parallel for private(nn, k1, k1max) reduction(+:energy)
#endif
    /* Note: this needs to be a for loop since OpenMP doesn't gracefully handle reductions over loops 
    with non-random-access iterators */
    for(nn = 0; nn < core->size(); ++nn)
    {
        auto it_n = core->begin();
        std::advance(it_n, nn);
        const OrbitalInfo& sn = it_n->first;
        const double En = it_n->second->Energy();

        auto it_alpha = excited->begin();
        while(it_alpha != excited->end())
        {
            const OrbitalInfo& salpha = it_alpha->first;
            const double Ealpha = it_alpha->second->Energy();

            double C_nalpha = 0.;
            if(InQSpace(sn, salpha) && ParityCheck(sn, salpha, k))
                C_nalpha = MathConstant::Instance()->Electron3j(sn.TwoJ(), salpha.TwoJ(), k);

            if(C_nalpha)
            {
                C_nalpha = C_nalpha * sn.MaxNumElectrons() * salpha.MaxNumElectrons();
                double energy_denominator = En - Ealpha + delta;

                k1 = kmin(sa, sn, salpha, sc);
                k1max = kmax(sa, sn, salpha, sc);

                while(k1 <= k1max)
                {
                    double coeff = MathConstant::Instance()->Electron3j(sa.TwoJ(), sn.TwoJ(), k1) *
                                   MathConstant::Instance()->Electron3j(salpha.TwoJ(), sc.TwoJ(), k1) *
                                   MathConstant::Instance()->Wigner6j(sa.J(), sc.J(), double(k), salpha.J(), sn.J(), double(k1)) *
                                   C_nalpha / coeff_ac;
                    if((k1 + k)%2)
                        coeff = -coeff;

                    if(coeff)
                    {
                        // R1 = R_k1 (a alpha, n c)
                        double R1 = two_body->GetTwoElectronIntegral(k1, sa, salpha, sn, sc);

                        // R2 = R_k (n b, alpha d)
                        double R2 = two_body->GetTwoElectronIntegral(k, sn, sb, salpha, sd);

                        energy += TermRatio(R1 * R2 * coeff, energy_denominator, sn, salpha);
                    }
                    k1 += kstep;
                }

                // Mirror diagram
                k1 = kmin(sb, sn, salpha, sd);
                k1max = kmax(sb, sn, salpha, sd);

                while(k1 <= k1max)
                {
                    double coeff = MathConstant::Instance()->Electron3j(sb.TwoJ(), sn.TwoJ(), k1) *
                                   MathConstant::Instance()->Electron3j(salpha.TwoJ(), sd.TwoJ(), k1) *
                                   MathConstant::Instance()->Wigner6j(sb.J(), sd.J(), double(k), salpha.J(), sn.J(), double(k1)) *
                                   C_nalpha / coeff_bd;
                    if((k1 + k)%2)
                        coeff = -coeff;

                    if(coeff)
                    {
                        // R1 = R_k1 (b alpha, n d)
                        double R1 = two_body->GetTwoElectronIntegral(k1, sb, salpha, sn, sd);

                        // R2 = R_k (a n, c alpha)
                        double R2 = two_body->GetTwoElectronIntegral(k, sa, sn, sc, salpha);

                        energy += TermRatio(R1 * R2 * coeff, energy_denominator, sn, salpha);
                    }
                    k1 += kstep;
                }
            }
            it_alpha++;
        }
    }

    if(debug)
        *logstream << "  " << std::setprecision(6) << energy * MathConstant::Instance()->HartreeEnergyInInvCm() << std::endl;
    return energy;
}

double CoreMBPTCalculator::CalculateTwoElectron3(unsigned int k, const OrbitalInfo& sa, const OrbitalInfo& sb, const OrbitalInfo& sc, const OrbitalInfo& sd) const
{
    return CalculateTwoElectron2(k, sc, sd, sa, sb);
}

double CoreMBPTCalculator::CalculateTwoElectron4(unsigned int k, const OrbitalInfo& sa, const OrbitalInfo& sb, const OrbitalInfo& sc, const OrbitalInfo& sd) const
{
    const bool debug = DebugOptions.LogMBPT();

    if(debug)
        *logstream << "TwoE 4/5: ";

    double energy = 0.;
    const double coeff_ac = MathConstant::Instance()->Electron3j(sa.TwoJ(), sc.TwoJ(), k);
    const double coeff_bd = MathConstant::Instance()->Electron3j(sb.TwoJ(), sd.TwoJ(), k);
    if(!coeff_ac || !coeff_bd)
        return energy;

    unsigned int k1, k1max;
    unsigned int k2, k2max;

    int nn;
#ifdef AMBIT_USE_OPENMP
    #pragma omp parallel for private(nn, k1, k1max, k2, k2max) reduction(+:energy)
#endif
    /* Note: this needs to be a for loop since OpenMP doesn't gracefully handle reductions over loops 
    with non-random-access iterators */
    for(nn = 0; nn < core->size(); ++nn)
    {
        auto it_n = core->begin();
        std::advance(it_n, nn);

        const OrbitalInfo& sn = it_n->first;
        const double En = it_n->second->Energy();

        auto it_alpha = excited->begin();
        while(it_alpha != excited->end())
        {
            const OrbitalInfo& salpha = it_alpha->first;
            const double Ealpha = it_alpha->second->Energy();

            double C_nalpha = 0.;
            if(InQSpace(sn, salpha) &&
               ((sa.L() + sn.L() + salpha.L() + sd.L())%2 == 0) &&
               ((sn.L() + sc.L() + sb.L() + salpha.L())%2 == 0))
                C_nalpha = double(sn.MaxNumElectrons()) * double(salpha.MaxNumElectrons()) * (2. * double(k) + 1.);

            if(C_nalpha)
            {
                C_nalpha = C_nalpha/(coeff_ac*coeff_bd);
                double energy_denominator = En - Ealpha + delta;

                unsigned int phase = (unsigned int)(sa.TwoJ() + sb.TwoJ() + sc.TwoJ() + sd.TwoJ() + sn.TwoJ() + salpha.TwoJ())/2;
                if(phase%2)
                    C_nalpha = -C_nalpha;

                k1 = kmin(sa, sn, salpha, sd);
                k1max = kmax(sa, sn, salpha, sd);

                while(k1 <= k1max)
                {
                    double coeff_ad = MathConstant::Instance()->Electron3j(sa.TwoJ(), sn.TwoJ(), k1) *
                                      MathConstant::Instance()->Electron3j(salpha.TwoJ(), sd.TwoJ(), k1);

                    if(coeff_ad)
                    {
                        // R1 = R_k1 (a alpha, n d)
                        double R1 = two_body->GetTwoElectronIntegral(k1, sa, salpha, sn, sd);

                        k2 = kmin(sn, sc, sb, salpha);
                        k2max = kmax(sn, sc, sb, salpha);

                        while(k2 <= k2max)
                        {
                            double coeff = MathConstant::Instance()->Electron3j(sn.TwoJ(), sc.TwoJ(), k2) *
                                           MathConstant::Instance()->Electron3j(sb.TwoJ(), salpha.TwoJ(), k2);
                            if(coeff)
                                coeff = coeff * MathConstant::Instance()->Wigner6j(sc.J(), sa.J(), k, k1, k2, sn.J())
                                              * MathConstant::Instance()->Wigner6j(sb.J(), sd.J(), k, k1, k2, salpha.J());

                            if(coeff)
                            {
                                coeff = coeff * coeff_ad * C_nalpha;

                                // R2 = R_k2 (2b, c4)
                                double R2 = two_body->GetTwoElectronIntegral(k2, sn, sb, sc, salpha);

                                energy += TermRatio(R1 * R2 * coeff, energy_denominator, sn, salpha);
                            }
                            k2 += kstep;
                        }
                    }
                    k1 += kstep;
                }
            }
            it_alpha++;
        }
    }

    if(debug)
        *logstream << "  " << std::setprecision(6) << energy * MathConstant::Instance()->HartreeEnergyInInvCm() << std::endl;
    return energy;
}

double CoreMBPTCalculator::CalculateTwoElectron5(unsigned int k, const OrbitalInfo& sa, const OrbitalInfo& sb, const OrbitalInfo& sc, const OrbitalInfo& sd) const
{
    return CalculateTwoElectron4(k, sb, sa, sd, sc);
}

double CoreMBPTCalculator::CalculateTwoElectron6(unsigned int k, const OrbitalInfo& sa, const OrbitalInfo& sb, const OrbitalInfo& sc, const OrbitalInfo& sd) const
{
    const bool debug = DebugOptions.LogMBPT();

    if(debug)
        *logstream << "TwoE 6:   ";

    double energy = 0.;
    const double coeff_ac = MathConstant::Instance()->Electron3j(sa.TwoJ(), sc.TwoJ(), k);
    const double coeff_bd = MathConstant::Instance()->Electron3j(sb.TwoJ(), sd.TwoJ(), k);
    if(!coeff_ac || !coeff_bd)
        return energy;

    const double ValenceEnergy = ValenceEnergies.find(sa.Kappa())->second + ValenceEnergies.find(sb.Kappa())->second;

    unsigned int k1, k1max;
    unsigned int k2, k2max;

    int mm;
#ifdef AMBIT_USE_OPENMP
    #pragma omp parallel for private(mm, k1, k1max, k2, k2max) reduction(+:energy)
#endif
    /* Note: this needs to be a for loop since OpenMP doesn't gracefully handle reductions over loops 
    with non-random-access iterators */
    for(mm = 0; mm < core->size(); ++mm)
    {
        auto it_m = core->begin();
        std::advance(it_m, mm);
        const OrbitalInfo& sm = it_m->first;
        const double Em = it_m->second->Energy();

        auto it_n = core->begin();
        while(it_n != core->end())
        {
            const OrbitalInfo& sn = it_n->first;
            const double En = it_n->second->Energy();

            double coeff_mn = 0.;
            if(InQSpace(sn, sm) &&
               ((sa.L() + sm.L() + sb.L() + sn.L())%2 == 0) &&
               ((sm.L() + sc.L() + sd.L() + sn.L())%2 == 0))
                coeff_mn = double(sm.MaxNumElectrons()) * double(sn.MaxNumElectrons()) * (2. * double(k) + 1.);

            if(coeff_mn)
            {
                coeff_mn = coeff_mn/(coeff_ac*coeff_bd);
                double energy_denominator = Em + En - ValenceEnergy + delta;

                unsigned int phase = (unsigned int)(sa.TwoJ() + sb.TwoJ() + sc.TwoJ() + sd.TwoJ() + sm.TwoJ() + sn.TwoJ())/2;
                if((phase + k + 1)%2)
                    coeff_mn = -coeff_mn;

                k1 = kmin(sa, sm, sb, sn);
                k1max = kmax(sa, sm, sb, sn);

                while(k1 <= k1max)
                {
                    double coeff_ab = MathConstant::Instance()->Electron3j(sa.TwoJ(), sm.TwoJ(), k1) *
                                      MathConstant::Instance()->Electron3j(sb.TwoJ(), sn.TwoJ(), k1);

                    if(coeff_ab)
                    {
                        // R1 = R_k1 (ab, mn)
                        double R1 = two_body->GetTwoElectronIntegral(k1, sa, sb, sm, sn);

                        k2 = kmin(sm, sc, sn, sd);
                        k2max = kmax(sm, sc, sn, sd);

                        while(k2 <= k2max)
                        {
                            double coeff = MathConstant::Instance()->Electron3j(sm.TwoJ(), sc.TwoJ(), k2) *
                                           MathConstant::Instance()->Electron3j(sn.TwoJ(), sd.TwoJ(), k2);
                            if(coeff)
                                coeff = coeff * MathConstant::Instance()->Wigner6j(sc.J(), sa.J(), k, k1, k2, sm.J())
                                              * MathConstant::Instance()->Wigner6j(sd.J(), sb.J(), k, k1, k2, sn.J());

                            if(coeff)
                            {
                                coeff = coeff * coeff_ab * coeff_mn;
                                if((k1 + k2)%2)
                                    coeff = -coeff;

                                // R2 = R_k2 (mn, cd)
                                double R2 = two_body->GetTwoElectronIntegral(k2, sm, sn, sc, sd);

                                energy += TermRatio(R1 * R2 * coeff, energy_denominator, sm, sn, sa, sb);
                            }
                            k2 += kstep;
                        }
                    }
                    k1 += kstep;
                }
            }
            it_n++;
        }
    }

    if(debug)
        *logstream << "  " << std::setprecision(6) << energy * MathConstant::Instance()->HartreeEnergyInInvCm() << std::endl;
    return energy;
}

double CoreMBPTCalculator::CalculateTwoElectronSub(unsigned int k, const OrbitalInfo& sa, const OrbitalInfo& sb, const OrbitalInfo& sc, const OrbitalInfo& sd) const
{
    const bool debug = DebugOptions.LogMBPT();

    if(debug)
        *logstream << "2eSub :   ";

    double energy = 0.;

    const double Ea = ValenceEnergies.find(sa.Kappa())->second;
    const double Eb = ValenceEnergies.find(sb.Kappa())->second;
    const double Ec = ValenceEnergies.find(sc.Kappa())->second;
    const double Ed = ValenceEnergies.find(sd.Kappa())->second;

    // Hole line is attached to sa or sc
    int nn;
#ifdef AMBIT_USE_OPENMP
    #pragma omp parallel for private(nn) reduction(+:energy)
#endif
    /* Note: this needs to be a for loop since OpenMP doesn't gracefully handle reductions over loops 
    with non-random-access iterators */
    for(nn = 0; nn < core->size(); ++nn)
    {
        auto it_n = core->begin();
        std::advance(it_n, nn);
        const OrbitalInfo& sn = it_n->first;
        if(InQSpace(sn))
        {
            const double En = it_n->second->Energy();

            if(sn.Kappa() == sa.Kappa())
            {
                double R1 = two_body->GetTwoElectronIntegral(k, sn, sb, sc, sd);
                double R2 = one_body->GetMatrixElement(sa, sn);
                double energy_denominator = (En - Ea + delta);

                energy -= TermRatio(R1 * R2, energy_denominator, sn, sa);
            }

            if(sn.Kappa() == sc.Kappa())
            {
                double R1 = two_body->GetTwoElectronIntegral(k, sa, sb, sn, sd);
                double R2 = one_body->GetMatrixElement(sn, sc);
                double energy_denominator = (En - Ec + delta);

                energy -= TermRatio(R1 * R2, energy_denominator, sn, sc);
            }

            if(sn.Kappa() == sb.Kappa())
            {
                double R1 = two_body->GetTwoElectronIntegral(k, sa, sn, sc, sd);
                double R2 = one_body->GetMatrixElement(sb, sn);
                double energy_denominator = (En - Eb + delta);

                energy -= TermRatio(R1 * R2, energy_denominator, sn, sb);
            }

            if(sn.Kappa() == sd.Kappa())
            {
                double R1 = two_body->GetTwoElectronIntegral(k, sa, sb, sc, sn);
                double R2 = one_body->GetMatrixElement(sn, sd);
                double energy_denominator = (En - Ed + delta);

                energy -= TermRatio(R1 * R2, energy_denominator, sn, sd);
            }
        }
    }

    if(debug)
        *logstream << "  " << std::setprecision(6) << energy * MathConstant::Instance()->HartreeEnergyInInvCm() << std::endl;
    return energy;
}

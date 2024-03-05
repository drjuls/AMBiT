#include "Include.h"
#include "ValenceMBPTCalculator.h"
#include "Universal/MathConstant.h"
#include "Universal/PhysicalConstant.h"

#ifdef AMBIT_USE_OPENMP
#include <omp.h>
#endif

namespace Ambit
{
ValenceMBPTCalculator::ValenceMBPTCalculator(pOrbitalManagerConst orbitals, pHFIntegrals one_body, pSlaterIntegrals two_body, const std::string& fermi_orbitals):
    MBPTCalculator(orbitals, fermi_orbitals, two_body->OffParityExists()), one_body(one_body), two_body(two_body), excited(orbitals->excited), high(orbitals->high)
{}

ValenceMBPTCalculator::~ValenceMBPTCalculator()
{
    one_body->clear();
    two_body->clear();
}

unsigned int ValenceMBPTCalculator::GetStorageSize()
{
    unsigned int total = one_body->CalculateOneElectronIntegrals(high, valence, true);
    total += two_body->CalculateTwoElectronIntegrals(valence, valence, excited, high, true);
    total += two_body->CalculateTwoElectronIntegrals(valence, valence, orbitals->hole, high, true);

    return total;
}

void ValenceMBPTCalculator::UpdateIntegrals()
{
    SetValenceEnergies();

    one_body->CalculateOneElectronIntegrals(high, valence);
    two_body->CalculateTwoElectronIntegrals(valence, valence, excited, high);
    two_body->CalculateTwoElectronIntegrals(valence, valence, valence, high);
}

double ValenceMBPTCalculator::GetOneElectronSubtraction(const OrbitalInfo& s1, const OrbitalInfo& s2)
{
    if(s1.Kappa() != s2.Kappa())
        return 0.;

    return CalculateOneElectronSub(s1, s2);
}

double ValenceMBPTCalculator::GetTwoElectronValence(unsigned int k, const OrbitalInfo& s1, const OrbitalInfo& s2, const OrbitalInfo& s3, const OrbitalInfo& s4)
{
    return CalculateTwoElectronValence(k, s1, s2, s3, s4);
}

double ValenceMBPTCalculator::GetTwoElectronSubtraction(unsigned int k, const OrbitalInfo& s1, const OrbitalInfo& s2, const OrbitalInfo& s3, const OrbitalInfo& s4)
{
    return CalculateTwoElectronSub(k, s1, s2, s3, s4);
}

double ValenceMBPTCalculator::GetTwoElectronBoxValence(unsigned int k, const OrbitalInfo& s1, const OrbitalInfo& s2, const OrbitalInfo& s3, const OrbitalInfo& s4)
{
    return CalculateTwoElectronValence(k, s1, s2, s3, s4);
}

double ValenceMBPTCalculator::CalculateOneElectronSub(const OrbitalInfo& sa, const OrbitalInfo& sb) const
{
    const bool debug = DebugOptions.LogMBPT();
    if(debug)
        *logstream << "Val 1.1:    ";

    double energy = 0.;
    const double ValenceEnergy = ValenceEnergies.find(sa.Kappa())->second;

    auto it_alpha = high->begin();
    while(it_alpha != high->end())
    {
        const OrbitalInfo& salpha = it_alpha->first;
        const double Ealpha = it_alpha->second->Energy();

        // InQSpace() usually will be true for states in "high" but test just in case this is overridden
        if(InQSpace(salpha) && (sa.Kappa() == salpha.Kappa()))
        {
            double term = one_body->GetMatrixElement(sa, salpha) * one_body->GetMatrixElement(salpha, sb);
            double energy_denominator = ValenceEnergy - Ealpha + delta;
            energy += TermRatio(term, energy_denominator, sa, salpha);
        }

        it_alpha++;
    }

    if(debug)
        *logstream << "  " << std::setprecision(6) << energy * MathConstant::Instance()->HartreeEnergyInInvCm() << std::endl;
    return energy;
}

double ValenceMBPTCalculator::CalculateTwoElectronValence(unsigned int k, const OrbitalInfo& sa, const OrbitalInfo& sb, const OrbitalInfo& sc, const OrbitalInfo& sd) const
{
    const bool debug = DebugOptions.LogMBPT();
    if(debug)
        *logstream << "Val 2.1:   ";

    MathConstant* constants = MathConstant::Instance();

    double energy = 0.;
    const double coeff_ac = constants->Electron3j(sa.TwoJ(), sc.TwoJ(), k);
    const double coeff_bd = constants->Electron3j(sb.TwoJ(), sd.TwoJ(), k);
    if(!coeff_ac || !coeff_bd)
        return energy;

    const double ValenceEnergy = ValenceEnergies.find(sa.Kappa())->second + ValenceEnergies.find(sb.Kappa())->second;
    int k1, k1max;
    int k2, k2max;

    int ii;
    for(ii = 0; ii < excited->size(); ++ii)
    {
        auto it_alpha = excited->begin();
        std::advance(it_alpha, ii);
        const OrbitalInfo& salpha = it_alpha->first;
        const double Ealpha = it_alpha->second->Energy();

        auto it_beta = excited->begin();
        while(it_beta != excited->end())
        {
            const OrbitalInfo& sbeta = it_beta->first;
            const double Ebeta = it_beta->second->Energy();

            if(InQSpace(salpha, sbeta) &&
               ((sa.L() + salpha.L() + sb.L() + sbeta.L())%2 == 0) &&
               ((salpha.L() + sc.L() + sbeta.L() + sd.L())%2 == 0))
            {
                double coeff_alphabeta = salpha.MaxNumElectrons() * sbeta.MaxNumElectrons() * (2 * k + 1);

                coeff_alphabeta = coeff_alphabeta/(coeff_ac * coeff_bd);
                double energy_denominator = (ValenceEnergy - Ealpha - Ebeta + delta);

                int exponent = (sa.TwoJ() + sb.TwoJ() + sc.TwoJ() + sd.TwoJ() + salpha.TwoJ() + sbeta.TwoJ())/2;
                coeff_alphabeta *= constants->minus_one_to_the_power(exponent + k + 1);

                k1 = kmin(sa, salpha);
                k1max = kmax(sa, salpha);

                while(k1 <= k1max)
                {
                    double coeff_ab = constants->Electron3j(sa.TwoJ(), salpha.TwoJ(), k1) *
                                      constants->Electron3j(sb.TwoJ(), sbeta.TwoJ(), k1);

                    if(coeff_ab)
                    {
                        double R1 = two_body->GetTwoElectronIntegral(k1, sa, sb, salpha, sbeta);

                        k2 = kmin(salpha, sc);
                        k2max = kmax(salpha, sc);

                        while(k2 <= k2max)
                        {
                            double coeff = constants->Electron3j(salpha.TwoJ(), sc.TwoJ(), k2) *
                                           constants->Electron3j(sbeta.TwoJ(), sd.TwoJ(), k2);

                            if(coeff)
                                coeff = coeff * constants->Wigner6j(sc.J(), sa.J(), k, k1, k2, salpha.J())
                                              * constants->Wigner6j(sd.J(), sb.J(), k, k1, k2, sbeta.J());
                            if(coeff)
                            {
                                coeff = coeff * coeff_ab * coeff_alphabeta;
                                if((k1 + k2)%2)
                                    coeff = -coeff;

                                double R2 = two_body->GetTwoElectronIntegral(k2, salpha, sbeta, sc, sd);
                                double numerator = R1 * R2 * coeff;

                                energy += TermRatio(numerator, energy_denominator, sa, sb, salpha, sbeta);
                            }
                            k2 += kstep;
                        }
                    }
                    k1 += kstep;
                }
            }
            it_beta++;
        }
    }

    if(debug)
        *logstream << "  " << std::setprecision(6) << energy * constants->HartreeEnergyInInvCm() << std::endl;
    return energy;
}

double ValenceMBPTCalculator::CalculateTwoElectronSub(unsigned int k, const OrbitalInfo& sa, const OrbitalInfo& sb, const OrbitalInfo& sc, const OrbitalInfo& sd) const
{
    const bool debug = DebugOptions.LogMBPT();
    if(debug)
        *logstream << "Val 2.2:   ";

    double energy = 0.;

    const double Ea = ValenceEnergies.find(sa.Kappa())->second;
    const double Eb = ValenceEnergies.find(sb.Kappa())->second;
    const double Ec = ValenceEnergies.find(sc.Kappa())->second;
    const double Ed = ValenceEnergies.find(sd.Kappa())->second;

    int ii;
    for(ii=0; ii < high->size(); ++ii)
    {
        auto it_alpha = high->begin();
        std::advance(it_alpha, ii);
        const OrbitalInfo& salpha = it_alpha->first;
        const double Ealpha = it_alpha->second->Energy();

        // InQSpace() usually will be true for states in "high" but test just in case this is overridden
        if(InQSpace(salpha))
        {
            if(sa.Kappa() == salpha.Kappa())
            {
                double R1 = two_body->GetTwoElectronIntegral(k, salpha, sb, sc, sd);
                double energy_denominator = Ea - Ealpha + delta;
                energy += TermRatio(R1 * one_body->GetMatrixElement(sa, salpha), energy_denominator, sa, salpha);
            }

            if(sc.Kappa() == salpha.Kappa())
            {
                double R1 = two_body->GetTwoElectronIntegral(k, sa, sb, salpha, sd);
                double energy_denominator = Ec - Ealpha + delta;
                energy += TermRatio(R1 * one_body->GetMatrixElement(salpha, sc), energy_denominator, sc, salpha);
            }

            if(sb.Kappa() == salpha.Kappa())
            {
                double R1 = two_body->GetTwoElectronIntegral(k, sa, salpha, sc, sd);
                double energy_denominator = Eb - Ealpha + delta;
                energy += TermRatio(R1 * one_body->GetMatrixElement(sb, salpha), energy_denominator, sb, salpha);
            }

            if(sd.Kappa() == salpha.Kappa())
            {
                double R1 = two_body->GetTwoElectronIntegral(k, sa, sb, sc, salpha);
                double energy_denominator = Ed - Ealpha + delta;
                energy += TermRatio(R1 * one_body->GetMatrixElement(salpha, sd), energy_denominator, sd, salpha);
            }
        }
    }

    if(debug)
        *logstream << "  " << std::setprecision(6) << energy * MathConstant::Instance()->HartreeEnergyInInvCm() << std::endl;
    return energy;
}
}

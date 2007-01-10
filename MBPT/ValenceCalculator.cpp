#include "Include.h"
#include "ValenceCalculator.h"
#include "Universal/Constant.h"
#include "Universal/CoulombIntegrator.h"
#include "HartreeFock/StateIntegrator.h"

#define MAX_K 12

ValenceCalculator::ValenceCalculator(Lattice* lattice, const Core* atom_core, const ExcitedStates* excited_states):
    MBPTCalculator(lattice, atom_core, excited_states)
{
    UseBrillouinWignerPT();
}

double ValenceCalculator::GetOneElectronValence(const State* s1, const State* s2)
{
    MaxStateSize = core->GetConstHFPotential().size();

    if(s1->Kappa() != s2->Kappa())
        return 0.;

    return CalculateOneElectronValence1(*s1, *s2);
}

double ValenceCalculator::GetTwoElectronValence(const State* s1, const State* s2, const State* s3, const State* s4, unsigned int k)
{
    MaxStateSize = core->GetConstHFPotential().size();

    double term = 0;
    term += CalculateTwoElectronValence1(*s1, *s2, *s3, *s4, k);
    term += CalculateTwoElectronValence2(*s1, *s2, *s3, *s4, k);

    return term;
}

double ValenceCalculator::GetTwoElectronBoxValence(const State* s1, const State* s2, const State* s3, const State* s4, unsigned int k)
{
    MaxStateSize = core->GetConstHFPotential().size();

    return CalculateTwoElectronValence1(*s1, *s2, *s3, *s4, k);
}

double ValenceCalculator::CalculateOneElectronValence1(const State& si, const State& sf) const
{
    const bool debug = DebugOptions.LogMBPT();
    StateIntegrator SI(lattice);

    if(debug)
        *outstream << "Val 1.1:    ";

    double energy = 0.;
    const double ValenceEnergy = ValenceEnergies.find(si.Kappa())->second;

    ConstStateIterator it4 = excited->GetConstStateIterator();
    while(!it4.AtEnd())
    {   const DiscreteState& s4 = *(it4.GetState());

        if(s4.RequiredPQN() >= 5 && s4.Kappa() == si.Kappa())
        {
            double term = SI.HamiltonianMatrixElement(s4, si, *core) * SI.HamiltonianMatrixElement(s4, sf, *core);
            if(BrillouinWignerPT)
                term = term/(ValenceEnergy - s4.Energy() + delta);
            else
                term = term/(si.Energy() - s4.Energy());
                
            energy = energy + term;
        }
        it4.Next();
    }

    if(debug)
        *outstream << "  " << energy * Constant::HartreeEnergy_cm << std::endl;
    return energy;
}

double ValenceCalculator::CalculateTwoElectronValence1(const State& sa, const State& sb, const State& sc, const State& sd, unsigned int k) const
{
    const bool debug = DebugOptions.LogMBPT();
    const double NuclearInverseMass = core->GetNuclearInverseMass();

    unsigned int Ja = (unsigned int)(sa.J() * 2.);
    unsigned int Jb = (unsigned int)(sb.J() * 2.);
    unsigned int Jc = (unsigned int)(sc.J() * 2.);
    unsigned int Jd = (unsigned int)(sd.J() * 2.);

    std::vector<double> density(MaxStateSize);
    std::vector<double> pot(MaxStateSize);

    CoulombIntegrator I(lattice);
    StateIntegrator SI(lattice);
    const double* dR = lattice->dR();

    if(debug)
        *outstream << "Val 2.1:   ";
    double spacing = 1./double(excited->NumStates() * excited->NumStates());
    double count = 0.;
    unsigned int i;

    double energy = 0.;
    const double coeff_ac = Constant::Electron3j(Ja, Jc, k, 1, -1);
    const double coeff_bd = Constant::Electron3j(Jb, Jd, k, 1, -1);
    if(!coeff_ac || !coeff_bd)
        return energy;

    const double ValenceEnergy = ValenceEnergies.find(sa.Kappa())->second + ValenceEnergies.find(sb.Kappa())->second;

    ConstStateIterator it3 = excited->GetConstStateIterator();
    while(!it3.AtEnd())
    {   const DiscreteState& s3 = *(it3.GetState());
        unsigned int J3 = (unsigned int)(s3.J() * 2.);

        ConstStateIterator it4 = excited->GetConstStateIterator();
        while(!it4.AtEnd())
        {   const DiscreteState& s4 = *(it4.GetState());
            unsigned int J4 = (unsigned int)(s4.J() * 2.);

            if(debug)
            {   count += spacing;
                if(count >= 0.02)
                {   *logstream << ".";
                    count -= 0.02;
                }
            }

            double coeff_34 = double(J3 + 1) * double(J4 + 1) * (2. * double(k) + 1.);
            if((sa.L() + s3.L())%2 != (sb.L() + s4.L())%2)
                coeff_34 = 0.;
            if((s3.L() + sc.L())%2 != (sd.L() + s4.L())%2)
                coeff_34 = 0.;

            if((s3.RequiredPQN() >= 5 || s4.RequiredPQN() >= 5) && coeff_34)
            {
                coeff_34 = coeff_34/(coeff_ac*coeff_bd);
                coeff_34 = coeff_34/(ValenceEnergy - s3.Energy() - s4.Energy() + delta);
                unsigned int exponent = (unsigned int)(Ja + Jb + Jc + Jd + J3 + J4)/2;
                if((exponent + k + 1)%2)
                    coeff_34 = -coeff_34;

                unsigned int start_k1 = (sa.L() + s3.L())%2;
                for(unsigned int k1 = start_k1; k1 <= MAX_K; k1+=2)
                {
                    double coeff_ab = Constant::Electron3j(Ja, J3, k1, 1, -1) *
                                      Constant::Electron3j(Jb, J4, k1, 1, -1);

                    if(coeff_ab)
                    {
                        // R1 = R_k1 (ab, mn)
                        double R1 = 0.;
                        for(i=0; i<mmin(sb.Size(), s4.Size()); i++)
                        {
                            density[i] = sb.f[i] * s4.f[i] + Constant::AlphaSquared * sb.g[i] * s4.g[i];
                        }
                        I.FastCoulombIntegrate(density, pot, k1, mmin(sb.Size(), s4.Size()));
                        
                        for(i=0; i < mmin(sa.Size(), s3.Size()); i++)
                            R1 = R1 + pot[i] * (sa.f[i] * s3.f[i] + Constant::AlphaSquared * sa.g[i] * s3.g[i]) * dR[i];

                        if(NuclearInverseMass && (k1 == 1))
                        {   R1 = R1 - NuclearInverseMass * SI.IsotopeShiftIntegral(s3, sa) * SI.IsotopeShiftIntegral(s4, sb);
                        }

                        unsigned int start_k2 = (s3.L() + sc.L())%2;
                        for(unsigned int k2 = start_k2; k2 <= MAX_K; k2+=2)
                        {
                            double coeff = Constant::Electron3j(J3, Jc, k2, 1, -1) *
                                           Constant::Electron3j(J4, Jd, k2, 1, -1);
                            if(coeff)
                                coeff = coeff * Constant::Wigner6j(sc.J(), sa.J(), k, k1, k2, s3.J())
                                              * Constant::Wigner6j(sd.J(), sb.J(), k, k1, k2, s4.J());

                            if(coeff)
                            {   
                                coeff = coeff * coeff_ab * coeff_34;
                                if((k1 + k2)%2)
                                    coeff = -coeff;

                                // R2 = R_k2 (mn, cd)
                                double R2 = 0.;
                                for(i=0; i<mmin(s4.Size(), sd.Size()); i++)
                                {
                                    density[i] = s4.f[i] * sd.f[i] + Constant::AlphaSquared * s4.g[i] * sd.g[i];
                                }
                                I.FastCoulombIntegrate(density, pot, k2, mmin(s4.Size(), sd.Size()));
                                
                                for(i=0; i < mmin(s3.Size(), sc.Size()); i++)
                                    R2 = R2 + pot[i] * (s3.f[i] * sc.f[i] + Constant::AlphaSquared * s3.g[i] * sc.g[i]) * dR[i];

                                if(NuclearInverseMass && (k2 == 1))
                                {   R2 = R2 - NuclearInverseMass * SI.IsotopeShiftIntegral(s3, sc) * SI.IsotopeShiftIntegral(s4, sd);
                                }

                                energy += R1 * R2 * coeff;
                            }
                        }
                    }
                }
            }
            it4.Next();
        }
        it3.Next();
    }

    if(debug)
        *outstream << "  " << energy * Constant::HartreeEnergy_cm << std::endl;
    return energy;
}

double ValenceCalculator::CalculateTwoElectronValence2(const State& sa, const State& sb, const State& sc, const State& sd, unsigned int k) const
{
    const bool debug = DebugOptions.LogMBPT();
    const double NuclearInverseMass = core->GetNuclearInverseMass();

    std::vector<double> density(MaxStateSize);
    std::vector<double> pot(MaxStateSize);

    CoulombIntegrator I(lattice);
    StateIntegrator SI(lattice);
    const double* dR = lattice->dR();

    if(debug)
        *outstream << "Val 2.2:   ";

    unsigned int i;
    double energy = 0.;

    // Hole line is attached to sa or sc
    for(i=0; i<mmin(sb.Size(), sd.Size()); i++)
    {
        density[i] = sb.f[i] * sd.f[i] + Constant::AlphaSquared * sb.g[i] * sd.g[i];
    }
    I.FastCoulombIntegrate(density, pot, k, mmin(sb.Size(), sd.Size()));

    double SMS_bd = 0.;
    if(NuclearInverseMass && (k == 1))
    {   SMS_bd = NuclearInverseMass * SI.IsotopeShiftIntegral(sb, sd);
    }

    ConstStateIterator it3 = excited->GetConstStateIterator();
    while(!it3.AtEnd())
    {   const DiscreteState& s3 = *(it3.GetState());

        if(s3.RequiredPQN() >= 5 && s3.Kappa() == sa.Kappa())
        {
            double R1 = 0.;
            for(i=0; i < mmin(s3.Size(), sc.Size()); i++)
                R1 += pot[i] * (s3.f[i] * sc.f[i] + Constant::AlphaSquared * s3.g[i] * sc.g[i]) * dR[i];

            if(SMS_bd)
            {   R1 = R1 - SMS_bd * SI.IsotopeShiftIntegral(s3, sc);
            }

            energy += R1 * SI.HamiltonianMatrixElement(sa, s3, *core) / (ValenceEnergies.find(sa.Kappa())->second - s3.Energy() + delta);
        }

        if(s3.RequiredPQN() >= 5 && s3.Kappa() == sc.Kappa())
        {
            double R1 = 0.;
            for(i=0; i < mmin(sa.Size(), s3.Size()); i++)
                R1 += pot[i] * (sa.f[i] * s3.f[i] + Constant::AlphaSquared * sa.g[i] * s3.g[i]) * dR[i];

            if(SMS_bd)
            {   R1 = R1 + SMS_bd * SI.IsotopeShiftIntegral(s3, sa);
            }

            energy += R1 * SI.HamiltonianMatrixElement(sc, s3, *core) / (ValenceEnergies.find(sc.Kappa())->second - s3.Energy() + delta);
        }

        it3.Next();
    }

    // Hole line is attached to sb or sd.
    for(i=0; i<mmin(sa.Size(), sc.Size()); i++)
    {
        density[i] = sa.f[i] * sc.f[i] + Constant::AlphaSquared * sa.g[i] * sc.g[i];
    }
    I.FastCoulombIntegrate(density, pot, k, mmin(sa.Size(), sc.Size()));

    double SMS_ac = 0.;
    if(NuclearInverseMass && (k == 1))
    {   SMS_ac = NuclearInverseMass * SI.IsotopeShiftIntegral(sa, sc);
    }

    it3.First();
    while(!it3.AtEnd())
    {   const DiscreteState& s3 = *(it3.GetState());

        if(s3.RequiredPQN() >= 5 && s3.Kappa() == sb.Kappa())
        {
            double R1 = 0.;
            for(i=0; i < mmin(s3.Size(), sd.Size()); i++)
                R1 += pot[i] * (s3.f[i] * sd.f[i] + Constant::AlphaSquared * s3.g[i] * sd.g[i]) * dR[i];

            if(SMS_ac)
            {   R1 = R1 - SMS_ac * SI.IsotopeShiftIntegral(s3, sd);
            }

            energy += R1 * SI.HamiltonianMatrixElement(sb, s3, *core) / (ValenceEnergies.find(sb.Kappa())->second - s3.Energy() + delta);
        }

        if(s3.RequiredPQN() >= 5 && s3.Kappa() == sd.Kappa())
        {
            double R1 = 0.;
            for(i=0; i < mmin(sb.Size(), s3.Size()); i++)
                R1 += pot[i] * (sb.f[i] * s3.f[i] + Constant::AlphaSquared * sb.g[i] * s3.g[i]) * dR[i];

            if(SMS_ac)
            {   R1 = R1 + SMS_ac * SI.IsotopeShiftIntegral(s3, sb);
            }

            energy += R1 * SI.HamiltonianMatrixElement(sd, s3, *core) / (ValenceEnergies.find(sd.Kappa())->second - s3.Energy() + delta);
        }

        it3.Next();
    }

    if(debug)
        *outstream << "  " << energy * Constant::HartreeEnergy_cm << std::endl;
    return energy;
}

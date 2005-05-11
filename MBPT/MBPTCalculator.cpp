#include "Include.h"
#include "MBPTCalculator.h"
#include "Universal/Constant.h"
#include "Universal/CoulombIntegrator.h"
#include "HartreeFock/StateIntegrator.h"

#define MAX_K 12

MBPTCalculator::MBPTCalculator(Lattice* lat, const Core* atom_core, const ExcitedStates* excited_states):
    lattice(lat), core(atom_core), excited(excited_states), BrillouinWignerPT(false)
{}

void MBPTCalculator::GetSecondOrderSigma(int kappa, SigmaPotential* sigma) const
{
    DiscreteState s(lattice, abs(kappa)+1, kappa);

    CalculateCorrelation1and3(s, s, sigma);
    CalculateCorrelation2(s, s, sigma);
    CalculateCorrelation4(s, s, sigma);
}

double MBPTCalculator::GetSecondOrderSigma(const State* s, SigmaPotential* sigma) const
{
    double energy = GetSecondOrderSigma(s, s, sigma);

    *outstream << "Total energy from direct summation = " << (s->Energy() + energy)*Constant::HartreeEnergy_cm << std::endl;

    if(sigma)
    {   double matrix_element = sigma->GetMatrixElement(s->f, s->f);
        *outstream << "Total energy via sigma potential   = " << (s->Energy() + matrix_element)*Constant::HartreeEnergy_cm << std::endl;
    }

    return (s->Energy() + energy);    
}

double MBPTCalculator::GetSecondOrderSigma(const State* s1, const State* s2, SigmaPotential* sigma) const
{
    if(sigma != NULL)
        sigma->Reset();

    if(s1->Kappa() != s2->Kappa())
        return 0.;

    double ten1 = CalculateCorrelation1and3(*s1, *s2, sigma);
    double ten2 = CalculateCorrelation2(*s1, *s2, sigma);
    double ten4 = CalculateCorrelation4(*s1, *s2, sigma);

    double energy = 0.;
    if(sigma != NULL)
    {   energy = sigma->GetMatrixElement(s1->f, s2->f);
        *outstream << "\tSecond order contribution: " << energy * Constant::HartreeEnergy_cm << std::endl;
    }

    return (ten1 + ten2 + ten4);
}

double MBPTCalculator::CalculateCorrelation1and3(const State& si, const State& sf, SigmaPotential* sigma) const
{
    bool calculate_matrix_element = (si.Size() && sf.Size());
    bool calculate_sigma_potential = (sigma != NULL);
    unsigned int min_size = mmin(si.Size(), sf.Size());

    const double NuclearInverseMass = core->GetNuclearInverseMass();

    std::vector<double> Pot24;
    CoulombIntegrator I(*lattice);
    StateIntegrator SI(*lattice);
    const double* dR = lattice->dR();

    *outstream << "Cor 1+3:  ";
    double spacing = 1./double(core->NumStates() * excited->NumStates());
    double count = 0;
    unsigned int i;

    // Value of the matrix elements <si | Sigma | sf>
    double energy1 = 0., energy3 = 0.;

    // Firstly, get the loop 24
    ConstStateIterator it2 = core->GetConstStateIterator();
    while(!it2.AtEnd())
    {   const State& s2 = *(it2.GetState());

        ConstStateIterator it4 = excited->GetConstStateIterator();
        while(!it4.AtEnd())
        {   const State& s4 = *(it4.GetState());
            unsigned int start_k = (s2.L() + s4.L())%2;

            count += spacing;
            if(count >= 0.02)
            {   *logstream << ".";
                count -= 0.02;
            }

            for(unsigned int k=start_k; k<=MAX_K; k+=2)
            {
                double coeff24 = Constant::Electron3j(s2.J(), s4.J(), k, 0.5, -0.5);
                if(coeff24)
                {
                    coeff24 = coeff24 * coeff24 * (2. * s2.J() + 1.) * (2. * s4.J() + 1.)
                                                / (2. * k + 1.);
                    coeff24 = coeff24 * it4.Weight();

                    std::vector<double> density(mmin(s2.Size(), s4.Size()));
                    for(i=0; i<density.size(); i++)
                    {
                        density[i] = s2.f[i] * s4.f[i] + Constant::AlphaSquared * s2.g[i] * s4.g[i];
                    }
                    density.resize(core->GetHFPotential().size());
                    I.FastCoulombIntegrate(density, Pot24, k);

                    double SMS_24 = 0.;
                    if(NuclearInverseMass && (k == 1))
                        SMS_24 = -SI.IsotopeShiftIntegral(s4, s2);

                    // Correlation 1 has excited state 3
                    ConstStateIterator it3_1 = excited->GetConstStateIterator();
                    while(!it3_1.AtEnd())
                    {   const State& s3 = *(it3_1.GetState());

                        if((si.L() + s3.L())%2 == start_k)
                        {
                            double coeff = Constant::Electron3j(si.J(), s3.J(), k, 0.5, -0.5);

                            if(coeff)
                            {
                                coeff = coeff * coeff * coeff24 * (2. * s3.J() + 1.);
                                coeff = coeff * it3_1.Weight();

                                if(BrillouinWignerPT)
                                    coeff = coeff/(ValenceEnergy1 + s2.Energy() - s3.Energy() - s4.Energy());
                                else
                                    coeff = coeff/(si.Energy() + s2.Energy() - s3.Energy() - s4.Energy());

                                // R1 = R_k (i2, 34)
                                // R2 = R_k (f2, 34)
                                double R1 = 0.;
                                double R2 = 0.;
                                std::vector<double> Y;

                                if(calculate_sigma_potential)
                                {   Y.resize(s3.Size());
                                    for(i=0; i < s3.Size(); i++)
                                        Y[i] = Pot24[i] * s3.f[i];
                                }

                                if(calculate_matrix_element)
                                {   for(i=0; i < mmin(s3.Size(), min_size); i++)
                                    {   R1 = R1 + Pot24[i] * (s3.f[i] * si.f[i] + Constant::AlphaSquared * s3.g[i] * si.g[i]) * dR[i];
                                        R2 = R2 + Pot24[i] * (s3.f[i] * sf.f[i] + Constant::AlphaSquared * s3.g[i] * sf.g[i]) * dR[i];
                                    }
                                    while(i < mmin(s3.Size(), si.Size()))
                                    {   R1 = R1 + Pot24[i] * (s3.f[i] * si.f[i] + Constant::AlphaSquared * s3.g[i] * si.g[i]) * dR[i];
                                        i++;
                                    }
                                    while(i < mmin(s3.Size(), sf.Size()))
                                    {   R2 = R2 + Pot24[i] * (s3.f[i] * sf.f[i] + Constant::AlphaSquared * s3.g[i] * sf.g[i]) * dR[i];
                                        i++;
                                    }
                                }

                                if(SMS_24)
                                {   double R1_sms;

                                    if(calculate_sigma_potential)
                                    {   std::vector<double> P(Y.size());
                                        R1_sms = SI.IsotopeShiftIntegral(si, s3, &P);
                                        for(i=0; i < Y.size(); i++)
                                            Y[i] = Y[i] - NuclearInverseMass * SMS_24 * P[i];
                                    }
                                    else
                                        R1_sms = SI.IsotopeShiftIntegral(si, s3);

                                    if(calculate_matrix_element)
                                    {   double R2_sms = SI.IsotopeShiftIntegral(sf, s3);
                                        R1 = R1 - R1_sms * SMS_24 * NuclearInverseMass;
                                        R2 = R2 - R2_sms * SMS_24 * NuclearInverseMass;
                                    }
                                }

                                if(calculate_sigma_potential)
                                    sigma->AddToSigma(Y, Y, coeff);

                                if(calculate_matrix_element)
                                    energy1 += R1 * R2 * coeff;
                            }
                        }
                        it3_1.Next();
                    }

                    // Correlation 3 has core state 3
                    ConstStateIterator it3_3 = core->GetConstStateIterator();
                    while(!it3_3.AtEnd())
                    {   const State& s3 = *(it3_3.GetState());

                        if((si.L() + s3.L())%2 == start_k)
                        {
                            double coeff = Constant::Electron3j(si.J(), s3.J(), k, 0.5, -0.5);

                            if(coeff)
                            {
                                coeff = coeff * coeff * coeff24 * (2. * s3.J() + 1.);

                                if(BrillouinWignerPT)
                                    coeff = coeff/(ValenceEnergy1 + s4.Energy() - s2.Energy() - s3.Energy());
                                else
                                    coeff = coeff/(sf.Energy() + s4.Energy() - s2.Energy() - s3.Energy());

                                // R1 = R_k (i2, 34)
                                // R2 = R_k (f2, 34)
                                double R1 = 0.;
                                double R2 = 0.;
                                std::vector<double> Y;

                                if(calculate_sigma_potential)
                                {   Y.resize(s3.Size());
                                    for(i=0; i < s3.Size(); i++)
                                        Y[i] = Pot24[i] * s3.f[i];
                                }

                                if(calculate_matrix_element)
                                {   for(i=0; i < mmin(s3.Size(), min_size); i++)
                                    {   R1 = R1 + Pot24[i] * (s3.f[i] * si.f[i] + Constant::AlphaSquared * s3.g[i] * si.g[i]) * dR[i];
                                        R2 = R2 + Pot24[i] * (s3.f[i] * sf.f[i] + Constant::AlphaSquared * s3.g[i] * sf.g[i]) * dR[i];
                                    }
                                    while(i < mmin(s3.Size(), si.Size()))
                                    {   R1 = R1 + Pot24[i] * (s3.f[i] * si.f[i] + Constant::AlphaSquared * s3.g[i] * si.g[i]) * dR[i];
                                        i++;
                                    }
                                    while(i < mmin(s3.Size(), sf.Size()))
                                    {   R2 = R2 + Pot24[i] * (s3.f[i] * sf.f[i] + Constant::AlphaSquared * s3.g[i] * sf.g[i]) * dR[i];
                                        i++;
                                    }
                                }

                                if(SMS_24)
                                {   double R1_sms;

                                    if(calculate_sigma_potential)
                                    {   std::vector<double> P(Y.size());
                                        R1_sms = -SI.IsotopeShiftIntegral(si, s3, &P);
                                        for(i=0; i < Y.size(); i++)
                                            // Plus sign is because P is opposite sign
                                            Y[i] = Y[i] + NuclearInverseMass * SMS_24 * P[i];
                                    }
                                    else
                                        R1_sms = -SI.IsotopeShiftIntegral(si, s3);

                                    if(calculate_matrix_element)
                                    {   double R2_sms = -SI.IsotopeShiftIntegral(sf, s3);
                                        R1 = R1 - R1_sms * SMS_24 * NuclearInverseMass;
                                        R2 = R2 - R2_sms * SMS_24 * NuclearInverseMass;
                                    }
                                }

                                if(calculate_sigma_potential)
                                    sigma->AddToSigma(Y, Y, coeff);

                                if(calculate_matrix_element)
                                    energy3 += R1 * R2 * coeff;
                            }
                        }
                        it3_3.Next();
                    }

                } // coeff24
            } // k
            it4.Next();
        }
        it2.Next();
    }

    if(calculate_matrix_element)
        *outstream << "  " << energy1 * Constant::HartreeEnergy_cm
                   << "  " << energy3 * Constant::HartreeEnergy_cm << std::endl;
    return energy1 + energy3;
}

double MBPTCalculator::CalculateCorrelation2(const State& si, const State& sf, SigmaPotential* sigma) const
{
    bool calculate_matrix_element = (si.Size() && sf.Size());
    bool calculate_sigma_potential = (sigma != NULL);

    const double NuclearInverseMass = core->GetNuclearInverseMass();

    std::vector<double> Pot24, Pot23;
    CoulombIntegrator I(*lattice);
    StateIntegrator SI(*lattice);
    const double* dR = lattice->dR();

    *outstream << "Cor 2:    ";
    double spacing = 1./double(core->NumStates() * excited->NumStates());
    double count = 0;
    unsigned int i;

    // Value of the matrix element <si | Sigma | sf>
    double energy = 0.;

    ConstStateIterator it2 = core->GetConstStateIterator();
    while(!it2.AtEnd())
    {   const State& s2 = *(it2.GetState());

        ConstStateIterator it4 = excited->GetConstStateIterator();
        while(!it4.AtEnd())
        {   const State& s4 = *(it4.GetState());
            unsigned int start_k1 = (s2.L() + s4.L())%2;

            count += spacing;
            if(count >= 0.02)
            {   *logstream << ".";
                count -= 0.02;
            }

            for(unsigned int k1=start_k1; k1<=MAX_K; k1+=2)
            {
                double coeff24 = Constant::Electron3j(s2.J(), s4.J(), k1, 0.5, -0.5);
                if(coeff24)
                {
                    coeff24 = coeff24 * (2. * s2.J() + 1.) * (2. * s4.J() + 1.);
                    coeff24 = coeff24 * it4.Weight();

                    std::vector<double> density(mmin(s2.Size(), s4.Size()));
                    for(i=0; i<density.size(); i++)
                    {
                        density[i] = s2.f[i] * s4.f[i];
                    }
                    density.resize(core->GetHFPotential().size());
                    I.FastCoulombIntegrate(density, Pot24, k1);

                    double SMS_24 = 0.;
                    if(NuclearInverseMass && (k1 == 1))
                        SMS_24 = -SI.IsotopeShiftIntegral(s4, s2);

                    ConstStateIterator it3 = excited->GetConstStateIterator();
                    while(!it3.AtEnd())
                    {   const State& s3 = *(it3.GetState());

                        double coeff13 = Constant::Electron3j(si.J(), s3.J(), k1, 0.5, -0.5);
                        unsigned int start_k2 = (s2.L() + s3.L())%2;
                        if(coeff13 && ((si.L() + s4.L())%2 == start_k2) && ((si.L() + s3.L())%2 == start_k1))
                        {
                            coeff13 = coeff13 * (2. * s3.J() + 1.) * it3.Weight();

                            if(BrillouinWignerPT)
                                coeff13 = coeff13/(ValenceEnergy1 + s2.Energy() - s3.Energy() - s4.Energy());
                            else
                                coeff13 = coeff13/(si.Energy() + s2.Energy() - s3.Energy() - s4.Energy());

                            // Sign
                            if((start_k1 + start_k2)%2)
                                coeff13 = -coeff13;

                            for(unsigned int k2=start_k2; k2<=MAX_K; k2+=2)
                            {
                                double coeff 
                                    = coeff13 * coeff24 * Constant::Electron3j(si.J(), s4.J(), k2, 0.5, -0.5)
                                    * Constant::Electron3j(s3.J(), s2.J(), k2, 0.5, -0.5)
                                    * Constant::Wigner6j(si.J(), s3.J(), k1, s2.J(), s4.J(), k2);

                                if(coeff)
                                {
                                    density.clear();
                                    density.resize(mmin(s2.Size(), s3.Size()));
                                    for(i=0; i<density.size(); i++)
                                    {
                                        density[i] = s2.f[i] * s3.f[i];
                                    }
                                    density.resize(core->GetHFPotential().size());
                                    I.FastCoulombIntegrate(density, Pot23, k2);

                                    double SMS_23 = 0.;
                                    if(NuclearInverseMass && (k2 == 1))
                                        SMS_23 = SI.IsotopeShiftIntegral(s3, s2);

                                    // R1 = R_k1 (i2, 34)
                                    double R1 = 0.;
                                    std::vector<double> Y1;

                                    if(calculate_sigma_potential)
                                    {   Y1.resize(s3.Size());
                                        for(i=0; i < s3.Size(); i++)
                                            Y1[i] = Pot24[i] * s3.f[i];
                                    }

                                    if(calculate_matrix_element)
                                    {   for(i=0; i < mmin(s3.Size(), si.Size()); i++)
                                        {   R1 = R1 + Pot24[i] * (s3.f[i] * si.f[i] + Constant::AlphaSquared * s3.g[i] * si.g[i]) * dR[i];
                                        }
                                    }
                               
                                    if(SMS_24)
                                    {   double R1_sms;

                                        if(calculate_sigma_potential)
                                        {   std::vector<double> P1(Y1.size());
                                            R1_sms = SI.IsotopeShiftIntegral(si, s3, &P1);
                                            for(i=0; i < Y1.size(); i++)
                                                Y1[i] = Y1[i] - NuclearInverseMass * SMS_24 * P1[i];
                                        }
                                        else
                                            R1_sms = SI.IsotopeShiftIntegral(si, s3);

                                        if(calculate_matrix_element)
                                            R1 = R1 - R1_sms * SMS_24 * NuclearInverseMass;
                                    }

                                    // R2 = R_k2 (34, 2f) = R_k2 (f2, 43)
                                    double R2 = 0.;
                                    std::vector<double> Y2;

                                    if(calculate_sigma_potential)
                                    {   Y2.resize(s4.Size());
                                        for(i=0; i < s4.Size(); i++)
                                            Y2[i] = Pot23[i] * s4.f[i];
                                    }

                                    if(calculate_matrix_element)
                                    {   for(i=0; i < mmin(s4.Size(), sf.Size()); i++)
                                        {   R2 = R2 + Pot23[i] * (s4.f[i] * sf.f[i] + Constant::AlphaSquared * s4.g[i] * sf.g[i]) * dR[i];
                                        }
                                    }

                                    if(SMS_23)
                                    {   double R2_sms;

                                        if(calculate_sigma_potential)
                                        {   std::vector<double> P2(Y2.size());
                                            R2_sms = -SI.IsotopeShiftIntegral(sf, s4, &P2);

                                            for(i=0; i<Y2.size(); i++)
                                                // P2 is wrong sign
                                                Y2[i] = Y2[i] + NuclearInverseMass * SMS_23 * P2[i];
                                        }
                                        else
                                            R2_sms = -SI.IsotopeShiftIntegral(sf, s4);

                                        if(calculate_matrix_element)
                                            R2 = R2 - R2_sms * SMS_23 * NuclearInverseMass;
                                    }

                                    if(calculate_sigma_potential)
                                        sigma->AddToSigma(Y1, Y2, coeff);

                                    if(calculate_matrix_element)
                                        energy += R1 * R2 * coeff;
                                }
                            } // k2
                        }
                        it3.Next();
                    } 
                } // coeff24
            } // k1
            it4.Next();
        }
        it2.Next();
    }

    if(calculate_matrix_element)
        *outstream << "  " << energy * Constant::HartreeEnergy_cm << std::endl;
    return energy;
}

double MBPTCalculator::CalculateCorrelation4(const State& si, const State& sf, SigmaPotential* sigma) const
{
    bool calculate_matrix_element = (si.Size() && sf.Size());
    bool calculate_sigma_potential = (sigma != NULL);

    const double NuclearInverseMass = core->GetNuclearInverseMass();

    std::vector<double> Pot24, Pot34;
    CoulombIntegrator I(*lattice);
    StateIntegrator SI(*lattice);
    const double* dR = lattice->dR();

    *outstream << "Cor 4:    ";
    double spacing = 1./double(core->NumStates() * excited->NumStates());
    double count = 0;
    unsigned int i;

    double energy = 0.;
    ConstStateIterator it4 = excited->GetConstStateIterator();
    while(!it4.AtEnd())
    {   const State& s4 = *(it4.GetState());

        ConstStateIterator it2 = core->GetConstStateIterator();
        while(!it2.AtEnd())
        {   const State& s2 = *(it2.GetState());
            unsigned int start_k1 = (s2.L() + s4.L())%2;

            count += spacing;
            if(count >= 0.02)
            {   *logstream << ".";
                count -= 0.02;
            }

            for(unsigned int k1=start_k1; k1<=MAX_K; k1+=2)
            {
                double coeff24 = Constant::Electron3j(s2.J(), s4.J(), k1, 0.5, -0.5);
                if(coeff24)
                {
                    coeff24 = coeff24 * (2. * s2.J() + 1.) * (2. * s4.J() + 1.);
                    coeff24 = coeff24 * it4.Weight();

                    std::vector<double> density(mmin(s2.Size(), s4.Size()));
                    for(i=0; i<density.size(); i++)
                    {
                        density[i] = s2.f[i] * s4.f[i];
                    }
                    density.resize(core->GetHFPotential().size());
                    I.FastCoulombIntegrate(density, Pot24, k1);

                    double SMS_24 = 0.;
                    if(NuclearInverseMass && (k1 == 1))
                        SMS_24 = -SI.IsotopeShiftIntegral(s4, s2);

                    ConstStateIterator it3 = core->GetConstStateIterator();
                    while(!it3.AtEnd())
                    {   const State& s3 = *(it3.GetState());

                        double coeff13 = Constant::Electron3j(si.J(), s3.J(), k1, 0.5, -0.5);
                        unsigned int start_k2 = (s3.L() + s4.L())%2;
                        if(coeff13 && ((si.L() + s2.L())%2 == start_k2) && ((si.L() + s3.L())%2 == start_k1))
                        {
                            coeff13 = coeff13 * (2. * s3.J() + 1.);

                            if(BrillouinWignerPT)
                                coeff13 = coeff13/(ValenceEnergy1 + s4.Energy() - s2.Energy() - s3.Energy());
                            else
                                coeff13 = coeff13/(sf.Energy() + s4.Energy() - s2.Energy() - s3.Energy());

                            // Sign
                            if((start_k1 + start_k2)%2)
                                coeff13 = -coeff13;

                            for(unsigned int k2= start_k2; k2<=MAX_K; k2+=2)
                            {
                                double coeff
                                    = coeff13 * coeff24 * Constant::Electron3j(si.J(), s2.J(), k2, 0.5, -0.5)
                                    * Constant::Electron3j(s3.J(), s4.J(), k2, 0.5, -0.5)
                                    * Constant::Wigner6j(si.J(), s3.J(), k1, s4.J(), s2.J(), k2);

                                if(coeff)
                                {
                                    density.clear();
                                    density.resize(mmin(s3.Size(), s4.Size()));
                                    for(i=0; i<density.size(); i++)
                                    {
                                        density[i] = s3.f[i] * s4.f[i];
                                    }
                                    density.resize(core->GetHFPotential().size());
                                    I.FastCoulombIntegrate(density, Pot34, k2);

                                    double SMS_34 = 0.;
                                    if(NuclearInverseMass && (k2 == 1))
                                        SMS_34 = SI.IsotopeShiftIntegral(s4, s3);

                                    // R1 = R_k1 (i4, 32)
                                    double R1 = 0.;
                                    std::vector<double> Y1;

                                    if(calculate_sigma_potential)
                                    {   Y1.resize(s3.Size());
                                        for(i=0; i < s3.Size(); i++)
                                            Y1[i] = Pot24[i] * s3.f[i];
                                    }

                                    if(calculate_matrix_element)
                                    {   for(i=0; i < mmin(s3.Size(), si.Size()); i++)
                                            R1 = R1 + Pot24[i] * (s3.f[i] * si.f[i] + Constant::AlphaSquared * s3.g[i] * si.g[i]) * dR[i];
                                    }

                                    if(SMS_24)
                                    {   double R1_sms;
                                        if(calculate_sigma_potential)
                                        {   std::vector<double> P1(Y1.size());
                                            R1_sms = -SI.IsotopeShiftIntegral(si, s3, &P1);
                                            for(i=0; i<Y1.size(); i++)
                                                Y1[i] = Y1[i] + NuclearInverseMass * SMS_24 * P1[i];
                                        }
                                        else
                                            R1_sms = -SI.IsotopeShiftIntegral(si, s3);

                                        if(calculate_matrix_element)
                                            R1 = R1 - R1_sms * SMS_24 * NuclearInverseMass;
                                    }

                                    // R2 = R_k2 (32, 4f) = R_k2 (f4, 23)
                                    double R2 = 0.;
                                    std::vector<double> Y2;

                                    if(calculate_sigma_potential)
                                    {   Y2.resize(s2.Size());
                                        for(i=0; i < s2.Size(); i++)
                                            Y2[i] = Pot34[i] * s2.f[i];
                                    }

                                    if(calculate_matrix_element)
                                        for(i=0; i < mmin(sf.Size(), s2.Size()); i++)
                                            R2 = R2 + Pot34[i] * (s2.f[i] * sf.f[i] + Constant::AlphaSquared * s2.g[i] * sf.g[i]) * dR[i];

                                    if(SMS_34)
                                    {   double R2_sms;

                                        if(calculate_sigma_potential)
                                        {   std::vector<double> P2(Y2.size());
                                            R2_sms = SI.IsotopeShiftIntegral(sf, s2, &P2);
                                            for(i=0; i < Y2.size(); i++)
                                                Y2[i] = Y2[i] - NuclearInverseMass * SMS_34 * P2[i];
                                        }
                                        else
                                            R2_sms = SI.IsotopeShiftIntegral(sf, s2);

                                        if(calculate_matrix_element)
                                            R2 = R2 - R2_sms * SMS_34 * NuclearInverseMass;
                                    }

                                    if(calculate_sigma_potential)
                                        sigma->AddToSigma(Y1, Y2, coeff);

                                    energy += R1 * R2 * coeff;
                                }
                            } // k2
                        }
                        it3.Next();
                    } 
                } // coeff24
            } // k1
            it2.Next();
        }
        it4.Next();
    }

    if(calculate_matrix_element)
        *outstream << "  " << energy * Constant::HartreeEnergy_cm << std::endl;
    return energy;
}

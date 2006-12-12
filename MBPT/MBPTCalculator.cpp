#include "Include.h"
#include "MBPTCalculator.h"
#include "Universal/Constant.h"
#include "Universal/CoulombIntegrator.h"
#include "HartreeFock/StateIntegrator.h"

#define MAX_K 12

MBPTCalculator::MBPTCalculator(Lattice* lat, const Core* atom_core, const ExcitedStates* excited_states):
    lattice(lat), core(atom_core), excited(excited_states), BrillouinWignerPT(false), delta(0.)
{
    MaxStateSize = core->GetConstHFPotential().size();
    SetValenceEnergies();
}

void MBPTCalculator::GetSecondOrderSigma(int kappa, SigmaPotential* sigma)
{
    MaxStateSize = core->GetConstHFPotential().size();

    DiscreteState s(lattice, abs(kappa)+1, kappa);

    CalculateCorrelation1and3(s, s, sigma);
    CalculateCorrelation2(s, s, sigma);
    CalculateCorrelation4(s, s, sigma);
}

double MBPTCalculator::GetSecondOrderSigma(const State* s, SigmaPotential* sigma)
{
    double energy = GetSecondOrderSigma(s, s, sigma);

    *outstream << "Total energy from direct summation = " << (s->Energy() + energy)*Constant::HartreeEnergy_cm << std::endl;

    if(sigma)
    {   double matrix_element = sigma->GetMatrixElement(s->f, s->f);
        *outstream << "Total energy via sigma potential   = " << (s->Energy() + matrix_element)*Constant::HartreeEnergy_cm << std::endl;
    }

    return (s->Energy() + energy);    
}

double MBPTCalculator::GetSecondOrderSigma(const State* s1, const State* s2, SigmaPotential* sigma)
{
    MaxStateSize = core->GetConstHFPotential().size();

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

double MBPTCalculator::GetSigmaSubtraction(const State* s1, const State* s2)
{
    MaxStateSize = core->GetConstHFPotential().size();

    if(s1->Kappa() != s2->Kappa())
        return 0.;

    double term1 = CalculateSubtraction1(*s1, *s2);
    double term2 = CalculateSubtraction2(*s1, *s2);
    double term3 = CalculateSubtraction3(*s1, *s2);
    
    return (term1 + term2 + term3);
}

double MBPTCalculator::GetTwoElectronDiagrams(const State* s1, const State* s2, const State* s3, const State* s4, unsigned int k)
{
    MaxStateSize = core->GetConstHFPotential().size();

    double term = 0;
    term += CalculateTwoElectron1(*s1, *s2, *s3, *s4, k);
    term += CalculateTwoElectron2(*s1, *s2, *s3, *s4, k);
    term += CalculateTwoElectron3(*s1, *s2, *s3, *s4, k);
    term += CalculateTwoElectron4(*s1, *s2, *s3, *s4, k);
    term += CalculateTwoElectron5(*s1, *s2, *s3, *s4, k);
    term += CalculateTwoElectron6(*s1, *s2, *s3, *s4, k);

    return term;
}

double MBPTCalculator::GetTwoElectronBoxDiagrams(const State* s1, const State* s2, const State* s3, const State* s4, unsigned int k)
{
    MaxStateSize = core->GetConstHFPotential().size();

    double term = 0;
    term += CalculateTwoElectron4(*s1, *s2, *s3, *s4, k);
    term += CalculateTwoElectron5(*s1, *s2, *s3, *s4, k);
    term += CalculateTwoElectron6(*s1, *s2, *s3, *s4, k);

    return term;
}

double MBPTCalculator::GetTwoElectronSubtraction(const State* s1, const State* s2, const State* s3, const State* s4, unsigned int k)
{
    MaxStateSize = core->GetConstHFPotential().size();

    return CalculateTwoElectronSub(*s1, *s2, *s3, *s4, k);
}

double MBPTCalculator::CalculateCorrelation1and3(const State& si, const State& sf, SigmaPotential* sigma) const
{
    bool calculate_matrix_element = (si.Size() && sf.Size());
    bool calculate_sigma_potential = (sigma != NULL);
    unsigned int min_size = mmin(si.Size(), sf.Size());
    const bool debug = DebugOptions.LogMBPT();
    const double NuclearInverseMass = core->GetNuclearInverseMass();

    unsigned int Ji = (unsigned int)(si.J() * 2.);

    std::vector<double> density(MaxStateSize);
    std::vector<double> Pot24(MaxStateSize);
    std::vector<double> Y(MaxStateSize); unsigned int Y_size;
    CoulombIntegrator I(lattice);
    StateIntegrator SI(lattice);
    const double* dR = lattice->dR();

    if(debug)
        *outstream << "Cor 1+3:  ";
    double spacing = 1./double(core->NumStates() * excited->NumStates());
    double count = 0.;
    unsigned int i;

    // Value of the matrix elements <si | Sigma | sf>
    double energy1 = 0., energy3 = 0.;
    const double ValenceEnergy = ValenceEnergies.find(si.Kappa())->second;

    // Firstly, get the loop 24
    ConstStateIterator it2 = core->GetConstStateIterator();
    while(!it2.AtEnd())
    {   const State& s2 = *(it2.GetState());
        unsigned int J2 = (unsigned int)(s2.J() * 2.);

        ConstStateIterator it4 = excited->GetConstStateIterator();
        while(!it4.AtEnd())
        {   const State& s4 = *(it4.GetState());
            unsigned int J4 = (unsigned int)(s4.J() * 2.);

            if(debug)
            {   count += spacing;
                if(count >= 0.02)
                {   *logstream << ".";
                    count -= 0.02;
                }
            }

            unsigned int start_k = (s2.L() + s4.L())%2;
            for(unsigned int k=start_k; k<=MAX_K; k+=2)
            {
                double coeff24 = Constant::Electron3j(J2, J4, k, 1, -1);
                if(coeff24)
                {
                    coeff24 = coeff24 * coeff24 * (J2 + 1) * (J4 + 1)
                                                / (2. * k + 1.);
                    coeff24 = coeff24 * it4.Weight();

                    for(i=0; i<mmin(s2.Size(), s4.Size()); i++)
                    {
                        density[i] = s2.f[i] * s4.f[i] + Constant::AlphaSquared * s2.g[i] * s4.g[i];
                    }
                    I.FastCoulombIntegrate(density, Pot24, k, mmin(s2.Size(), s4.Size()));
    
                    double SMS_24 = 0.;
                    if(NuclearInverseMass && (k == 1))
                        SMS_24 = -SI.IsotopeShiftIntegral(s4, s2);

                    // Correlation 1 has excited state 3
                    ConstStateIterator it3_1 = excited->GetConstStateIterator();
                    while(!it3_1.AtEnd())
                    {   const State& s3 = *(it3_1.GetState());
                        unsigned int J3 = (unsigned int)(s3.J() * 2.);

                        if((si.L() + s3.L())%2 == start_k)
                        {
                            double coeff = Constant::Electron3j(Ji, J3, k, 1, -1);

                            if(coeff)
                            {
                                coeff = coeff * coeff * coeff24 * (J3 + 1);

                                if(BrillouinWignerPT)
                                    coeff = coeff/(ValenceEnergy + s2.Energy() - s3.Energy() - s4.Energy() + delta);
                                else
                                    coeff = coeff/(si.Energy() + s2.Energy() - s3.Energy() - s4.Energy());

                                // R1 = R_k (i2, 34)
                                // R2 = R_k (f2, 34)
                                double R1 = 0.;
                                double R2 = 0.;

                                if(calculate_sigma_potential)
                                {   Y_size = s3.Size();
                                    for(i=0; i < Y_size; i++)
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
                                    {   std::vector<double> P(Y_size);
                                        R1_sms = SI.IsotopeShiftIntegral(si, s3, &P);
                                        for(i=0; i < Y_size; i++)
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
                                    sigma->AddToSigma(Y, Y, coeff, Y_size, Y_size);

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
                        unsigned int J3 = (unsigned int)(s3.J() * 2.);

                        if((si.L() + s3.L())%2 == start_k)
                        {
                            double coeff = Constant::Electron3j(Ji, J3, k, 1, -1);

                            if(coeff)
                            {
                                coeff = coeff * coeff * coeff24 * (J3 + 1);

                                if(BrillouinWignerPT)
                                    coeff = coeff/(ValenceEnergy + s4.Energy() - s2.Energy() - s3.Energy() - delta);
                                else
                                    coeff = coeff/(sf.Energy() + s4.Energy() - s2.Energy() - s3.Energy());

                                // R1 = R_k (i2, 34)
                                // R2 = R_k (f2, 34)
                                double R1 = 0.;
                                double R2 = 0.;

                                if(calculate_sigma_potential)
                                {   Y_size = s3.Size();
                                    for(i=0; i < Y_size; i++)
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
                                    {   std::vector<double> P(Y_size);
                                        R1_sms = -SI.IsotopeShiftIntegral(si, s3, &P);
                                        for(i=0; i < Y_size; i++)
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
                                    sigma->AddToSigma(Y, Y, coeff, Y_size, Y_size);

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

    if(calculate_matrix_element && debug)
        *outstream << "  " << energy1 * Constant::HartreeEnergy_cm
                   << "  " << energy3 * Constant::HartreeEnergy_cm << std::endl;
    return energy1 + energy3;
}

double MBPTCalculator::CalculateCorrelation2(const State& si, const State& sf, SigmaPotential* sigma) const
{
    bool calculate_matrix_element = (si.Size() && sf.Size());
    bool calculate_sigma_potential = (sigma != NULL);
    const bool debug = DebugOptions.LogMBPT();
    const double NuclearInverseMass = core->GetNuclearInverseMass();

    unsigned int Ji = (unsigned int)(si.J() * 2.);

    std::vector<double> density(MaxStateSize);
    std::vector<double> Pot24(MaxStateSize);
    std::vector<double> Pot23(MaxStateSize);
    std::vector<double> Y1(MaxStateSize); unsigned int Y1_size;
    std::vector<double> Y2(MaxStateSize); unsigned int Y2_size;
    CoulombIntegrator I(lattice);
    StateIntegrator SI(lattice);
    const double* dR = lattice->dR();

    if(debug)
        *outstream << "Cor 2:    ";
    double spacing = 1./double(core->NumStates() * excited->NumStates());
    double count = 0.;
    unsigned int i;

    // Value of the matrix element <si | Sigma | sf>
    double energy = 0.;
    const double ValenceEnergy = ValenceEnergies.find(si.Kappa())->second;

    ConstStateIterator it2 = core->GetConstStateIterator();
    while(!it2.AtEnd())
    {   const State& s2 = *(it2.GetState());
        unsigned int J2 = (unsigned int)(s2.J() * 2.);

        ConstStateIterator it4 = excited->GetConstStateIterator();
        while(!it4.AtEnd())
        {   const State& s4 = *(it4.GetState());
            unsigned int J4 = (unsigned int)(s4.J() * 2.);

            if(debug)
            {   count += spacing;
                if(count >= 0.02)
                {   *logstream << ".";
                    count -= 0.02;
                }
            }

            unsigned int start_k1 = (s2.L() + s4.L())%2;
            for(unsigned int k1=start_k1; k1<=MAX_K; k1+=2)
            {
                double coeff24 = Constant::Electron3j(J2, J4, k1, 1, -1);
                if(coeff24)
                {
                    coeff24 = coeff24 * (J2 + 1) * (J4 + 1);

                    for(i=0; i<mmin(s2.Size(), s4.Size()); i++)
                    {
                        density[i] = s2.f[i] * s4.f[i] + Constant::AlphaSquared * s2.g[i] * s4.g[i];
                    }
                    I.FastCoulombIntegrate(density, Pot24, k1, mmin(s2.Size(), s4.Size()));

                    double SMS_24 = 0.;
                    if(NuclearInverseMass && (k1 == 1))
                        SMS_24 = -SI.IsotopeShiftIntegral(s4, s2);

                    ConstStateIterator it3 = excited->GetConstStateIterator();
                    while(!it3.AtEnd())
                    {   const State& s3 = *(it3.GetState());
                        unsigned int J3 = (unsigned int)(s3.J() * 2.);

                        double coeff13 = Constant::Electron3j(Ji, J3, k1, 1, -1);
                        unsigned int start_k2 = (s2.L() + s3.L())%2;
                        if(coeff13 && ((si.L() + s4.L())%2 == start_k2) && ((si.L() + s3.L())%2 == start_k1))
                        {
                            coeff13 = coeff13 * (J3 + 1);

                            if(BrillouinWignerPT)
                                coeff13 = coeff13/(ValenceEnergy + s2.Energy() - s3.Energy() - s4.Energy() + delta);
                            else
                                coeff13 = coeff13/(si.Energy() + s2.Energy() - s3.Energy() - s4.Energy());

                            // Sign
                            if((start_k1 + start_k2)%2)
                                coeff13 = -coeff13;

                            for(unsigned int k2=start_k2; k2<=MAX_K; k2+=2)
                            {
                                double coeff 
                                    = coeff13 * coeff24 * Constant::Electron3j(Ji, J4, k2, 1, -1)
                                    * Constant::Electron3j(J3, J2, k2, 1, -1)
                                    * Constant::Wigner6j(si.J(), s3.J(), k1, s2.J(), s4.J(), k2);

                                if(coeff)
                                {
                                    for(i=0; i<mmin(s2.Size(), s3.Size()); i++)
                                    {
                                        density[i] = s2.f[i] * s3.f[i] + Constant::AlphaSquared * s2.g[i] * s3.g[i];
                                    }
                                    I.FastCoulombIntegrate(density, Pot23, k2, mmin(s2.Size(), s3.Size()));

                                    double SMS_23 = 0.;
                                    if(NuclearInverseMass && (k2 == 1))
                                        SMS_23 = SI.IsotopeShiftIntegral(s3, s2);

                                    // R1 = R_k1 (i2, 34)
                                    double R1 = 0.;

                                    if(calculate_sigma_potential)
                                    {   Y1_size = s3.Size();
                                        for(i=0; i < Y1_size; i++)
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
                                        {   std::vector<double> P1(Y1_size);
                                            R1_sms = SI.IsotopeShiftIntegral(si, s3, &P1);
                                            for(i=0; i < Y1_size; i++)
                                                Y1[i] = Y1[i] - NuclearInverseMass * SMS_24 * P1[i];
                                        }
                                        else
                                            R1_sms = SI.IsotopeShiftIntegral(si, s3);

                                        if(calculate_matrix_element)
                                            R1 = R1 - R1_sms * SMS_24 * NuclearInverseMass;
                                    }

                                    // R2 = R_k2 (34, 2f) = R_k2 (f2, 43)
                                    double R2 = 0.;

                                    if(calculate_sigma_potential)
                                    {   Y2_size = s4.Size();
                                        for(i=0; i < Y2_size; i++)
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
                                        {   std::vector<double> P2(Y2_size);
                                            R2_sms = -SI.IsotopeShiftIntegral(sf, s4, &P2);

                                            for(i=0; i<Y2_size; i++)
                                                // P2 is wrong sign
                                                Y2[i] = Y2[i] + NuclearInverseMass * SMS_23 * P2[i];
                                        }
                                        else
                                            R2_sms = -SI.IsotopeShiftIntegral(sf, s4);

                                        if(calculate_matrix_element)
                                            R2 = R2 - R2_sms * SMS_23 * NuclearInverseMass;
                                    }

                                    if(calculate_sigma_potential)
                                        sigma->AddToSigma(Y1, Y2, coeff, Y1_size, Y2_size);

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

    if(calculate_matrix_element && debug)
        *outstream << "  " << energy * Constant::HartreeEnergy_cm << std::endl;
    return energy;
}

double MBPTCalculator::CalculateCorrelation4(const State& si, const State& sf, SigmaPotential* sigma) const
{
    bool calculate_matrix_element = (si.Size() && sf.Size());
    bool calculate_sigma_potential = (sigma != NULL);
    const bool debug = DebugOptions.LogMBPT();
    const double NuclearInverseMass = core->GetNuclearInverseMass();

    unsigned int Ji = (unsigned int)(si.J() * 2.);

    std::vector<double> density(MaxStateSize);
    std::vector<double> Pot24(MaxStateSize);
    std::vector<double> Pot34(MaxStateSize);
    std::vector<double> Y1(MaxStateSize); unsigned int Y1_size;
    std::vector<double> Y2(MaxStateSize); unsigned int Y2_size;
    CoulombIntegrator I(lattice);
    StateIntegrator SI(lattice);
    const double* dR = lattice->dR();

    if(debug)
        *outstream << "Cor 4:    ";
    double spacing = 1./double(core->NumStates() * excited->NumStates());
    double count = 0.;
    unsigned int i;

    double energy = 0.;
    const double ValenceEnergy = ValenceEnergies.find(si.Kappa())->second;

    ConstStateIterator it4 = excited->GetConstStateIterator();
    while(!it4.AtEnd())
    {   const State& s4 = *(it4.GetState());
        unsigned int J4 = (unsigned int)(s4.J() * 2.);

        ConstStateIterator it2 = core->GetConstStateIterator();
        while(!it2.AtEnd())
        {   const State& s2 = *(it2.GetState());
            unsigned int J2 = (unsigned int)(s2.J() * 2.);

            if(debug)
            {   count += spacing;
                if(count >= 0.02)
                {   *logstream << ".";
                    count -= 0.02;
                }
            }

            unsigned int start_k1 = (s2.L() + s4.L())%2;
            for(unsigned int k1=start_k1; k1<=MAX_K; k1+=2)
            {
                double coeff24 = Constant::Electron3j(J2, J4, k1, 1, -1);
                if(coeff24)
                {
                    coeff24 = coeff24 * (J2 + 1) * (J4 + 1);

                    for(i=0; i<mmin(s2.Size(), s4.Size()); i++)
                    {
                        density[i] = s2.f[i] * s4.f[i] + Constant::AlphaSquared * s2.g[i] * s4.g[i];
                    }
                    I.FastCoulombIntegrate(density, Pot24, k1, mmin(s2.Size(), s4.Size()));

                    double SMS_24 = 0.;
                    if(NuclearInverseMass && (k1 == 1))
                        SMS_24 = -SI.IsotopeShiftIntegral(s4, s2);

                    ConstStateIterator it3 = core->GetConstStateIterator();
                    while(!it3.AtEnd())
                    {   const State& s3 = *(it3.GetState());
                        unsigned int J3 = (unsigned int)(s3.J() * 2.);

                        double coeff13 = Constant::Electron3j(Ji, J3, k1, 1, -1);
                        unsigned int start_k2 = (s3.L() + s4.L())%2;
                        if(coeff13 && ((si.L() + s2.L())%2 == start_k2) && ((si.L() + s3.L())%2 == start_k1))
                        {
                            coeff13 = coeff13 * (J3 + 1);

                            if(BrillouinWignerPT)
                                coeff13 = coeff13/(ValenceEnergy + s4.Energy() - s2.Energy() - s3.Energy() - delta);
                            else
                                coeff13 = coeff13/(sf.Energy() + s4.Energy() - s2.Energy() - s3.Energy());

                            // Sign
                            if((start_k1 + start_k2)%2)
                                coeff13 = -coeff13;

                            for(unsigned int k2= start_k2; k2<=MAX_K; k2+=2)
                            {
                                double coeff
                                    = coeff13 * coeff24 * Constant::Electron3j(Ji, J2, k2, 1, -1)
                                    * Constant::Electron3j(J3, J4, k2, 1, -1)
                                    * Constant::Wigner6j(si.J(), s3.J(), k1, s4.J(), s2.J(), k2);

                                if(coeff)
                                {
                                    for(i=0; i<mmin(s3.Size(), s4.Size()); i++)
                                    {
                                        density[i] = s3.f[i] * s4.f[i] + Constant::AlphaSquared * s3.g[i] * s4.g[i];
                                    }
                                    I.FastCoulombIntegrate(density, Pot34, k2, mmin(s3.Size(), s4.Size()));

                                    double SMS_34 = 0.;
                                    if(NuclearInverseMass && (k2 == 1))
                                        SMS_34 = SI.IsotopeShiftIntegral(s4, s3);

                                    // R1 = R_k1 (i4, 32)
                                    double R1 = 0.;

                                    if(calculate_sigma_potential)
                                    {   Y1_size = s3.Size();
                                        for(i=0; i < Y1_size; i++)
                                            Y1[i] = Pot24[i] * s3.f[i];
                                    }

                                    if(calculate_matrix_element)
                                    {   for(i=0; i < mmin(s3.Size(), si.Size()); i++)
                                            R1 = R1 + Pot24[i] * (s3.f[i] * si.f[i] + Constant::AlphaSquared * s3.g[i] * si.g[i]) * dR[i];
                                    }

                                    if(SMS_24)
                                    {   double R1_sms;
                                        if(calculate_sigma_potential)
                                        {   std::vector<double> P1(Y1_size);
                                            R1_sms = -SI.IsotopeShiftIntegral(si, s3, &P1);
                                            for(i=0; i<Y1_size; i++)
                                                Y1[i] = Y1[i] + NuclearInverseMass * SMS_24 * P1[i];
                                        }
                                        else
                                            R1_sms = -SI.IsotopeShiftIntegral(si, s3);

                                        if(calculate_matrix_element)
                                            R1 = R1 - R1_sms * SMS_24 * NuclearInverseMass;
                                    }

                                    // R2 = R_k2 (32, 4f) = R_k2 (f4, 23)
                                    double R2 = 0.;

                                    if(calculate_sigma_potential)
                                    {   Y2_size = s2.Size();
                                        for(i=0; i < Y2_size; i++)
                                            Y2[i] = Pot34[i] * s2.f[i];
                                    }

                                    if(calculate_matrix_element)
                                        for(i=0; i < mmin(sf.Size(), s2.Size()); i++)
                                            R2 = R2 + Pot34[i] * (s2.f[i] * sf.f[i] + Constant::AlphaSquared * s2.g[i] * sf.g[i]) * dR[i];

                                    if(SMS_34)
                                    {   double R2_sms;

                                        if(calculate_sigma_potential)
                                        {   std::vector<double> P2(Y2_size);
                                            R2_sms = SI.IsotopeShiftIntegral(sf, s2, &P2);
                                            for(i=0; i < Y2_size; i++)
                                                Y2[i] = Y2[i] - NuclearInverseMass * SMS_34 * P2[i];
                                        }
                                        else
                                            R2_sms = SI.IsotopeShiftIntegral(sf, s2);

                                        if(calculate_matrix_element)
                                            R2 = R2 - R2_sms * SMS_34 * NuclearInverseMass;
                                    }

                                    if(calculate_sigma_potential)
                                        sigma->AddToSigma(Y1, Y2, coeff, Y1_size, Y2_size);

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

    if(calculate_matrix_element && debug)
        *outstream << "  " << energy * Constant::HartreeEnergy_cm << std::endl;
    return energy;
}

double MBPTCalculator::CalculateSubtraction1(const State& si, const State& sf, SigmaPotential* sigma) const
{
    bool calculate_matrix_element = (si.Size() && sf.Size());
    bool calculate_sigma_potential = (sigma != NULL);
    const bool debug = DebugOptions.LogMBPT();

    unsigned int Ji = (unsigned int)(si.J() * 2.);

    std::vector<double> density(MaxStateSize);
    std::vector<double> Pot24(MaxStateSize);
    CoulombIntegrator I(lattice);
    StateIntegrator SI(lattice);
    const double* dR = lattice->dR();

    if(debug)
        *outstream << "Sub 1:    ";
    double spacing = 1./double(core->NumStates() * excited->NumStates());
    double count = 0.;
    unsigned int i;

    double energy = 0.;
    ConstStateIterator it2 = core->GetConstStateIterator();
    while(!it2.AtEnd())
    {   const State& s2 = *(it2.GetState());
        unsigned int J2 = (unsigned int)(s2.J() * 2.);

        ConstStateIterator it4 = excited->GetConstStateIterator();
        while(!it4.AtEnd())
        {   const State& s4 = *(it4.GetState());

            if(debug)
            {   count += spacing;
                if(count >= 0.02)
                {   *logstream << ".";
                    count -= 0.02;
                }
            }

            if(s2.Kappa() == s4.Kappa())
            {
                double coeff = SI.HamiltonianMatrixElement(s4, s2, *core);
                if(coeff)
                {
                    coeff = coeff * (J2 + 1);

                    for(i=0; i<mmin(s2.Size(), s4.Size()); i++)
                    {
                        density[i] = s2.f[i] * s4.f[i] + Constant::AlphaSquared * s2.g[i] * s4.g[i];
                    }
                    I.FastCoulombIntegrate(density, Pot24, 0, mmin(s2.Size(), s4.Size()));

                    coeff = coeff/(s2.Energy() - s4.Energy() + delta);

                    // R1 = R_0 (i2, f4)
                    double R1 = 0.;

                    if(calculate_matrix_element)
                    {   for(i=0; i < mmin(si.Size(), sf.Size()); i++)
                            R1 = R1 + Pot24[i] * (si.f[i] * sf.f[i] + Constant::AlphaSquared * si.g[i] * sf.g[i]) * dR[i];
                    }

                    if(calculate_sigma_potential)
                        sigma->AddDiagonal(Pot24, coeff);

                    energy += 2. * R1 * coeff;
                }
            }
            it4.Next();
        }
        it2.Next();
    }

    if(calculate_matrix_element && debug)
        *outstream << "  " << energy * Constant::HartreeEnergy_cm << std::endl;
    return energy;
}

double MBPTCalculator::CalculateSubtraction2(const State& si, const State& sf, SigmaPotential* sigma) const
{
    bool calculate_matrix_element = (si.Size() && sf.Size());
    bool calculate_sigma_potential = (sigma != NULL);
    const bool debug = DebugOptions.LogMBPT();
    const double NuclearInverseMass = core->GetNuclearInverseMass();

    unsigned int Ji = (unsigned int)(si.J() * 2.);

    std::vector<double> density(MaxStateSize);
    std::vector<double> Pot2f(MaxStateSize);
    std::vector<double> Pot4f(MaxStateSize);
    CoulombIntegrator I(lattice);
    StateIntegrator SI(lattice);
    const double* dR = lattice->dR();

    if(debug)
        *outstream << "Sub 2:    ";
    double spacing = 1./double(core->NumStates() * excited->NumStates());
    double count = 0.;
    unsigned int i;

    double energy = 0.;
    ConstStateIterator it2 = core->GetConstStateIterator();
    while(!it2.AtEnd())
    {   const State& s2 = *(it2.GetState());
        unsigned int J2 = (unsigned int)(s2.J() * 2.);

        ConstStateIterator it4 = excited->GetConstStateIterator();
        while(!it4.AtEnd())
        {   const State& s4 = *(it4.GetState());

            count += spacing;
            if(count >= 0.02)
            {   *logstream << ".";
                count -= 0.02;
            }

            if(s2.Kappa() == s4.Kappa())
            {
                double coeff24 = SI.HamiltonianMatrixElement(s4, s2, *core);
                if(coeff24)
                {
                    coeff24 = coeff24 * (J2 + 1);
                    coeff24 = coeff24/(s2.Energy() - s4.Energy() + delta);
            
                    unsigned int start_k = (si.L() + s2.L())%2;

                    for(unsigned int k=start_k; k<=MAX_K; k+=2)
                    {
                        double coeff = Constant::Electron3j(Ji, J2, k, 1, -1);
                        coeff = - coeff * coeff * coeff24;

                        if(coeff)
                        {
                            for(i=0; i<mmin(s2.Size(), sf.Size()); i++)
                            {
                                density[i] = s2.f[i] * sf.f[i] + Constant::AlphaSquared * s2.g[i] * sf.g[i];
                            }
                            I.FastCoulombIntegrate(density, Pot2f, k, mmin(s2.Size(), sf.Size()));

                            for(i=0; i<mmin(s4.Size(), sf.Size()); i++)
                            {
                                density[i] = s4.f[i] * sf.f[i] + Constant::AlphaSquared * s4.g[i] * sf.g[i];
                            }
                            I.FastCoulombIntegrate(density, Pot4f, k, mmin(s4.Size(), sf.Size()));
                            
                            // R1 = R_k (i4, 2f), R2 = R_k (i2, 4f)
                            double R1 = 0., R2 = 0.;        
                            if(calculate_matrix_element)
                            {   unsigned int limit = mmin(si.Size(), sf.Size());
                                for(i=0; i < mmin(limit, mmin(s2.Size(), s4.Size())); i++)
                                {
                                    R1 = R1 + Pot4f[i] * (si.f[i] * s2.f[i] + Constant::AlphaSquared * si.g[i] * s2.g[i]) * dR[i];
                                    R2 = R2 + Pot2f[i] * (si.f[i] * s4.f[i] + Constant::AlphaSquared * si.g[i] * s4.g[i]) * dR[i];
                                }
                                while(i < mmin(limit, s2.Size()))
                                {   R1 = R1 + Pot4f[i] * (si.f[i] * s2.f[i] + Constant::AlphaSquared * si.g[i] * s2.g[i]) * dR[i];
                                    i++;
                                }
                                while(i < mmin(limit, s4.Size()))
                                {   R2 = R2 + Pot2f[i] * (si.f[i] * s4.f[i] + Constant::AlphaSquared * si.g[i] * s4.g[i]) * dR[i];
                                    i++;
                                }
                            }

                            if(NuclearInverseMass && (k == 1))
                            {
                                R1 = R1 - NuclearInverseMass * SI.IsotopeShiftIntegral(s4, sf) * SI.IsotopeShiftIntegral(si, s2);
                                R2 = R2 - NuclearInverseMass * SI.IsotopeShiftIntegral(sf, s2) * SI.IsotopeShiftIntegral(s4, si);
                            }

                            if(calculate_matrix_element)
                                energy += (R1 + R2) * coeff;
                        }
                    }
                }
            }
            it4.Next();
        }
        it2.Next();
    }

    if(calculate_matrix_element && debug)
        *outstream << "  " << energy * Constant::HartreeEnergy_cm << std::endl;
    return energy;
}

double MBPTCalculator::CalculateSubtraction3(const State& si, const State& sf, SigmaPotential* sigma) const
{
    const bool debug = DebugOptions.LogMBPT();
    StateIntegrator SI(lattice);    

    if(debug)
        *outstream << "Sub 3:    ";

    double energy = 0.;
    const double ValenceEnergy = ValenceEnergies.find(si.Kappa())->second;

    ConstStateIterator it2 = core->GetConstStateIterator();
    while(!it2.AtEnd())
    {   const State& s2 = *(it2.GetState());

        if(s2.Kappa() == si.Kappa())
        {
            double term = SI.HamiltonianMatrixElement(si, s2, *core) * SI.HamiltonianMatrixElement(sf, s2, *core);
            if(BrillouinWignerPT)
                term = term/(s2.Energy() - ValenceEnergy + delta);
            else
                term = term/(s2.Energy() - si.Energy());
                
            energy = energy - term;
        }
        it2.Next();
    }

    if(debug)
        *outstream << "  " << energy * Constant::HartreeEnergy_cm << std::endl;
    return energy;
}

double MBPTCalculator::CalculateTwoElectron1(const State& sa, const State& sb, const State& sc, const State& sd, unsigned int k) const
{
    const bool debug = DebugOptions.LogMBPT();

    std::vector<double> Pot24(MaxStateSize);
    std::vector<double> density(MaxStateSize);

    CoulombIntegrator I(lattice);
    const double* dR = lattice->dR();

    if(debug)
        *outstream << "TwoE 1:   ";
    double spacing = 1./double(core->NumStates() * excited->NumStates());
    double count = 0.;
    unsigned int i;

    double energy = 0.;
    ConstStateIterator it2 = core->GetConstStateIterator();
    while(!it2.AtEnd())
    {   const State& s2 = *(it2.GetState());
        unsigned int J2 = (unsigned int)(s2.J() * 2.);

        ConstStateIterator it4 = excited->GetConstStateIterator();
        while(!it4.AtEnd())
        {   const State& s4 = *(it4.GetState());
            unsigned int J4 = (unsigned int)(s4.J() * 2.);

            if(debug)
            {   count += spacing;
                if(count >= 0.02)
                {   *logstream << ".";
                    count -= 0.02;
                }
            }

            double coeff = Constant::Electron3j(J2, J4, k, 1, -1);
            if((s2.L() + s4.L() + k)%2)
                coeff = 0.;

            if(coeff)
            {
                coeff = coeff * coeff * (J2 + 1) * (J4 + 1)
                                        / (2. * k + 1.);
                coeff = coeff/(s2.Energy() - s4.Energy() + delta);

                for(i=0; i<mmin(s2.Size(), s4.Size()); i++)
                {
                    density[i] = s2.f[i] * s4.f[i] + Constant::AlphaSquared * s2.g[i] * s4.g[i];
                }
                I.FastCoulombIntegrate(density, Pot24, k, mmin(s2.Size(), s4.Size()));

                // R1 = R_k (a2, c4)
                // R2 = R_k (b2, d4)
                // There is no first order SMS for this diagram.
                double R1 = 0.;
                double R2 = 0.;

                for(i=0; i < mmin(sa.Size(), sc.Size()); i++)
                    R1 = R1 + Pot24[i] * (sa.f[i] * sc.f[i] + Constant::AlphaSquared * sa.g[i] * sc.g[i]) * dR[i];
                
                for(i=0; i < mmin(sb.Size(), sd.Size()); i++)
                    R2 = R2 + Pot24[i] * (sb.f[i] * sd.f[i] + Constant::AlphaSquared * sb.g[i] * sd.g[i]) * dR[i];

                energy += 2. * R1 * R2 * coeff;
            }
            it4.Next();
        }
        it2.Next();
    }

    if(debug)
        *outstream << "  " << energy * Constant::HartreeEnergy_cm << std::endl;
    return energy;
}

double MBPTCalculator::CalculateTwoElectron2(const State& sa, const State& sb, const State& sc, const State& sd, unsigned int k) const
{
    const bool debug = DebugOptions.LogMBPT();
    const double NuclearInverseMass = core->GetNuclearInverseMass();

    unsigned int Ja = (unsigned int)(sa.J() * 2.);
    unsigned int Jb = (unsigned int)(sb.J() * 2.);
    unsigned int Jc = (unsigned int)(sc.J() * 2.);
    unsigned int Jd = (unsigned int)(sd.J() * 2.);

    std::vector<double> density(MaxStateSize);
    std::vector<double> Pot24(MaxStateSize);
    std::vector<double> Pota2(MaxStateSize);
    std::vector<double> Potb2(MaxStateSize);
    CoulombIntegrator I(lattice);
    StateIntegrator SI(lattice);
    const double* dR = lattice->dR();

    if(debug)
        *outstream << "TwoE 2:   ";
    double spacing = 1./double(core->NumStates() * excited->NumStates());
    double count = 0.;
    unsigned int i;

    double energy = 0.;
    const double coeff_ac = Constant::Electron3j(Ja, Jc, k, 1, -1);
    const double coeff_bd = Constant::Electron3j(Jb, Jd, k, 1, -1);
    if(!coeff_ac || !coeff_bd)
        return energy;

    ConstStateIterator it2 = core->GetConstStateIterator();
    while(!it2.AtEnd())
    {   const State& s2 = *(it2.GetState());
        unsigned int J2 = (unsigned int)(s2.J() * 2.);

        ConstStateIterator it4 = excited->GetConstStateIterator();
        while(!it4.AtEnd())
        {   const State& s4 = *(it4.GetState());
            unsigned int J4 = (unsigned int)(s4.J() * 2.);

            if(debug)
            {   count += spacing;
                if(count >= 0.02)
                {   *logstream << ".";
                    count -= 0.02;
                }
            }

            double coeff24 = Constant::Electron3j(J2, J4, k, 1, -1);
            if((s2.L() + s4.L() + k)%2)
                coeff24 = 0.;

            if(coeff24)
            {
                coeff24 = coeff24 * (J2 + 1) * (J4 + 1);
                coeff24 = coeff24/(s2.Energy() - s4.Energy() + delta);

                for(i=0; i<mmin(s2.Size(), s4.Size()); i++)
                {
                    density[i] = s2.f[i] * s4.f[i] + Constant::AlphaSquared * s2.g[i] * s4.g[i];
                }
                I.FastCoulombIntegrate(density, Pot24, k, mmin(s2.Size(), s4.Size()));

                double SMS_24 = 0.;
                if(NuclearInverseMass && (k == 1))
                    SMS_24 = -SI.IsotopeShiftIntegral(s4, s2);

                unsigned int start_k1 = (sa.L() + s2.L())%2;
                for(unsigned int k1 = start_k1; k1 <= MAX_K; k1+=2)
                {
                    double coeff = Constant::Electron3j(Ja, J2, k1, 1, -1) *
                                   Constant::Electron3j(J4, Jc, k1, 1, -1) *
                                   Constant::Wigner6j(sa.J(), sc.J(), double(k), s4.J(), s2.J(), double(k1)) *
                                   coeff24 / coeff_ac;
                    if((k1 + k)%2)
                        coeff = -coeff;
                    if((s4.L() + sc.L() + k1)%2)
                        coeff = 0.;

                    if(coeff)
                    {
                        // R1 = R_k1 (a4, 2c)
                        double R1 = 0.;
                        for(i=0; i<mmin(sa.Size(), s2.Size()); i++)
                        {
                            density[i] = sa.f[i] * s2.f[i] + Constant::AlphaSquared * sa.g[i] * s2.g[i];
                        }
                        I.FastCoulombIntegrate(density, Pota2, k1, mmin(sa.Size(), s2.Size()));

                        for(i=0; i < mmin(s4.Size(), sc.Size()); i++)
                            R1 = R1 + Pota2[i] * (s4.f[i] * sc.f[i] + Constant::AlphaSquared * s4.g[i] * sc.g[i]) * dR[i];

                        if(NuclearInverseMass && (k1 == 1))
                        {   R1 = R1 - NuclearInverseMass * SI.IsotopeShiftIntegral(sa, s2) * SI.IsotopeShiftIntegral(s4, sc);
                        }
                        
                        // R2 = R_k (b2, d4)
                        double R2 = 0.;
                        for(i=0; i < mmin(sb.Size(), sd.Size()); i++)
                            R2 = R2 + Pot24[i] * (sb.f[i] * sd.f[i] + Constant::AlphaSquared * sb.g[i] * sd.g[i]) * dR[i];

                        if(NuclearInverseMass && (k == 1))
                        {   R2 = R2 - NuclearInverseMass * SI.IsotopeShiftIntegral(sb, sd) * SMS_24;
                        }

                        energy += R1 * R2 * coeff;
                    }
                }

                // Mirror diagram
                start_k1 = (sb.L() + s2.L())%2;
                for(unsigned int k1 = start_k1; k1 <= MAX_K; k1+=2)
                {
                    double coeff = Constant::Electron3j(Jb, J2, k1, 1, -1) *
                                   Constant::Electron3j(J4, Jd, k1, 1, -1) *
                                   Constant::Wigner6j(sb.J(), sd.J(), double(k), s4.J(), s2.J(), double(k1)) *
                                   coeff24 / coeff_bd;
                    if((k1 + k)%2)
                        coeff = -coeff;
                    if((s4.L() + sd.L() + k1)%2)
                        coeff = 0.;

                    if(coeff)
                    {
                        // R1 = R_k1 (b4, 2d)
                        double R1 = 0.;
                        for(i=0; i<mmin(sb.Size(), s2.Size()); i++)
                        {
                            density[i] = sb.f[i] * s2.f[i] + Constant::AlphaSquared * sb.g[i] * s2.g[i];
                        }
                        I.FastCoulombIntegrate(density, Potb2, k1, mmin(sb.Size(), s2.Size()));

                        for(i=0; i < mmin(s4.Size(), sd.Size()); i++)
                            R1 = R1 + Potb2[i] * (s4.f[i] * sd.f[i] + Constant::AlphaSquared * s4.g[i] * sd.g[i]) * dR[i];

                        if(NuclearInverseMass && (k1 == 1))
                        {   R1 = R1 - NuclearInverseMass * SI.IsotopeShiftIntegral(sb, s2) * SI.IsotopeShiftIntegral(s4, sd);
                        }
                        
                        // R2 = R_k (a2, c4)
                        double R2 = 0.;
                        for(i=0; i < mmin(sa.Size(), sc.Size()); i++)
                            R2 = R2 + Pot24[i] * (sa.f[i] * sc.f[i] + Constant::AlphaSquared * sa.g[i] * sc.g[i]) * dR[i];

                        if(NuclearInverseMass && (k == 1))
                        {   R2 = R2 - NuclearInverseMass * SI.IsotopeShiftIntegral(sa, sc) * SMS_24;
                        }

                        energy += R1 * R2 * coeff;
                    }
                }
            }
            it4.Next();
        }
        it2.Next();
    }

    if(debug)
        *outstream << "  " << energy * Constant::HartreeEnergy_cm << std::endl;
    return energy;
}

double MBPTCalculator::CalculateTwoElectron3(const State& sa, const State& sb, const State& sc, const State& sd, unsigned int k) const
{
    const bool debug = DebugOptions.LogMBPT();
    const double NuclearInverseMass = core->GetNuclearInverseMass();

    unsigned int Ja = (unsigned int)(sa.J() * 2.);
    unsigned int Jb = (unsigned int)(sb.J() * 2.);
    unsigned int Jc = (unsigned int)(sc.J() * 2.);
    unsigned int Jd = (unsigned int)(sd.J() * 2.);

    std::vector<double> density(MaxStateSize);
    std::vector<double> Pot24(MaxStateSize);
    std::vector<double> Pota4(MaxStateSize);
    std::vector<double> Potb4(MaxStateSize);
    CoulombIntegrator I(lattice);
    StateIntegrator SI(lattice);
    const double* dR = lattice->dR();

    if(debug)
        *outstream << "TwoE 3:   ";
    double spacing = 1./double(core->NumStates() * excited->NumStates());
    double count = 0.;
    unsigned int i;

    double energy = 0.;
    const double coeff_ac = Constant::Electron3j(Ja, Jc, k, 1, -1);
    const double coeff_bd = Constant::Electron3j(Jb, Jd, k, 1, -1);
    if(!coeff_ac || !coeff_bd)
        return energy;

    ConstStateIterator it2 = core->GetConstStateIterator();
    while(!it2.AtEnd())
    {   const State& s2 = *(it2.GetState());
        unsigned int J2 = (unsigned int)(s2.J() * 2.);

        ConstStateIterator it4 = excited->GetConstStateIterator();
        while(!it4.AtEnd())
        {   const State& s4 = *(it4.GetState());
            unsigned int J4 = (unsigned int)(s4.J() * 2.);

            if(debug)
            {   count += spacing;
                if(count >= 0.02)
                {   *logstream << ".";
                    count -= 0.02;
                }
            }

            double coeff24 = Constant::Electron3j(J4, J2, k, 1, -1);
            if((s2.L() + s4.L() + k)%2)
                coeff24 = 0.;

            if(coeff24)
            {
                coeff24 = coeff24 * (J2 + 1) * (J4 + 1);
                coeff24 = coeff24/(s2.Energy() - s4.Energy() + delta);

                for(i=0; i<mmin(s2.Size(), s4.Size()); i++)
                {
                    density[i] = s2.f[i] * s4.f[i] + Constant::AlphaSquared * s2.g[i] * s4.g[i];
                }
                I.FastCoulombIntegrate(density, Pot24, k, mmin(s2.Size(), s4.Size()));

                double SMS_42 = 0.;
                if(NuclearInverseMass && (k == 1))
                    SMS_42 = SI.IsotopeShiftIntegral(s4, s2);

                unsigned int start_k1 = (sa.L() + s4.L())%2;
                for(unsigned int k1 = start_k1; k1 <= MAX_K; k1+=2)
                {
                    double coeff = Constant::Electron3j(Ja, J4, k1, 1, -1) *
                                   Constant::Electron3j(J2, Jc, k1, 1, -1) *
                                   Constant::Wigner6j(sa.J(), sc.J(), double(k), s2.J(), s4.J(), double(k1)) *
                                   coeff24 / coeff_ac;
                    if((k1 + k)%2)
                        coeff = -coeff;
                    if((s2.L() + sc.L() + k1)%2)
                        coeff = 0.;

                    if(coeff)
                    {
                        // R1 = R_k1 (a2, 4c)
                        double R1 = 0.;
                        for(i=0; i<mmin(sa.Size(), s4.Size()); i++)
                        {
                            density[i] = sa.f[i] * s4.f[i] + Constant::AlphaSquared * sa.g[i] * s4.g[i];
                        }
                        I.FastCoulombIntegrate(density, Pota4, k1, mmin(sa.Size(), s4.Size()));

                        for(i=0; i < mmin(s2.Size(), sc.Size()); i++)
                            R1 = R1 + Pota4[i] * (s2.f[i] * sc.f[i] + Constant::AlphaSquared * s2.g[i] * sc.g[i]) * dR[i];

                        if(NuclearInverseMass && (k1 == 1))
                        {   R1 = R1 - NuclearInverseMass * SI.IsotopeShiftIntegral(s4, sa) * SI.IsotopeShiftIntegral(sc, s2);
                        }
                        
                        // R2 = R_k (4b, 2d)
                        double R2 = 0.;
                        for(i=0; i < mmin(sb.Size(), sd.Size()); i++)
                            R2 = R2 + Pot24[i] * (sb.f[i] * sd.f[i] + Constant::AlphaSquared * sb.g[i] * sd.g[i]) * dR[i];

                        if(NuclearInverseMass && (k == 1))
                        {   R2 = R2 - NuclearInverseMass * SI.IsotopeShiftIntegral(sb, sd) * SMS_42;
                        }

                        energy += R1 * R2 * coeff;
                    }
                }

                // Mirror diagram
                start_k1 = (sb.L() + s4.L())%2;
                for(unsigned int k1 = start_k1; k1 <= MAX_K; k1+=2)
                {
                    double coeff = Constant::Electron3j(Jb, J4, k1, 1, -1) *
                                   Constant::Electron3j(J2, Jd, k1, 1, -1) *
                                   Constant::Wigner6j(sb.J(), sd.J(), double(k), s2.J(), s4.J(), double(k1)) *
                                   coeff24 / coeff_bd;
                    if((k1 + k)%2)
                        coeff = -coeff;
                    if((s2.L() + sd.L() + k1)%2)
                        coeff = 0.;

                    if(coeff)
                    {
                        // R1 = R_k1 (b2, 4d)
                        double R1 = 0.;
                        for(i=0; i<mmin(sb.Size(), s4.Size()); i++)
                        {
                            density[i] = sb.f[i] * s4.f[i] + Constant::AlphaSquared * sb.g[i] * s4.g[i];
                        }
                        I.FastCoulombIntegrate(density, Potb4, k1, mmin(sb.Size(), s4.Size()));

                        for(i=0; i < mmin(s2.Size(), sd.Size()); i++)
                            R1 = R1 + Potb4[i] * (s2.f[i] * sd.f[i] + Constant::AlphaSquared * s2.g[i] * sd.g[i]) * dR[i];

                        if(NuclearInverseMass && (k1 == 1))
                        {   R1 = R1 - NuclearInverseMass * SI.IsotopeShiftIntegral(s4, sb) * SI.IsotopeShiftIntegral(sd, s2);
                        }
                        
                        // R2 = R_k (4a, 2c)
                        double R2 = 0.;
                        for(i=0; i < mmin(sa.Size(), sc.Size()); i++)
                            R2 = R2 + Pot24[i] * (sa.f[i] * sc.f[i] + Constant::AlphaSquared * sa.g[i] * sc.g[i]) * dR[i];

                        if(NuclearInverseMass && (k == 1))
                        {   R2 = R2 - NuclearInverseMass * SI.IsotopeShiftIntegral(sa, sc) * SMS_42;
                        }

                        energy += R1 * R2 * coeff;
                    }
                }
            }
            it4.Next();
        }
        it2.Next();
    }

    if(debug)
        *outstream << "  " << energy * Constant::HartreeEnergy_cm << std::endl;
    return energy;
}

double MBPTCalculator::CalculateTwoElectron4(const State& sa, const State& sb, const State& sc, const State& sd, unsigned int k) const
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
        *outstream << "TwoE 4/5: ";
    double spacing = 1./double(core->NumStates() * excited->NumStates());
    double count = 0.;
    unsigned int i;

    double energy = 0.;
    const double coeff_ac = Constant::Electron3j(Ja, Jc, k, 1, -1);
    const double coeff_bd = Constant::Electron3j(Jb, Jd, k, 1, -1);
    if(!coeff_ac || !coeff_bd)
        return energy;

    ConstStateIterator it2 = core->GetConstStateIterator();
    while(!it2.AtEnd())
    {   const State& s2 = *(it2.GetState());
        unsigned int J2 = (unsigned int)(s2.J() * 2.);

        ConstStateIterator it4 = excited->GetConstStateIterator();
        while(!it4.AtEnd())
        {   const State& s4 = *(it4.GetState());
            unsigned int J4 = (unsigned int)(s4.J() * 2.);

            if(debug)
            {   count += spacing;
                if(count >= 0.02)
                {   *logstream << ".";
                    count -= 0.02;
                }
            }

            double coeff24 = double(J2 + 1) * double(J4 + 1) * (2. * double(k) + 1.);
            if((sa.L() + s2.L())%2 != (s4.L() + sd.L())%2)
                coeff24 = 0.;
            if((s2.L() + sc.L())%2 != (sb.L() + s4.L())%2)
                coeff24 = 0.;

            if(coeff24)
            {
                coeff24 = coeff24/(coeff_ac*coeff_bd);
                coeff24 = coeff24/(s2.Energy() - s4.Energy() + delta);
                unsigned int exponent = (unsigned int)(Ja + Jb + Jc + Jd + J2 + J4)/2;
                if(exponent%2)
                    coeff24 = -coeff24;

                unsigned int start_k1 = (sa.L() + s2.L())%2;
                for(unsigned int k1 = start_k1; k1 <= MAX_K; k1+=2)
                {
                    double coeff_ad = Constant::Electron3j(Ja, J2, k1, 1, -1) *
                                      Constant::Electron3j(J4, Jd, k1, 1, -1);

                    if(coeff_ad)
                    {
                        // R1 = R_k1 (a4, 2d)
                        double R1 = 0.;
                        for(i=0; i<mmin(s4.Size(), sd.Size()); i++)
                        {
                            density[i] = s4.f[i] * sd.f[i] + Constant::AlphaSquared * s4.g[i] * sd.g[i];
                        }

                        I.FastCoulombIntegrate(density, pot, k1, mmin(s4.Size(), sd.Size()));
                        for(i=0; i < mmin(sa.Size(), s2.Size()); i++)
                            R1 = R1 + pot[i] * (sa.f[i] * s2.f[i] + Constant::AlphaSquared * sa.g[i] * s2.g[i]) * dR[i];

                        if(NuclearInverseMass && (k1 == 1))
                        {   R1 = R1 - NuclearInverseMass * SI.IsotopeShiftIntegral(sa, s2) * SI.IsotopeShiftIntegral(s4, sd);
                        }

                        unsigned int start_k2 = (s2.L() + sc.L())%2;
                        for(unsigned int k2 = start_k2; k2 <= MAX_K; k2+=2)
                        {
                            double coeff = Constant::Electron3j(J2, Jc, k2, 1, -1) *
                                           Constant::Electron3j(Jb, J4, k2, 1, -1);
                            if(coeff)
                                coeff = coeff * Constant::Wigner6j(sc.J(), sa.J(), k, k1, k2, s2.J())
                                              * Constant::Wigner6j(sb.J(), sd.J(), k, k1, k2, s4.J());

                            if(coeff)
                            {   
                                coeff = coeff * coeff_ad * coeff24;

                                // R2 = R_k2 (2b, c4)
                                double R2 = 0.;
                                for(i=0; i<mmin(sb.Size(), s4.Size()); i++)
                                {
                                    density[i] = sb.f[i] * s4.f[i] + Constant::AlphaSquared * sb.g[i] * s4.g[i];
                                }
                                I.FastCoulombIntegrate(density, pot, k2, mmin(sb.Size(), s4.Size()));
                                
                                for(i=0; i < mmin(s2.Size(), sc.Size()); i++)
                                    R2 = R2 + pot[i] * (s2.f[i] * sc.f[i] + Constant::AlphaSquared * s2.g[i] * sc.g[i]) * dR[i];

                                if(NuclearInverseMass && (k2 == 1))
                                {   R2 = R2 - NuclearInverseMass * SI.IsotopeShiftIntegral(sc, s2) * SI.IsotopeShiftIntegral(s4, sb);
                                }

                                energy += R1 * R2 * coeff;
                            }
                        }
                    }
                }
            }
            it4.Next();
        }
        it2.Next();
    }

    if(debug)
        *outstream << "  " << energy * Constant::HartreeEnergy_cm << std::endl;
    return energy;
}

double MBPTCalculator::CalculateTwoElectron5(const State& sa, const State& sb, const State& sc, const State& sd, unsigned int k) const
{
    return CalculateTwoElectron4(sb, sa, sd, sc, k);
}

double MBPTCalculator::CalculateTwoElectron6(const State& sa, const State& sb, const State& sc, const State& sd, unsigned int k) const
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
        *outstream << "TwoE 6:   ";
    double spacing = 1./double(core->NumStates() * core->NumStates());
    double count = 0.;
    unsigned int i;

    double energy = 0.;
    const double coeff_ac = Constant::Electron3j(Ja, Jc, k, 1, -1);
    const double coeff_bd = Constant::Electron3j(Jb, Jd, k, 1, -1);
    if(!coeff_ac || !coeff_bd)
        return energy;

    const double ValenceEnergy = ValenceEnergies.find(sa.Kappa())->second + ValenceEnergies.find(sb.Kappa())->second;

    ConstStateIterator itm = core->GetConstStateIterator();
    while(!itm.AtEnd())
    {   const State& sm = *(itm.GetState());
        unsigned int Jm = (unsigned int)(sm.J() * 2.);

        ConstStateIterator itn = core->GetConstStateIterator();
        while(!itn.AtEnd())
        {   const State& sn = *(itn.GetState());
            unsigned int Jn = (unsigned int)(sn.J() * 2.);

            if(debug)
            {   count += spacing;
                if(count >= 0.02)
                {   *logstream << ".";
                    count -= 0.02;
                }
            }

            double coeff_mn = double(Jm + 1) * double(Jn + 1) * (2. * double(k) + 1.);
            if((sa.L() + sm.L())%2 != (sb.L() + sn.L())%2)
                coeff_mn = 0.;
            if((sm.L() + sc.L())%2 != (sd.L() + sn.L())%2)
                coeff_mn = 0.;

            if(coeff_mn)
            {
                coeff_mn = coeff_mn/(coeff_ac*coeff_bd);
                coeff_mn = coeff_mn/(sm.Energy() + sn.Energy() - ValenceEnergy + delta);
                unsigned int exponent = (unsigned int)(Ja + Jb + Jc + Jd + Jm + Jn)/2;
                if((exponent + k + 1)%2)
                    coeff_mn = -coeff_mn;

                unsigned int start_k1 = (sa.L() + sm.L())%2;
                for(unsigned int k1 = start_k1; k1 <= MAX_K; k1+=2)
                {
                    double coeff_ab = Constant::Electron3j(Ja, Jm, k1, 1, -1) *
                                      Constant::Electron3j(Jb, Jn, k1, 1, -1);

                    if(coeff_ab)
                    {
                        // R1 = R_k1 (ab, mn)
                        double R1 = 0.;
                        for(i=0; i<mmin(sb.Size(), sn.Size()); i++)
                        {
                            density[i] = sb.f[i] * sn.f[i] + Constant::AlphaSquared * sb.g[i] * sn.g[i];
                        }
                        I.FastCoulombIntegrate(density, pot, k1, mmin(sb.Size(), sn.Size()));
                        
                        for(i=0; i < mmin(sa.Size(), sm.Size()); i++)
                            R1 = R1 + pot[i] * (sa.f[i] * sm.f[i] + Constant::AlphaSquared * sa.g[i] * sm.g[i]) * dR[i];

                        if(NuclearInverseMass && (k1 == 1))
                        {   R1 = R1 - NuclearInverseMass * SI.IsotopeShiftIntegral(sa, sm) * SI.IsotopeShiftIntegral(sb, sn);
                        }

                        unsigned int start_k2 = (sm.L() + sc.L())%2;
                        for(unsigned int k2 = start_k2; k2 <= MAX_K; k2+=2)
                        {
                            double coeff = Constant::Electron3j(Jm, Jc, k2, 1, -1) *
                                           Constant::Electron3j(Jn, Jd, k2, 1, -1);
                            if(coeff)
                                coeff = coeff * Constant::Wigner6j(sc.J(), sa.J(), k, k1, k2, sm.J())
                                              * Constant::Wigner6j(sd.J(), sb.J(), k, k1, k2, sn.J());

                            if(coeff)
                            {   
                                coeff = coeff * coeff_ab * coeff_mn;
                                if((k1 + k2)%2)
                                    coeff = -coeff;

                                // R2 = R_k2 (mn, cd)
                                double R2 = 0.;
                                for(i=0; i<mmin(sn.Size(), sd.Size()); i++)
                                {
                                    density[i] = sn.f[i] * sd.f[i] + Constant::AlphaSquared * sn.g[i] * sd.g[i];
                                }
                                I.FastCoulombIntegrate(density, pot, k2, mmin(sn.Size(), sd.Size()));
                                
                                for(i=0; i < mmin(sm.Size(), sc.Size()); i++)
                                    R2 = R2 + pot[i] * (sm.f[i] * sc.f[i] + Constant::AlphaSquared * sm.g[i] * sc.g[i]) * dR[i];

                                if(NuclearInverseMass && (k2 == 1))
                                {   R2 = R2 - NuclearInverseMass * SI.IsotopeShiftIntegral(sc, sm) * SI.IsotopeShiftIntegral(sd, sn);
                                }

                                energy += R1 * R2 * coeff;
                            }
                        }
                    }
                }
            }
            itn.Next();
        }
        itm.Next();
    }

    if(debug)
        *outstream << "  " << energy * Constant::HartreeEnergy_cm << std::endl;
    return energy;
}

double MBPTCalculator::CalculateTwoElectronSub(const State& sa, const State& sb, const State& sc, const State& sd, unsigned int k) const
{
    const bool debug = DebugOptions.LogMBPT();
    const double NuclearInverseMass = core->GetNuclearInverseMass();

    std::vector<double> density(MaxStateSize);
    std::vector<double> pot(MaxStateSize);

    CoulombIntegrator I(lattice);
    StateIntegrator SI(lattice);
    const double* dR = lattice->dR();

    if(debug)
        *outstream << "2eSub :   ";

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

    ConstStateIterator itn = core->GetConstStateIterator();
    while(!itn.AtEnd())
    {   const State& sn = *(itn.GetState());

        if(sn.Kappa() == sa.Kappa())
        {
            double R1 = 0.;
            for(i=0; i < mmin(sn.Size(), sc.Size()); i++)
                R1 += pot[i] * (sn.f[i] * sc.f[i] + Constant::AlphaSquared * sn.g[i] * sc.g[i]) * dR[i];

            if(SMS_bd)
            {   R1 = R1 + SMS_bd * SI.IsotopeShiftIntegral(sc, sn);
            }

            energy -= R1 * SI.HamiltonianMatrixElement(sa, sn, *core) / (sn.Energy() - ValenceEnergies.find(sa.Kappa())->second + delta);
        }

        if(sn.Kappa() == sc.Kappa())
        {
            double R1 = 0.;
            for(i=0; i < mmin(sa.Size(), sn.Size()); i++)
                R1 += pot[i] * (sa.f[i] * sn.f[i] + Constant::AlphaSquared * sa.g[i] * sn.g[i]) * dR[i];

            if(SMS_bd)
            {   R1 = R1 - SMS_bd * SI.IsotopeShiftIntegral(sa, sn);
            }

            energy -= R1 * SI.HamiltonianMatrixElement(sc, sn, *core) / (sn.Energy() - ValenceEnergies.find(sc.Kappa())->second + delta);
        }

        itn.Next();
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

    itn.First();
    while(!itn.AtEnd())
    {   const State& sn = *(itn.GetState());

        if(sn.Kappa() == sb.Kappa())
        {
            double R1 = 0.;
            for(i=0; i < mmin(sn.Size(), sd.Size()); i++)
                R1 += pot[i] * (sn.f[i] * sd.f[i] + Constant::AlphaSquared * sn.g[i] * sd.g[i]) * dR[i];

            if(SMS_ac)
            {   R1 = R1 + SMS_ac * SI.IsotopeShiftIntegral(sd, sn);
            }

            energy -= R1 * SI.HamiltonianMatrixElement(sb, sn, *core) / (sn.Energy() - ValenceEnergies.find(sb.Kappa())->second + delta);
        }

        if(sn.Kappa() == sd.Kappa())
        {
            double R1 = 0.;
            for(i=0; i < mmin(sb.Size(), sn.Size()); i++)
                R1 += pot[i] * (sb.f[i] * sn.f[i] + Constant::AlphaSquared * sb.g[i] * sn.g[i]) * dR[i];

            if(SMS_ac)
            {   R1 = R1 - SMS_ac * SI.IsotopeShiftIntegral(sb, sn);
            }

            energy -= R1 * SI.HamiltonianMatrixElement(sd, sn, *core) / (sn.Energy() - ValenceEnergies.find(sd.Kappa())->second + delta);
        }

        itn.Next();
    }

    if(debug)
        *outstream << "  " << energy * Constant::HartreeEnergy_cm << std::endl;
    return energy;
}

void MBPTCalculator::SetValenceEnergies()
{
    ConstStateIterator it_i = excited->GetConstStateIterator();
    ValenceEnergies.clear();

    // Get maximum angular momentum in excited states
    unsigned int max_l = 0;
    it_i.First();
    while(!it_i.AtEnd())
    {   max_l = mmax(it_i.GetState()->L(), max_l);
        it_i.Next();
    }

    for(int kappa = - (int)max_l - 1; kappa <= (int)max_l; kappa++)
        if(kappa != 0)
        {
            double valence_energy = 0.;
            unsigned int pqn = 10;

            // Get leading state (for energy denominator)
            it_i.First();
            while(!it_i.AtEnd())
            {   const DiscreteState* ds = it_i.GetState();
                if((ds->Kappa() == kappa) && (ds->RequiredPQN() < pqn))
                {   pqn = ds->RequiredPQN();
                    valence_energy = ds->Energy();
                }
                it_i.Next();
            }

            ValenceEnergies.insert(std::pair<int, double>(kappa, valence_energy));
        }
}

#include "Include.h"
#include "MBPTCalculator.h"
#include "Universal/Constant.h"
#include "Universal/CoulombIntegrator.h"
#include "HartreeFock/StateIntegrator.h"

MBPTCalculator::MBPTCalculator(Lattice* lat, const Core* atom_core, ExcitedStates* excited_states):
    lattice(lat), core(atom_core), excited(excited_states)
{}

double MBPTCalculator::GetSecondOrderSigma(const State* s, SigmaPotential* sigma)
{
    if(sigma != NULL)
        sigma->Reset();

    double ten1 = CalculateCorrelation1and3(s, sigma);
    double ten2 = CalculateCorrelation2(s, sigma);
    double ten4 = CalculateCorrelation4(s, sigma);

    double energy = 0.;
    if(sigma != NULL)
    {   std::vector<double> sum = sigma->GetPotential(s->f);
        for(unsigned int i=0; i<sum.size(); i++)
            energy += s->f[i] * sum[i] * lattice->dR(i);

        *outstream << "\tSecond order contribution: " << energy * Constant::HartreeEnergy_cm << std::endl;
    }

    *outstream << "Total energy from direct summation = " << (s->Energy() + ten1 + ten2 + ten4)*Constant::HartreeEnergy_cm << std::endl;
    *outstream << "Total energy via sigma potential   = " << (s->Energy() + energy) * Constant::HartreeEnergy_cm << std::endl;

    return (s->Energy() + ten1 + ten2 + ten4);
}

double MBPTCalculator::CalculateCorrelation1and3(const State* s, SigmaPotential* sigma)
{
    double NuclearInverseMass = core->GetNuclearInverseMass();

    std::vector<double> Pot24;
    CoulombIntegrator I(*lattice);
    StateIntegrator SI(*lattice);

    *outstream << "Cor 1+3:  ";
    double spacing = 1./double(core->NumStates() * excited->NumStates());
    double count = 0;
    unsigned int i;

    double energy1 = 0., energy3 = 0.;

    /* Firstly, get the loop 24 */
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

            for(unsigned int k=start_k; k<=8; k+=2)
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
                        density[i] = s2.f[i] * s4.f[i];
                    }
                    density.resize(core->GetHFPotential().size());
                    I.FastCoulombIntegrate(density, Pot24, k);

                    // Isotope shift
                    double SMS_24 = 0.;
                    if(NuclearInverseMass && (k == 1))
                        SMS_24 = -SI.IsotopeShiftIntegral(s4, s2);

                    // Correlation 1 has excited state 3
                    ConstStateIterator it3_1 = excited->GetConstStateIterator();
                    while(!it3_1.AtEnd())
                    {   const State& s3 = *(it3_1.GetState());

                        if((s->L() + s3.L())%2 == start_k)
                        {
                            double coeff = Constant::Electron3j(s->J(), s3.J(), k, 0.5, -0.5);

                            if(coeff)
                            {
                                coeff = coeff * coeff * coeff24 * (2. * s3.J() + 1.);
                                coeff = coeff/(s->Energy() + s2.Energy() - s3.Energy() - s4.Energy());
                                coeff = coeff * it3_1.Weight();

                                double de = 0.;
                                unsigned int upper_bound = mmin(s3.Size(), s->Size());
                                std::vector<double> Y(upper_bound);

                                for(i=0; i<upper_bound; i++)
                                {   Y[i] = Pot24[i] * s3.f[i];
                                    de = de + Y[i] * s->f[i] * lattice->dR(i);
                                }

                                if(SMS_24)
                                {
                                    std::vector<double> P(upper_bound);
                                    double de_sms = SI.IsotopeShiftIntegral(*s, s3, &P);

                                    for(i=0; i<upper_bound; i++)
                                        Y[i] = Y[i] - NuclearInverseMass * SMS_24 * P[i];

                                    de_sms = de_sms * SMS_24 * NuclearInverseMass;
                                    de = de - de_sms;
                                }

                                if(sigma)
                                    sigma->AddToSigma(Y, Y, coeff);
                                energy1 = energy1 + de * de * coeff;
                            }
                        }
                        it3_1.Next();
                    }

                    // Correlation 3 has core state 3
                    ConstStateIterator it3_3 = core->GetConstStateIterator();
                    while(!it3_3.AtEnd())
                    {   const State& s3 = *(it3_3.GetState());

                        if((s->L() + s3.L())%2 == start_k)
                        {
                            double coeff = Constant::Electron3j(s->J(), s3.J(), k, 0.5, -0.5);

                            if(coeff)
                            {
                                coeff = coeff * coeff * coeff24 * (2. * s3.J() + 1.);
                                coeff = coeff/(s->Energy() + s4.Energy() - s2.Energy() - s3.Energy());

                                double de = 0.;
                                unsigned int upper_bound = mmin(s3.Size(), s->Size());
                                std::vector<double> Y(upper_bound);

                                for(i=0; i<upper_bound; i++)
                                {   Y[i] = Pot24[i] * s3.f[i];
                                    de = de + Y[i] * s->f[i] * lattice->dR(i);
                                }

                                if(SMS_24)
                                {
                                    std::vector<double> P(upper_bound);
                                    double de_sms = -SI.IsotopeShiftIntegral(*s, s3, &P);

                                    for(i=0; i<upper_bound; i++)
                                        // Plus sign is because P is opposite sign
                                        Y[i] = Y[i] + NuclearInverseMass * SMS_24 * P[i];
                                    
                                    de_sms = de_sms * SMS_24 * NuclearInverseMass;
                                    de = de - de_sms;
                                }

                                if(sigma)
                                    sigma->AddToSigma(Y, Y, coeff);
                                energy3 = energy3 + de * de * coeff;
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

    *outstream << "  " << energy1 * Constant::HartreeEnergy_cm
              << "  " << energy3 * Constant::HartreeEnergy_cm << std::endl;

    return energy1 + energy3;
}

double MBPTCalculator::CalculateCorrelation2(const State* s, SigmaPotential* sigma)
{
    double NuclearInverseMass = core->GetNuclearInverseMass();

    std::vector<double> Pot24, Pot23;
    CoulombIntegrator I(*lattice);
    StateIntegrator SI(*lattice);
    const double* dR = lattice->dR();

    *outstream << "Cor 2:    ";
    double spacing = 1./double(core->NumStates() * excited->NumStates());
    double count = 0;
    unsigned int i;

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

            for(unsigned int k1=start_k1; k1<=8; k1+=2)
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

                    // Isotope shift
                    double SMS_24 = 0.;
                    if(NuclearInverseMass && (k1 == 1))
                        SMS_24 = -SI.IsotopeShiftIntegral(s4, s2);

                    ConstStateIterator it3 = excited->GetConstStateIterator();
                    while(!it3.AtEnd())
                    {   const State& s3 = *(it3.GetState());

                        double coeff13 = Constant::Electron3j(s->J(), s3.J(), k1, 0.5, -0.5);
                        unsigned int start_k2 = (s2.L() + s3.L())%2;
                        if(coeff13 && ((s->L() + s4.L())%2 == start_k2) && ((s->L() + s3.L())%2 == start_k1))
                        {
                            coeff13 = coeff13 * (2. * s3.J() + 1.) * it3.Weight()
                                      /(s->Energy() + s2.Energy() - s3.Energy() - s4.Energy());
                            // Sign
                            if((start_k1 + start_k2)%2)
                                coeff13 = -coeff13;

                            for(unsigned int k2=start_k2; k2<=8; k2+=2)
                            {
                                double coeff 
                                    = coeff13 * coeff24 * Constant::Electron3j(s->J(), s4.J(), k2, 0.5, -0.5)
                                    * Constant::Electron3j(s3.J(), s2.J(), k2, 0.5, -0.5)
                                    * Constant::Wigner6j(s->J(), s3.J(), k1, s2.J(), s4.J(), k2);

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

                                    // Isotope shift
                                    double SMS_23 = 0.;
                                    if(NuclearInverseMass && (k2 == 1))
                                        SMS_23 = SI.IsotopeShiftIntegral(s3, s2);

                                    double de1 = 0.;
                                    unsigned int upper_bound = mmin(s3.Size(), s->Size());
                                    std::vector<double> Y1(upper_bound);
                                    for(i=0; i<upper_bound; i++)
                                    {
                                        Y1[i] = Pot24[i] * s3.f[i];
                                        de1 = de1 + Y1[i] * s->f[i] * dR[i];
                                    }
                                
                                    if(SMS_24)
                                    {   
                                        std::vector<double> P1(upper_bound);
                                        double de1_sms = SI.IsotopeShiftIntegral(*s, s3, &P1);
                                        
                                        for(i=0; i<upper_bound; i++)
                                            Y1[i] = Y1[i] - NuclearInverseMass * SMS_24 * P1[i];

                                        de1_sms = de1_sms * SMS_24 * NuclearInverseMass;
                                        de1 = de1 - de1_sms;
                                    }

                                    double de2 = 0.;
                                    upper_bound = mmin(s4.Size(), s->Size());
                                    std::vector<double> Y2(upper_bound);
                                    for(i=0; i<upper_bound; i++)
                                    {
                                        Y2[i] = Pot23[i] * s4.f[i];
                                        de2 = de2 + Y2[i] * s->f[i] * dR[i];
                                    }

                                    if(SMS_23)
                                    {   
                                        std::vector<double> P2(upper_bound);
                                        double de2_sms = -SI.IsotopeShiftIntegral(*s, s4, &P2);

                                        for(i=0; i<upper_bound; i++)
                                            // P2 is wrong sign
                                            Y2[i] = Y2[i] + NuclearInverseMass * SMS_23 * P2[i];

                                        de2_sms = de2_sms * SMS_23 * NuclearInverseMass;
                                        de2 = de2 - de2_sms;
                                    }

                                    if(sigma)
                                        sigma->AddToSigma(Y1, Y2, coeff);
                                    energy = energy + de1 * de2 * coeff;
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

    *outstream << "  " << energy * Constant::HartreeEnergy_cm << std::endl;
    return energy;
}

double MBPTCalculator::CalculateCorrelation4(const State* s, SigmaPotential* sigma)
{
    double NuclearInverseMass = core->GetNuclearInverseMass();

    std::vector<double> Pot24, Pot34;
    CoulombIntegrator I(*lattice);
    StateIntegrator SI(*lattice);

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

            for(unsigned int k1=start_k1; k1<=8; k1+=2)
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

                    // Isotope shift
                    double SMS_24 = 0.;
                    if(NuclearInverseMass && (k1 == 1))
                        SMS_24 = -SI.IsotopeShiftIntegral(s4, s2);

                    ConstStateIterator it3 = core->GetConstStateIterator();
                    while(!it3.AtEnd())
                    {   const State& s3 = *(it3.GetState());

                        double coeff13 = Constant::Electron3j(s->J(), s3.J(), k1, 0.5, -0.5);
                        unsigned int start_k2 = (s3.L() + s4.L())%2;
                        if(coeff13 && ((s->L() + s2.L())%2 == start_k2) && ((s->L() + s3.L())%2 == start_k1))
                        {
                            coeff13 = coeff13 * (2. * s3.J() + 1.)
                                      /(s->Energy() + s4.Energy() - s2.Energy() - s3.Energy());
                            // Sign
                            if((start_k1 + start_k2)%2)
                                coeff13 = -coeff13;

                            for(unsigned int k2= start_k2; k2<=8; k2+=2)
                            {
                                double coeff 
                                    = coeff13 * coeff24 * Constant::Electron3j(s->J(), s2.J(), k2, 0.5, -0.5)
                                    * Constant::Electron3j(s3.J(), s4.J(), k2, 0.5, -0.5)
                                    * Constant::Wigner6j(s->J(), s3.J(), k1, s4.J(), s2.J(), k2);

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

                                    // Isotope shift
                                    double SMS_34 = 0.;
                                    if(NuclearInverseMass && (k2 == 1))
                                        SMS_34 = SI.IsotopeShiftIntegral(s4, s3);

                                    double de1 = 0.;
                                    unsigned int upper_bound = mmin(s3.Size(), s->Size());
                                    std::vector<double> Y1(upper_bound);
                                    for(i=0; i<upper_bound; i++)
                                    {
                                        Y1[i] = Pot24[i] * s3.f[i];
                                        de1 = de1 + Y1[i] * s->f[i] * lattice->dR(i);
                                    }

                                    if(SMS_24)
                                    {   
                                        std::vector<double> P1(upper_bound);
                                        double de1_sms = -SI.IsotopeShiftIntegral(*s, s3, &P1);

                                        for(i=0; i<upper_bound; i++)
                                            Y1[i] = Y1[i] + NuclearInverseMass * SMS_24 * P1[i];

                                        de1_sms = de1_sms * SMS_24 * NuclearInverseMass;
                                        de1 = de1 - de1_sms;
                                    }

                                    double de2 = 0.;
                                    upper_bound = mmin(s2.Size(), s->Size());
                                    std::vector<double> Y2(upper_bound);
                                    for(i=0; i<upper_bound; i++)
                                    {
                                        Y2[i] = Pot34[i] * s2.f[i];
                                        de2 = de2 + Y2[i] * s->f[i] * lattice->dR(i);
                                    }

                                    if(SMS_34)
                                    {   
                                        std::vector<double> P2(upper_bound);
                                        double de2_sms = SI.IsotopeShiftIntegral(*s, s2, &P2);

                                        for(i=0; i<upper_bound; i++)
                                            Y2[i] = Y2[i] - NuclearInverseMass * SMS_34 * P2[i];

                                        de2_sms = de2_sms * SMS_34 * NuclearInverseMass;
                                        de2 = de2 - de2_sms;
                                    }

                                    if(sigma)
                                        sigma->AddToSigma(Y1, Y2, coeff);
                                    energy = energy + de1 * de2 * coeff;
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

    *outstream << "  " << energy * Constant::HartreeEnergy_cm << std::endl;
    return energy;
}

#include "MBPTCalculator.h"
#include "Universal/Constant.h"
#include "Universal/CoulombIntegrator.h"
#include "HartreeFock/StateIntegrator.h"

MBPTCalculator::MBPTCalculator(Lattice* lat, const Core* atom_core, const ExcitedStates* excited_states):
    lattice(lat), core(atom_core), excited(excited_states), integrals(NULL), delta(0.)
{
    SetValenceEnergies();
}

MBPTCalculator::~MBPTCalculator(void)
{
    if(integrals)
        delete integrals;
}

unsigned int MBPTCalculator::GetStorageSize(const ExcitedStates* valence_states)
{
    if(!integrals)
        integrals = new SlaterIntegrals(excited);
 
    return integrals->GetStorageSize(*valence_states);
}

void MBPTCalculator::UpdateIntegrals(const ExcitedStates* valence_states)
{
    SetValenceEnergies();
    if(integrals)
        delete integrals;
    
    integrals = new SlaterIntegrals(excited);
    integrals->Update(*valence_states);
}

void MBPTCalculator::GetSecondOrderSigma(int kappa, SigmaPotential* sigma)
{
    if(DebugOptions.LogMBPT())
        *outstream << "\nkappa = " << kappa << std::endl;

    sigma->Reset();

    CalculateCorrelation1and3(kappa, sigma);
    CalculateCorrelation2(kappa, sigma);
    CalculateCorrelation4(kappa, sigma);
}

double MBPTCalculator::GetOneElectronDiagrams(const StateInfo& s1, const StateInfo& s2) const
{
    if(s1.Kappa() != s2.Kappa())
        return 0.;

    if(DebugOptions.LogMBPT())
        *outstream << "\n < " << s1.Name() << " | " << s2.Name() << " > :" << std::endl;

    double term1 = CalculateCorrelation1and3(s1, s2);
    double term2 = CalculateCorrelation2(s1, s2);
    double term4 = CalculateCorrelation4(s1, s2);

    return (term1 + term2 + term4);
}

double MBPTCalculator::GetOneElectronSubtraction(const StateInfo& s1, const StateInfo& s2) const
{
    if(s1.Kappa() != s2.Kappa())
        return 0.;

    if(DebugOptions.LogMBPT())
        *outstream << "\n < " << s1.Name() << " | " << s2.Name() << " > :" << std::endl;

    double term1 = CalculateSubtraction1(s1, s2);
    double term2 = CalculateSubtraction2(s1, s2);
    double term3 = CalculateSubtraction3(s1, s2);
    
    return (term1 + term2 + term3);
}

double MBPTCalculator::GetTwoElectronDiagrams(unsigned int k, const StateInfo& s1, const StateInfo& s2, const StateInfo& s3, const StateInfo& s4) const
{
    if(DebugOptions.LogMBPT())
        *outstream << "\n R^" << k << " ( " << s1.Name() << " " << s2.Name()
                   << ", " << s3.Name() << " " << s4.Name() << ") :" << std::endl;

    double term = 0;
    term += CalculateTwoElectron1(k, s1, s2, s3, s4);
    term += CalculateTwoElectron2(k, s1, s2, s3, s4);
    term += CalculateTwoElectron3(k, s1, s2, s3, s4);
    term += CalculateTwoElectron4(k, s1, s2, s3, s4);
    term += CalculateTwoElectron5(k, s1, s2, s3, s4);
    term += CalculateTwoElectron6(k, s1, s2, s3, s4);

    return term;
}

double MBPTCalculator::GetTwoElectronBoxDiagrams(unsigned int k, const StateInfo& s1, const StateInfo& s2, const StateInfo& s3, const StateInfo& s4) const
{
    if(DebugOptions.LogMBPT())
        *outstream << "\n R^" << k << " ( " << s1.Name() << " " << s2.Name()
                   << ", " << s3.Name() << " " << s4.Name() << ") :" << std::endl;

    double term = 0;
    term += CalculateTwoElectron4(k, s1, s2, s3, s4);
    term += CalculateTwoElectron5(k, s1, s2, s3, s4);
    term += CalculateTwoElectron6(k, s1, s2, s3, s4);

    return term;
}

double MBPTCalculator::GetTwoElectronSubtraction(unsigned int k, const StateInfo& s1, const StateInfo& s2, const StateInfo& s3, const StateInfo& s4) const
{
    if(DebugOptions.LogMBPT())
        *outstream << "\n R^" << k << " ( " << s1.Name() << " " << s2.Name()
                   << ", " << s3.Name() << " " << s4.Name() << ") :" << std::endl;

    return CalculateTwoElectronSub(k, s1, s2, s3, s4);
}

void MBPTCalculator::CalculateCorrelation1and3(int kappa, SigmaPotential* sigma) const
{
    const bool debug = DebugOptions.LogMBPT();
    unsigned int MaxStateSize = core->GetConstHFPotential().size();
    const double NuclearInverseMass = core->GetNuclearInverseMass();

    StateInfo sa(100, kappa);

    std::vector<double> density(MaxStateSize);
    std::vector<double> P_nalpha(MaxStateSize);
    std::vector<double> Y(MaxStateSize); unsigned int Y_size;
    CoulombIntegrator I(lattice);
    StateIntegrator SI(lattice);

    if(debug)
        *outstream << "Cor 1+3:  ";
    double spacing = 1./double(core->NumStates() * excited->NumStates());
    double count = 0.;
    unsigned int i;

    const double ValenceEnergy = ValenceEnergies.find(kappa)->second;

    unsigned int k1, k1max;

    // Firstly, get the loop 24
    ConstStateIterator it_n = core->GetConstStateIterator();
    while(!it_n.AtEnd())
    {
        const State& sn = *(it_n.GetState());

        ConstStateIterator it_alpha = excited->GetConstStateIterator();
        while(!it_alpha.AtEnd())
        {
            const State& salpha = *(it_alpha.GetState());

            if(debug)
            {   count += spacing;
                if(count >= 0.02)
                {   *logstream << ".";
                    count -= 0.02;
                }
            }

            k1 = absdiff(sn.L(), salpha.L());
            if(absdiff(sn.TwoJ(), salpha.TwoJ()) > 2 * k1)
                k1 += 2;

            k1max = (sn.TwoJ() + salpha.TwoJ())/2;

            while(k1 <= k1max)
            {
                double C_nalpha = Constant::Electron3j(sn.TwoJ(), salpha.TwoJ(), k1, 1, -1);

                if(C_nalpha)
                {
                    C_nalpha = C_nalpha * C_nalpha * (sn.TwoJ() + 1) * (salpha.TwoJ() + 1)
                                                / (2. * k1 + 1.);
                    C_nalpha = C_nalpha * it_alpha.Weight();

                    for(i=0; i<mmin(sn.Size(), salpha.Size()); i++)
                    {
                        density[i] = sn.f[i] * salpha.f[i] + Constant::AlphaSquared * sn.g[i] * salpha.g[i];
                    }
                    I.FastCoulombIntegrate(density, P_nalpha, k1, mmin(sn.Size(), salpha.Size()));

                    double SMS_nalpha = 0.;
                    if(NuclearInverseMass && (k1 == 1))
                        SMS_nalpha = -SI.IsotopeShiftIntegral(salpha, sn);

                    // Correlation 1 has excited state beta
                    ConstStateIterator it_beta = excited->GetConstStateIterator();
                    while(!it_beta.AtEnd())
                    {
                        const State& sbeta = *(it_beta.GetState());

                        double coeff;
                        if((sa.L() + sbeta.L())%2 == k1%2)
                            coeff = Constant::Electron3j(sa.TwoJ(), sbeta.TwoJ(), k1, 1, -1);
                        else
                            coeff = 0.;

                        if(coeff)
                        {
                            coeff = coeff * coeff * C_nalpha * (sbeta.TwoJ() + 1);
                            coeff = coeff/(ValenceEnergy + sn.Energy() - sbeta.Energy() - salpha.Energy() + delta);

                            // R1 = R_k1 (a n, beta alpha)
                            // R2 = R_k1 (b n, beta alpha)
                            Y_size = sbeta.Size();
                            for(i=0; i < Y_size; i++)
                                Y[i] = P_nalpha[i] * sbeta.f[i];

                            if(SMS_nalpha)
                            {
                                std::vector<double> P(Y_size);
                                SI.IsotopeShiftIntegral(sa.L(), sbeta, &P);
                                for(i=0; i < Y_size; i++)
                                    Y[i] = Y[i] - NuclearInverseMass * SMS_nalpha * P[i];
                            }

                            sigma->AddToSigma(Y, Y, coeff, Y_size, Y_size);
                        }
                        it_beta.Next();
                    }

                    // Correlation 3 has core state m
                    ConstStateIterator it_m = core->GetConstStateIterator();
                    while(!it_m.AtEnd())
                    {
                        const State& sm = *(it_m.GetState());

                        double coeff;
                        if((sa.L() + sm.L())%2 == k1%2)
                            coeff =  Constant::Electron3j(sa.TwoJ(), sm.TwoJ(), k1, 1, -1);
                        else
                            coeff = 0.;

                        if(coeff)
                        {
                            coeff = coeff * coeff * C_nalpha * (sm.TwoJ() + 1);
                            coeff = coeff/(ValenceEnergy + salpha.Energy() - sn.Energy() - sm.Energy() - delta);

                            // R1 = R_k1 (a alpha, m n)
                            // R2 = R_k1 (b alpha, m n)
                            Y_size = sm.Size();
                            for(i=0; i < Y_size; i++)
                                Y[i] = P_nalpha[i] * sm.f[i];

                            if(SMS_nalpha)
                            {   std::vector<double> P(Y_size);
                                SI.IsotopeShiftIntegral(sa.L(), sm, &P);
                                for(i=0; i < Y_size; i++)
                                    Y[i] = Y[i] + NuclearInverseMass * SMS_nalpha * P[i];
                            }

                            sigma->AddToSigma(Y, Y, coeff, Y_size, Y_size);
                        }
                        it_m.Next();
                    }
                } // C_nalpha
                k1 += 2;
            }
            it_alpha.Next();
        }
        it_n.Next();
    }
}

void MBPTCalculator::CalculateCorrelation2(int kappa, SigmaPotential* sigma) const
{
    const bool debug = DebugOptions.LogMBPT();
    unsigned int MaxStateSize = core->GetConstHFPotential().size();
    const double NuclearInverseMass = core->GetNuclearInverseMass();

    StateInfo sa(100, kappa);

    std::vector<double> density(MaxStateSize);
    std::vector<double> P_nalpha(MaxStateSize);
    std::vector<double> P_nbeta(MaxStateSize);
    std::vector<double> Y1(MaxStateSize); unsigned int Y1_size;
    std::vector<double> Y2(MaxStateSize); unsigned int Y2_size;
    CoulombIntegrator I(lattice);
    StateIntegrator SI(lattice);

    if(debug)
        *outstream << "Cor 2:    ";
    double spacing = 1./double(core->NumStates() * excited->NumStates());
    double count = 0.;
    unsigned int i;

    const double ValenceEnergy = ValenceEnergies.find(kappa)->second;

    unsigned int k1, k1max;
    unsigned int k2, k2max;

    ConstStateIterator it_n = core->GetConstStateIterator();
    while(!it_n.AtEnd())
    {
        const State& sn = *(it_n.GetState());

        ConstStateIterator it_alpha = excited->GetConstStateIterator();
        while(!it_alpha.AtEnd())
        {
            const State& salpha = *(it_alpha.GetState());

            if(debug)
            {   count += spacing;
                if(count >= 0.02)
                {   *logstream << ".";
                    count -= 0.02;
                }
            }

            k1 = absdiff(sn.L(), salpha.L());
            if(absdiff(sn.TwoJ(), salpha.TwoJ()) > 2 * k1)
                k1 += 2;

            k1max = (sn.TwoJ() + salpha.TwoJ())/2;

            while(k1 <= k1max)
            {
                double C_nalpha = Constant::Electron3j(sn.TwoJ(), salpha.TwoJ(), k1, 1, -1);

                if(C_nalpha)
                {
                    C_nalpha = C_nalpha * (sn.TwoJ() + 1) * (salpha.TwoJ() + 1);

                    for(i=0; i<mmin(sn.Size(), salpha.Size()); i++)
                    {
                        density[i] = sn.f[i] * salpha.f[i] + Constant::AlphaSquared * sn.g[i] * salpha.g[i];
                    }
                    I.FastCoulombIntegrate(density, P_nalpha, k1, mmin(sn.Size(), salpha.Size()));

                    double SMS_nalpha = 0.;
                    if(NuclearInverseMass && (k1 == 1))
                        SMS_nalpha = -SI.IsotopeShiftIntegral(salpha, sn);

                    ConstStateIterator it_beta = excited->GetConstStateIterator();
                    while(!it_beta.AtEnd())
                    {
                        const State& sbeta = *(it_beta.GetState());

                        double C_abeta;
                        if((sa.L() + sbeta.L() + k1)%2 == 0)
                            C_abeta = Constant::Electron3j(sa.TwoJ(), sbeta.TwoJ(), k1, 1, -1);
                        else
                            C_abeta = 0.;

                        if(C_abeta && ((sa.L() + salpha.L())%2 == (sn.L() + sbeta.L())%2))
                        {
                            C_abeta = C_abeta * (sbeta.TwoJ() + 1);
                            C_abeta = C_abeta/(ValenceEnergy + sn.Energy() - sbeta.Energy() - salpha.Energy() + delta);

                            k2 = absdiff(sn.L(), sbeta.L());
                            if(absdiff(sn.TwoJ(), sbeta.TwoJ()) > 2 * k2)
                                k2 += 2;
                            k2max = (sn.TwoJ() + sbeta.TwoJ())/2;

                            // Sign
                            if((k1 + k2)%2)
                                C_abeta = -C_abeta;

                            while(k2 <= k2max)
                            {
                                double coeff
                                    = C_abeta * C_nalpha * Constant::Electron3j(sa.TwoJ(), salpha.TwoJ(), k2, 1, -1)
                                    * Constant::Electron3j(sbeta.TwoJ(), sn.TwoJ(), k2, 1, -1)
                                    * Constant::Wigner6j(sa.J(), sbeta.J(), k1, sn.J(), salpha.J(), k2);
                                    // Note: The 6j symbol is given incorrectly in Berengut et al. PRA 73, 012504 (2006)

                                if(coeff)
                                {
                                    for(i=0; i<mmin(sn.Size(), sbeta.Size()); i++)
                                    {
                                        density[i] = sn.f[i] * sbeta.f[i] + Constant::AlphaSquared * sn.g[i] * sbeta.g[i];
                                    }
                                    I.FastCoulombIntegrate(density, P_nbeta, k2, mmin(sn.Size(), sbeta.Size()));

                                    double SMS_nbeta = 0.;
                                    if(NuclearInverseMass && (k2 == 1))
                                        SMS_nbeta = -SI.IsotopeShiftIntegral(sbeta, sn);

                                    // R1 = R_k1 (a n, beta alpha)
                                    Y1_size = sbeta.Size();
                                    for(i = 0; i < Y1_size; i++)
                                        Y1[i] = P_nalpha[i] * sbeta.f[i];
        
                                    if(SMS_nalpha)
                                    {
                                        std::vector<double> P1(Y1_size);
                                        SI.IsotopeShiftIntegral(sa.L(), sbeta, &P1);
                                        for(i = 0; i < Y1_size; i++)
                                            Y1[i] = Y1[i] - NuclearInverseMass * SMS_nalpha * P1[i];
                                    }

                                    // R2 = R_k2 (beta alpha, n b) = R_k2 (b n, alpha beta)
                                    Y2_size = salpha.Size();
                                    for(i = 0; i < Y2_size; i++)
                                        Y2[i] = P_nbeta[i] * salpha.f[i];

                                    if(SMS_nbeta)
                                    {   
                                        std::vector<double> P2(Y2_size);
                                        SI.IsotopeShiftIntegral(sa.L(), salpha, &P2);
                                        for(i=0; i<Y2_size; i++)
                                            Y2[i] = Y2[i] - NuclearInverseMass * SMS_nbeta * P2[i];
                                    }

                                    sigma->AddToSigma(Y1, Y2, coeff, Y1_size, Y2_size);
                                }
                                k2 += 2;
                            }
                        }
                        it_beta.Next();
                    } 
                } // C_nalpha
                k1 += 2;
            }
            it_alpha.Next();
        }
        it_n.Next();
    }
}

void MBPTCalculator::CalculateCorrelation4(int kappa, SigmaPotential* sigma) const
{
    const bool debug = DebugOptions.LogMBPT();
    unsigned int MaxStateSize = core->GetConstHFPotential().size();
    const double NuclearInverseMass = core->GetNuclearInverseMass();

    StateInfo sa(100, kappa);

    std::vector<double> density(MaxStateSize);
    std::vector<double> P_nalpha(MaxStateSize);
    std::vector<double> P_malpha(MaxStateSize);
    std::vector<double> Y1(MaxStateSize); unsigned int Y1_size;
    std::vector<double> Y2(MaxStateSize); unsigned int Y2_size;
    CoulombIntegrator I(lattice);
    StateIntegrator SI(lattice);

    if(debug)
        *outstream << "Cor 4:    ";
    double spacing = 1./double(core->NumStates() * excited->NumStates());
    double count = 0.;
    unsigned int i;

    const double ValenceEnergy = ValenceEnergies.find(kappa)->second;

    unsigned int k1, k1max;
    unsigned int k2, k2max;

    ConstStateIterator it_alpha = excited->GetConstStateIterator();
    while(!it_alpha.AtEnd())
    {
        const State& salpha = *(it_alpha.GetState());

        ConstStateIterator it_n = core->GetConstStateIterator();
        while(!it_n.AtEnd())
        {
            const State& sn = *(it_n.GetState());

            if(debug)
            {   count += spacing;
                if(count >= 0.02)
                {   *logstream << ".";
                    count -= 0.02;
                }
            }

            k1 = absdiff(sn.L(), salpha.L());
            if(absdiff(sn.TwoJ(), salpha.TwoJ()) > 2 * k1)
                k1 += 2;

            k1max = (sn.TwoJ() + salpha.TwoJ())/2;

            while(k1 <= k1max)
            {
                double C_nalpha = Constant::Electron3j(sn.TwoJ(), salpha.TwoJ(), k1, 1, -1);

                if(C_nalpha)
                {
                    C_nalpha = C_nalpha * (sn.TwoJ() + 1) * (salpha.TwoJ() + 1);

                    for(i=0; i<mmin(sn.Size(), salpha.Size()); i++)
                    {
                        density[i] = sn.f[i] * salpha.f[i] + Constant::AlphaSquared * sn.g[i] * salpha.g[i];
                    }
                    I.FastCoulombIntegrate(density, P_nalpha, k1, mmin(sn.Size(), salpha.Size()));

                    double SMS_nalpha = 0.;
                    if(NuclearInverseMass && (k1 == 1))
                        SMS_nalpha = -SI.IsotopeShiftIntegral(salpha, sn);

                    ConstStateIterator it_m = core->GetConstStateIterator();
                    while(!it_m.AtEnd())
                    {
                        const State& sm = *(it_m.GetState());

                        double C_am;
                        if((sa.L() + sm.L() + k1)%2 == 0)
                            C_am = Constant::Electron3j(sa.TwoJ(), sm.TwoJ(), k1, 1, -1);
                        else
                            C_am = 0.;

                        if(C_am && ((sa.L() + sn.L())%2 == (sm.L() + salpha.L())%2))
                        {
                            C_am = C_am * (sm.TwoJ() + 1);

                            C_am = C_am/(ValenceEnergy + salpha.Energy() - sn.Energy() - sm.Energy() - delta);

                            k2 = absdiff(sm.L(), salpha.L());
                            if(absdiff(sm.TwoJ(), salpha.TwoJ()) > 2 * k2)
                                k2 += 2;

                            k2max = (sm.TwoJ() + salpha.TwoJ())/2;

                            // Sign
                            if((k1 + k2)%2)
                                C_am = -C_am;

                            while(k2 <= k2max)
                            {
                                double coeff
                                    = C_am * C_nalpha * Constant::Electron3j(sa.TwoJ(), sn.TwoJ(), k2, 1, -1)
                                    * Constant::Electron3j(sm.TwoJ(), salpha.TwoJ(), k2, 1, -1)
                                    * Constant::Wigner6j(sa.J(), sm.J(), k1, salpha.J(), sn.J(), k2);

                                if(coeff)
                                {
                                    for(i=0; i<mmin(sm.Size(), salpha.Size()); i++)
                                    {
                                        density[i] = sm.f[i] * salpha.f[i] + Constant::AlphaSquared * sm.g[i] * salpha.g[i];
                                    }
                                    I.FastCoulombIntegrate(density, P_malpha, k2, mmin(sm.Size(), salpha.Size()));

                                    double SMS_malpha = 0.;
                                    if(NuclearInverseMass && (k2 == 1))
                                        SMS_malpha = -SI.IsotopeShiftIntegral(salpha, sm);

                                    // R1 = R_k1 (a alpha, m n)
                                    Y1_size = sm.Size();
                                    for(i=0; i < Y1_size; i++)
                                        Y1[i] = P_nalpha[i] * sm.f[i];

                                    if(SMS_nalpha)
                                    {
                                        std::vector<double> P1(Y1_size);
                                        SI.IsotopeShiftIntegral(sa.L(), sm, &P1);
                                        for(i=0; i<Y1_size; i++)
                                            Y1[i] = Y1[i] + NuclearInverseMass * SMS_nalpha * P1[i];
                                    }

                                    // R2 = R_k2 (m n, alpha b) = R_k2 (b alpha, n m)
                                    Y2_size = sn.Size();
                                    for(i=0; i < Y2_size; i++)
                                        Y2[i] = P_malpha[i] * sn.f[i];

                                    if(SMS_malpha)
                                    {
                                        std::vector<double> P2(Y2_size);
                                        SI.IsotopeShiftIntegral(sa.L(), sn, &P2);
                                        for(i=0; i < Y2_size; i++)
                                            Y2[i] = Y2[i] + NuclearInverseMass * SMS_malpha * P2[i];
                                    }

                                    sigma->AddToSigma(Y1, Y2, coeff, Y1_size, Y2_size);
                                }
                                k2 += 2;
                            }
                        }
                        it_m.Next();
                    } 
                } // C_nalpha
                k1 += 2;
            }
            it_n.Next();
        }
        it_alpha.Next();
    }
}

double MBPTCalculator::CalculateCorrelation1and3(const StateInfo& sa, const StateInfo& sb) const
{
    const bool debug = DebugOptions.LogMBPT();

    if(debug)
        *outstream << "Cor 1+3:  ";

    const double ValenceEnergy = ValenceEnergies.find(sa.Kappa())->second;

    unsigned int k1, k1max;
    double energy1 = 0., energy3 = 0.;

    // Firstly, get the loop 24
    ConstStateIterator it_n = core->GetConstStateIterator();
    while(!it_n.AtEnd())
    {
        const StateInfo sn = it_n.GetStateInfo();
        const double En = it_n.GetState()->Energy();

        ConstStateIterator it_alpha = excited->GetConstStateIterator();
        while(!it_alpha.AtEnd())
        {
            const StateInfo salpha = it_alpha.GetStateInfo();
            const double Ealpha = it_alpha.GetState()->Energy();

            k1 = kmin(sn, salpha);
            k1max = kmax(sn, salpha);

            while(k1 <= k1max)
            {
                double C_nalpha = Constant::Electron3j(sn.TwoJ(), salpha.TwoJ(), k1, 1, -1);

                if(C_nalpha)
                {
                    C_nalpha = C_nalpha * C_nalpha * sn.MaxNumElectrons() * salpha.MaxNumElectrons()
                                                / (2. * k1 + 1.);
                    C_nalpha = C_nalpha * it_alpha.Weight();

                    // Correlation 1 has excited state beta
                    ConstStateIterator it_beta = excited->GetConstStateIterator();
                    while(!it_beta.AtEnd())
                    {
                        const StateInfo sbeta = it_beta.GetStateInfo();
                        const double Ebeta = it_beta.GetState()->Energy();

                        double coeff;
                        if((sa.L() + sbeta.L() + k1)%2 == 0)
                            coeff = Constant::Electron3j(sa.TwoJ(), sbeta.TwoJ(), k1, 1, -1);
                        else
                            coeff = 0.;

                        if(coeff)
                        {
                            coeff = coeff * coeff * C_nalpha * sbeta.MaxNumElectrons();
                            coeff = coeff/(ValenceEnergy + En - Ebeta - Ealpha + delta);

                            // R1 = R_k1 (a n, beta alpha)
                            // R2 = R_k1 (b n, beta alpha)
                            double R1 = integrals->GetTwoElectronIntegral(k1, sa, sn, sbeta, salpha);
                            double R2 = integrals->GetTwoElectronIntegral(k1, sb, sn, sbeta, salpha);
                            
                            energy1 += R1 * R2 * coeff;
                        }
                        it_beta.Next();
                    }

                    // Correlation 3 has core state m
                    ConstStateIterator it_m = core->GetConstStateIterator();
                    while(!it_m.AtEnd())
                    {
                        const StateInfo sm = it_m.GetStateInfo();
                        const double Em = it_m.GetState()->Energy();

                        double coeff;
                        if((sa.L() + sm.L() + k1)%2 == 0)
                            coeff =  Constant::Electron3j(sa.TwoJ(), sm.TwoJ(), k1, 1, -1);
                        else
                            coeff = 0.;

                        if(coeff)
                        {
                            coeff = coeff * coeff * C_nalpha * sm.MaxNumElectrons();
                            coeff = coeff/(ValenceEnergy + Ealpha - En - Em - delta);

                            // R1 = R_k1 (a alpha, m n)
                            // R2 = R_k1 (b alpha, m n)
                            double R1 = integrals->GetTwoElectronIntegral(k1, sa, salpha, sm, sn);
                            double R2 = integrals->GetTwoElectronIntegral(k1, sb, salpha, sm, sn);

                            energy3 += R1 * R2 * coeff;
                        }
                        it_m.Next();
                    }
                }
                k1 += 2;
            }
            it_alpha.Next();
        }
        it_n.Next();
    }
    
    if(debug)
        *outstream << "  " << energy1 * Constant::HartreeEnergy_cm
                   << "  " << energy3 * Constant::HartreeEnergy_cm << std::endl;
    return energy1 + energy3;
}

double MBPTCalculator::CalculateCorrelation2(const StateInfo& sa, const StateInfo& sb) const
{
    const bool debug = DebugOptions.LogMBPT();

    if(debug)
        *outstream << "Cor 2:    ";

    const double ValenceEnergy = ValenceEnergies.find(sa.Kappa())->second;

    unsigned int k1, k1max;
    unsigned int k2, k2max;
    
    double energy = 0.;

    ConstStateIterator it_n = core->GetConstStateIterator();
    while(!it_n.AtEnd())
    {
        const StateInfo sn = it_n.GetStateInfo();
        const double En = it_n.GetState()->Energy();

        ConstStateIterator it_alpha = excited->GetConstStateIterator();
        while(!it_alpha.AtEnd())
        {
            const StateInfo salpha = it_alpha.GetStateInfo();
            const double Ealpha = it_alpha.GetState()->Energy();

            k1 = kmin(sn, salpha);
            k1max = kmax(sn, salpha);

            while(k1 <= k1max)
            {
                double C_nalpha = Constant::Electron3j(sn.TwoJ(), salpha.TwoJ(), k1, 1, -1);

                if(C_nalpha)
                {
                    C_nalpha = C_nalpha * sn.MaxNumElectrons() * salpha.MaxNumElectrons();

                    ConstStateIterator it_beta = excited->GetConstStateIterator();
                    while(!it_beta.AtEnd())
                    {
                        const StateInfo sbeta = it_beta.GetStateInfo();
                        const double Ebeta = it_beta.GetState()->Energy();

                        double C_abeta;
                        if((sa.L() + sbeta.L() + k1)%2 == 0)
                            C_abeta = Constant::Electron3j(sa.TwoJ(), sbeta.TwoJ(), k1, 1, -1);
                        else
                            C_abeta = 0.;

                        if(C_abeta && ((sa.L() + salpha.L())%2 == (sn.L() + sbeta.L())%2))
                        {
                            C_abeta = C_abeta * sbeta.MaxNumElectrons();
                            C_abeta = C_abeta/(ValenceEnergy + En - Ebeta - Ealpha + delta);

                            k2 = kmin(sn, sbeta, sa, salpha);
                            k2max = kmax(sn, sbeta, sa, salpha);

                            // Sign
                            if((k1 + k2)%2)
                                C_abeta = -C_abeta;

                            while(k2 <= k2max)
                            {
                                double coeff
                                    = C_abeta * C_nalpha * Constant::Electron3j(sa.TwoJ(), salpha.TwoJ(), k2, 1, -1)
                                    * Constant::Electron3j(sbeta.TwoJ(), sn.TwoJ(), k2, 1, -1)
                                    * Constant::Wigner6j(sa.J(), sbeta.J(), k1, sn.J(), salpha.J(), k2);
                                    // Note: The 6j symbol is given incorrectly in Berengut et al. PRA 73, 012504 (2006)

                                if(coeff)
                                {
                                    // R1 = R_k1 (a n, beta alpha)
                                    double R1 = integrals->GetTwoElectronIntegral(k1, sa, sn, sbeta, salpha);

                                    // R2 = R_k2 (beta alpha, n b) = R_k2 (n b, beta alpha)
                                    double R2 = integrals->GetTwoElectronIntegral(k2, sn, sb, sbeta, salpha);
                                    
                                    energy += R1 * R2 * coeff;
                                }
                                k2 += 2;
                            }
                        }
                        it_beta.Next();
                    } 
                } // C_nalpha
                k1 += 2;
            }
            it_alpha.Next();
        }
        it_n.Next();
    }
    
    if(debug)
        *outstream << "  " << energy * Constant::HartreeEnergy_cm << std::endl;
    return energy;    
}

double MBPTCalculator::CalculateCorrelation4(const StateInfo& sa, const StateInfo& sb) const
{
    const bool debug = DebugOptions.LogMBPT();

    if(debug)
        *outstream << "Cor 4:    ";

    const double ValenceEnergy = ValenceEnergies.find(sa.Kappa())->second;

    unsigned int k1, k1max;
    unsigned int k2, k2max;

    double energy = 0.;

    ConstStateIterator it_n = core->GetConstStateIterator();
    while(!it_n.AtEnd())
    {
        const StateInfo sn = it_n.GetStateInfo();
        const double En = it_n.GetState()->Energy();

        ConstStateIterator it_alpha = excited->GetConstStateIterator();
        while(!it_alpha.AtEnd())
        {
            const StateInfo salpha = it_alpha.GetStateInfo();
            const double Ealpha = it_alpha.GetState()->Energy();

            k1 = kmin(sn, salpha);
            k1max = kmax(sn, salpha);

            while(k1 <= k1max)
            {
                double C_nalpha = Constant::Electron3j(sn.TwoJ(), salpha.TwoJ(), k1, 1, -1);

                if(C_nalpha)
                {
                    C_nalpha = C_nalpha * sn.MaxNumElectrons() * salpha.MaxNumElectrons();

                    ConstStateIterator it_m = core->GetConstStateIterator();
                    while(!it_m.AtEnd())
                    {
                        const StateInfo sm = it_m.GetStateInfo();
                        const double Em = it_m.GetState()->Energy();

                        double C_am;
                        if((sa.L() + sm.L() + k1)%2 == 0)
                            C_am = Constant::Electron3j(sa.TwoJ(), sm.TwoJ(), k1, 1, -1);
                        else
                            C_am = 0.;

                        if(C_am && ((sa.L() + sn.L())%2 == (sm.L() + salpha.L())%2))
                        {
                            C_am = C_am * sm.MaxNumElectrons();
                            C_am = C_am/(ValenceEnergy + Ealpha - En - Em - delta);

                            k2 = kmin(sm, salpha, sa, sn);
                            k2max = kmax(sm, salpha, sa, sn);

                            // Sign
                            if((k1 + k2)%2)
                                C_am = -C_am;

                            while(k2 <= k2max)
                            {
                                double coeff
                                    = C_am * C_nalpha * Constant::Electron3j(sa.TwoJ(), sn.TwoJ(), k2, 1, -1)
                                    * Constant::Electron3j(sm.TwoJ(), salpha.TwoJ(), k2, 1, -1)
                                    * Constant::Wigner6j(sa.J(), sm.J(), k1, salpha.J(), sn.J(), k2);

                                if(coeff)
                                {   // R1 = R_k1 (a alpha, m n)
                                    double R1 = integrals->GetTwoElectronIntegral(k1, sa, salpha, sm, sn);

                                    // R2 = R_k2 (m n, alpha b)
                                    double R2 = integrals->GetTwoElectronIntegral(k2, sm, sn, salpha, sb);
                                    
                                    energy += R1 * R2 * coeff;
                                }
                                k2 += 2;
                            }
                        }
                        it_m.Next();
                    } 
                } // C_nalpha
                k1 += 2;
            }
            it_alpha.Next();
        }
        it_n.Next();
    }

    if(debug)
        *outstream << "  " << energy * Constant::HartreeEnergy_cm << std::endl;
    return energy;    
}

double MBPTCalculator::CalculateSubtraction1(const StateInfo& sa, const StateInfo& sb) const
{
    const bool debug = DebugOptions.LogMBPT();

    if(debug)
        *outstream << "Sub 1:    ";

    double energy = 0.;

    ConstStateIterator it_n = core->GetConstStateIterator();
    while(!it_n.AtEnd())
    {
        const StateInfo sn = it_n.GetStateInfo();
        const double En = it_n.GetState()->Energy();

        ConstStateIterator it_alpha = excited->GetConstStateIterator();
        while(!it_alpha.AtEnd())
        {
            const StateInfo salpha = it_alpha.GetStateInfo();
            const double Ealpha = it_alpha.GetState()->Energy();

            if(sn.Kappa() == salpha.Kappa())
            {
                double coeff = integrals->GetOneElectronIntegral(sn, salpha);
                coeff = coeff * sn.MaxNumElectrons();

                coeff = coeff/(En - Ealpha + delta);

                // R1 = R_0 (a n, b alpha)
                double R1 = integrals->GetTwoElectronIntegral(0, sa, sn, sb, salpha);

                // Factor of 2 from identical mirror diagram:
                //   R_0 (a n, b alpha) = R_0 (a alpha, b n)
                // and no SMS since k == 0.
                energy += 2. * R1 * coeff;
            }
            it_alpha.Next();
        }
        it_n.Next();
    }

    if(debug)
        *outstream << "  " << energy * Constant::HartreeEnergy_cm << std::endl;
    return energy;
}

double MBPTCalculator::CalculateSubtraction2(const StateInfo& sa, const StateInfo& sb) const
{
    const bool debug = DebugOptions.LogMBPT();
    if(debug)
        *outstream << "Sub 2:    ";

    double energy = 0.;
    
    unsigned int k1, k1max;

    ConstStateIterator it_n = core->GetConstStateIterator();
    while(!it_n.AtEnd())
    {
        const StateInfo sn = it_n.GetStateInfo();
        const double En = it_n.GetState()->Energy();

        ConstStateIterator it_alpha = excited->GetConstStateIterator();
        while(!it_alpha.AtEnd())
        {
            const StateInfo salpha = it_alpha.GetStateInfo();
            const double Ealpha = it_alpha.GetState()->Energy();

            if(sn.Kappa() == salpha.Kappa())
            {
                double C_nalpha = integrals->GetOneElectronIntegral(sn, salpha);
                C_nalpha = C_nalpha * sn.MaxNumElectrons();
                C_nalpha = C_nalpha/(En - Ealpha + delta);

                k1 = kmin(sa, sn);
                k1max = kmax(sa, sn);

                while(k1 <= k1max)
                {
                    double coeff = Constant::Electron3j(sa.TwoJ(), sn.TwoJ(), k1, 1, -1);
                    coeff = - coeff * coeff * C_nalpha;

                    if(coeff)
                    {
                        // R1 = R_k1 (a alpha, n b), R2 = R_k1 (a n, alpha b)
                        double R1 = integrals->GetTwoElectronIntegral(k1, sa, salpha, sn, sb);
                        double R2 = integrals->GetTwoElectronIntegral(k1, sa, sn, salpha, sb);

                        energy += (R1 + R2) * coeff;
                    }
                    k1 += 2;
                }
            }
            it_alpha.Next();
        }
        it_n.Next();
    }

    if(debug)
        *outstream << "  " << energy * Constant::HartreeEnergy_cm << std::endl;
    return energy;
}

double MBPTCalculator::CalculateSubtraction3(const StateInfo& sa, const StateInfo& sb) const
{
    const bool debug = DebugOptions.LogMBPT();

    if(debug)
        *outstream << "Sub 3:    ";

    double energy = 0.;
    double ValenceEnergy = ValenceEnergies.find(sa.Kappa())->second;

    ConstStateIterator it_n = core->GetConstStateIterator();
    while(!it_n.AtEnd())
    {
        const StateInfo sn = it_n.GetStateInfo();
        const double En = it_n.GetState()->Energy();

        if(sn.Kappa() == sa.Kappa())
        {
            double term = integrals->GetOneElectronIntegral(sa, sn) * integrals->GetOneElectronIntegral(sn, sb);
            term = term/(En - ValenceEnergy + delta);

            energy = energy - term;
        }
        it_n.Next();
    }

    if(debug)
        *outstream << "  " << energy * Constant::HartreeEnergy_cm << std::endl;
    return energy;
}

double MBPTCalculator::CalculateTwoElectron1(unsigned int k, const StateInfo& sa, const StateInfo& sb, const StateInfo& sc, const StateInfo& sd) const
{
    const bool debug = DebugOptions.LogMBPT();

    if(debug)
        *outstream << "TwoE 1:   ";

    double energy = 0.;

    ConstStateIterator it_n = core->GetConstStateIterator();
    while(!it_n.AtEnd())
    {
        const StateInfo sn = it_n.GetStateInfo();
        const double En = it_n.GetState()->Energy();

        ConstStateIterator it_alpha = excited->GetConstStateIterator();
        while(!it_alpha.AtEnd())
        {
            const StateInfo salpha = it_alpha.GetStateInfo();
            const double Ealpha = it_alpha.GetState()->Energy();

            double coeff;
            if((sn.L() + salpha.L() + k)%2)
                coeff = 0.;
            else
                coeff = Constant::Electron3j(sn.TwoJ(), salpha.TwoJ(), k, 1, -1);

            if(coeff)
            {
                coeff = coeff * coeff * sn.MaxNumElectrons() * salpha.MaxNumElectrons()
                                        / (2. * k + 1.);
                coeff = coeff/(En - Ealpha + delta);

                // There are two diagrams:
                //  1. R_k(a n, c alpha) * R_k(alpha b, n d)
                //  2. R_k(a alpha, c n) * R_k(n b, alpha d)
                // The first order SMS cancels for these two diagrams.
                double R1 = integrals->GetTwoElectronIntegral(k, sa, sn, sc, salpha)
                            * integrals->GetTwoElectronIntegral(k, salpha, sb, sn, sd);
                double R2 = integrals->GetTwoElectronIntegral(k, sa, salpha, sc, sn)
                            * integrals->GetTwoElectronIntegral(k, sn, sb, salpha, sd);

                energy += (R1 + R2) * coeff;
            }
            it_alpha.Next();
        }
        it_n.Next();
    }

    if(debug)
        *outstream << "  " << energy * Constant::HartreeEnergy_cm << std::endl;
    return energy;
}

double MBPTCalculator::CalculateTwoElectron2(unsigned int k, const StateInfo& sa, const StateInfo& sb, const StateInfo& sc, const StateInfo& sd) const
{
    const bool debug = DebugOptions.LogMBPT();

    if(debug)
        *outstream << "TwoE 2/3: ";

    double energy = 0.;
    const double coeff_ac = Constant::Electron3j(sa.TwoJ(), sc.TwoJ(), k, 1, -1);
    const double coeff_bd = Constant::Electron3j(sb.TwoJ(), sd.TwoJ(), k, 1, -1);
    if(!coeff_ac || !coeff_bd)
        return energy;

    unsigned int k1, k1max;

    ConstStateIterator it_n = core->GetConstStateIterator();
    while(!it_n.AtEnd())
    {
        const StateInfo sn = it_n.GetStateInfo();
        const double En = it_n.GetState()->Energy();

        ConstStateIterator it_alpha = excited->GetConstStateIterator();
        while(!it_alpha.AtEnd())
        {
            const StateInfo salpha = it_alpha.GetStateInfo();
            const double Ealpha = it_alpha.GetState()->Energy();

            double C_nalpha;
            if((sn.L() + salpha.L() + k)%2)
                C_nalpha = 0.;
            else
                C_nalpha = Constant::Electron3j(sn.TwoJ(), salpha.TwoJ(), k, 1, -1);                

            if(C_nalpha)
            {
                C_nalpha = C_nalpha * sn.MaxNumElectrons() * salpha.MaxNumElectrons();
                C_nalpha = C_nalpha/(En - Ealpha + delta);

                k1 = kmin(sa, sn, salpha, sc);
                k1max = kmax(sa, sn, salpha, sc);

                while(k1 <= k1max)
                {
                    double coeff = Constant::Electron3j(sa.TwoJ(), sn.TwoJ(), k1, 1, -1) *
                                   Constant::Electron3j(salpha.TwoJ(), sc.TwoJ(), k1, 1, -1) *
                                   Constant::Wigner6j(sa.J(), sc.J(), double(k), salpha.J(), sn.J(), double(k1)) *
                                   C_nalpha / coeff_ac;
                    if((k1 + k)%2)
                        coeff = -coeff;
    // REMOVE:
                    if((salpha.L() + sc.L() + k1)%2)
                    {   *errstream << "TwoElectron2: messed up" << std::endl;
                        coeff = 0;
                    }

                    if(coeff)
                    {
                        // R1 = R_k1 (a alpha, n c)
                        double R1 = integrals->GetTwoElectronIntegral(k1, sa, salpha, sn, sc);

                        // R2 = R_k (n b, alpha d)
                        double R2 = integrals->GetTwoElectronIntegral(k, sn, sb, salpha, sd);

                        energy += R1 * R2 * coeff;
                    }
                    k1 += 2;
                }

                // Mirror diagram
                k1 = kmin(sb, sn, salpha, sd);
                k1max = kmax(sb, sn, salpha, sd);

                while(k1 <= k1max)
                {
                    double coeff = Constant::Electron3j(sb.TwoJ(), sn.TwoJ(), k1, 1, -1) *
                                   Constant::Electron3j(salpha.TwoJ(), sd.TwoJ(), k1, 1, -1) *
                                   Constant::Wigner6j(sb.J(), sd.J(), double(k), salpha.J(), sn.J(), double(k1)) *
                                   C_nalpha / coeff_bd;
                    if((k1 + k)%2)
                        coeff = -coeff;
    // REMOVE:
                    if((salpha.L() + sd.L() + k1)%2)
                    {   *errstream << "TwoElectron2: messed up mirror" << std::endl;
                        coeff = 0;
                    }

                    if(coeff)
                    {
                        // R1 = R_k1 (b alpha, n d)
                        double R1 = integrals->GetTwoElectronIntegral(k1, sb, salpha, sn, sd);
                        
                        // R2 = R_k (a n, c alpha)
                        double R2 = integrals->GetTwoElectronIntegral(k, sa, sn, sc, salpha);

                        energy += R1 * R2 * coeff;
                    }
                    k1 += 2;
                }
            }
            it_alpha.Next();
        }
        it_n.Next();
    }

    if(debug)
        *outstream << "  " << energy * Constant::HartreeEnergy_cm << std::endl;
    return energy;
}

double MBPTCalculator::CalculateTwoElectron3(unsigned int k, const StateInfo& sa, const StateInfo& sb, const StateInfo& sc, const StateInfo& sd) const
{
    return CalculateTwoElectron2(k, sc, sd, sa, sb);
}

double MBPTCalculator::CalculateTwoElectron4(unsigned int k, const StateInfo& sa, const StateInfo& sb, const StateInfo& sc, const StateInfo& sd) const
{
    const bool debug = DebugOptions.LogMBPT();

    if(debug)
        *outstream << "TwoE 4/5: ";

    double energy = 0.;
    const double coeff_ac = Constant::Electron3j(sa.TwoJ(), sc.TwoJ(), k, 1, -1);
    const double coeff_bd = Constant::Electron3j(sb.TwoJ(), sd.TwoJ(), k, 1, -1);
    if(!coeff_ac || !coeff_bd)
        return energy;

    unsigned int k1, k1max;
    unsigned int k2, k2max;

    ConstStateIterator it_n = core->GetConstStateIterator();
    while(!it_n.AtEnd())
    {
        const StateInfo sn = it_n.GetStateInfo();
        const double En = it_n.GetState()->Energy();

        ConstStateIterator it_alpha = excited->GetConstStateIterator();
        while(!it_alpha.AtEnd())
        {
            const StateInfo salpha = it_alpha.GetStateInfo();
            const double Ealpha = it_alpha.GetState()->Energy();

            double C_nalpha;
            if((sa.L() + sn.L())%2 != (salpha.L() + sd.L())%2)
                C_nalpha = 0.;
            else if((sn.L() + sc.L())%2 != (sb.L() + salpha.L())%2)
                C_nalpha = 0.;
            else
                C_nalpha = double(sn.MaxNumElectrons()) * double(salpha.MaxNumElectrons()) * (2. * double(k) + 1.);

            if(C_nalpha)
            {
                C_nalpha = C_nalpha/(coeff_ac*coeff_bd);
                C_nalpha = C_nalpha/(En - Ealpha + delta);
                unsigned int phase = (unsigned int)(sa.TwoJ() + sb.TwoJ() + sc.TwoJ() + sd.TwoJ() + sn.TwoJ() + salpha.TwoJ())/2;
                if(phase%2)
                    C_nalpha = -C_nalpha;

                k1 = kmin(sa, sn, salpha, sd);
                k1max = kmax(sa, sn, salpha, sd);

                while(k1 <= k1max)
                {
                    double coeff_ad = Constant::Electron3j(sa.TwoJ(), sn.TwoJ(), k1, 1, -1) *
                                      Constant::Electron3j(salpha.TwoJ(), sd.TwoJ(), k1, 1, -1);

                    if(coeff_ad)
                    {
                        // R1 = R_k1 (a alpha, n d)
                        double R1 = integrals->GetTwoElectronIntegral(k1, sa, salpha, sn, sd);

                        k2 = kmin(sn, sc, sb, salpha);
                        k2max = kmax(sn, sc, sb, salpha);

                        while(k2 <= k2max)
                        {
                            double coeff = Constant::Electron3j(sn.TwoJ(), sc.TwoJ(), k2, 1, -1) *
                                           Constant::Electron3j(sb.TwoJ(), salpha.TwoJ(), k2, 1, -1);
                            if(coeff)
                                coeff = coeff * Constant::Wigner6j(sc.J(), sa.J(), k, k1, k2, sn.J())
                                              * Constant::Wigner6j(sb.J(), sd.J(), k, k1, k2, salpha.J());

                            if(coeff)
                            {   
                                coeff = coeff * coeff_ad * C_nalpha;

                                // R2 = R_k2 (2b, c4)
                                double R2 = integrals->GetTwoElectronIntegral(k2, sn, sb, sc, salpha);

                                energy += R1 * R2 * coeff;
                            }
                            k2 += 2;
                        }
                    }
                    k1 += 2;
                }
            }
            it_alpha.Next();
        }
        it_n.Next();
    }

    if(debug)
        *outstream << "  " << energy * Constant::HartreeEnergy_cm << std::endl;
    return energy;
}

double MBPTCalculator::CalculateTwoElectron5(unsigned int k, const StateInfo& sa, const StateInfo& sb, const StateInfo& sc, const StateInfo& sd) const
{
    return CalculateTwoElectron4(k, sb, sa, sd, sc);
}

double MBPTCalculator::CalculateTwoElectron6(unsigned int k, const StateInfo& sa, const StateInfo& sb, const StateInfo& sc, const StateInfo& sd) const
{
    const bool debug = DebugOptions.LogMBPT();

    if(debug)
        *outstream << "TwoE 6:   ";

    double energy = 0.;
    const double coeff_ac = Constant::Electron3j(sa.TwoJ(), sc.TwoJ(), k, 1, -1);
    const double coeff_bd = Constant::Electron3j(sb.TwoJ(), sd.TwoJ(), k, 1, -1);
    if(!coeff_ac || !coeff_bd)
        return energy;

    const double ValenceEnergy = ValenceEnergies.find(sa.Kappa())->second + ValenceEnergies.find(sb.Kappa())->second;

    unsigned int k1, k1max;
    unsigned int k2, k2max;

    ConstStateIterator it_m = core->GetConstStateIterator();
    while(!it_m.AtEnd())
    {
        const StateInfo sm = it_m.GetStateInfo();
        const double Em = it_m.GetState()->Energy();

        ConstStateIterator it_n = core->GetConstStateIterator();
        while(!it_n.AtEnd())
        {
            const StateInfo sn = it_n.GetStateInfo();
            const double En = it_n.GetState()->Energy();

            double coeff_mn;
            if((sa.L() + sm.L())%2 != (sb.L() + sn.L())%2)
                coeff_mn = 0.;
            else if((sm.L() + sc.L())%2 != (sd.L() + sn.L())%2)
                coeff_mn = 0.;
            else
                coeff_mn = double(sm.MaxNumElectrons()) * double(sn.MaxNumElectrons()) * (2. * double(k) + 1.);

            if(coeff_mn)
            {
                coeff_mn = coeff_mn/(coeff_ac*coeff_bd);
                coeff_mn = coeff_mn/(Em + En - ValenceEnergy + delta);
                unsigned int phase = (unsigned int)(sa.TwoJ() + sb.TwoJ() + sc.TwoJ() + sd.TwoJ() + sm.TwoJ() + sn.TwoJ())/2;
                if((phase + k + 1)%2)
                    coeff_mn = -coeff_mn;

                k1 = kmin(sa, sm, sb, sn);
                k1max = kmax(sa, sm, sb, sn);

                while(k1 <= k1max)
                {
                    double coeff_ab = Constant::Electron3j(sa.TwoJ(), sm.TwoJ(), k1, 1, -1) *
                                      Constant::Electron3j(sb.TwoJ(), sn.TwoJ(), k1, 1, -1);

                    if(coeff_ab)
                    {
                        // R1 = R_k1 (ab, mn)
                        double R1 = integrals->GetTwoElectronIntegral(k1, sa, sb, sm, sn);
                        
                        k2 = kmin(sm, sc, sn, sd);
                        k2max = kmax(sm, sc, sn, sd);

                        while(k2 <= k2max)
                        {
                            double coeff = Constant::Electron3j(sm.TwoJ(), sc.TwoJ(), k2, 1, -1) *
                                           Constant::Electron3j(sn.TwoJ(), sd.TwoJ(), k2, 1, -1);
                            if(coeff)
                                coeff = coeff * Constant::Wigner6j(sc.J(), sa.J(), k, k1, k2, sm.J())
                                              * Constant::Wigner6j(sd.J(), sb.J(), k, k1, k2, sn.J());

                            if(coeff)
                            {   
                                coeff = coeff * coeff_ab * coeff_mn;
                                if((k1 + k2)%2)
                                    coeff = -coeff;

                                // R2 = R_k2 (mn, cd)
                                double R2 = integrals->GetTwoElectronIntegral(k2, sm, sn, sc, sd);

                                energy += R1 * R2 * coeff;
                            }
                            k2 += 2;
                        }
                    }
                    k1 += 2;
                }
            }
            it_n.Next();
        }
        it_m.Next();
    }

    if(debug)
        *outstream << "  " << energy * Constant::HartreeEnergy_cm << std::endl;
    return energy;
}

double MBPTCalculator::CalculateTwoElectronSub(unsigned int k, const StateInfo& sa, const StateInfo& sb, const StateInfo& sc, const StateInfo& sd) const
{
    const bool debug = DebugOptions.LogMBPT();

    if(debug)
        *outstream << "2eSub :   ";

    double energy = 0.;

    const double Ea = ValenceEnergies.find(sa.Kappa())->second;
    const double Eb = ValenceEnergies.find(sb.Kappa())->second;
    const double Ec = ValenceEnergies.find(sc.Kappa())->second;
    const double Ed = ValenceEnergies.find(sd.Kappa())->second;

    // Hole line is attached to sa or sc
    ConstStateIterator it_n = core->GetConstStateIterator();
    while(!it_n.AtEnd())
    {
        const StateInfo sn = it_n.GetStateInfo();
        const double En = it_n.GetState()->Energy();

        if(sn.Kappa() == sa.Kappa())
        {
            double R1 = integrals->GetTwoElectronIntegral(k, sn, sb, sc, sd);
            energy -= R1 * integrals->GetOneElectronIntegral(sa, sn) / (En - Ea + delta);
        }

        if(sn.Kappa() == sc.Kappa())
        {
            double R1 = integrals->GetTwoElectronIntegral(k, sa, sb, sn, sd);
            energy -= R1 * integrals->GetOneElectronIntegral(sn, sc) / (En - Ec + delta);
        }

        it_n.Next();
    }

    // Hole line is attached to sb or sd.
    it_n.First();
    while(!it_n.AtEnd())
    {
        const StateInfo sn = it_n.GetStateInfo();
        const double En = it_n.GetState()->Energy();

        if(sn.Kappa() == sb.Kappa())
        {
            double R1 = integrals->GetTwoElectronIntegral(k, sa, sn, sc, sd);
            energy -= R1 * integrals->GetOneElectronIntegral(sb, sn) / (En - Eb + delta);
        }

        if(sn.Kappa() == sd.Kappa())
        {
            double R1 = integrals->GetTwoElectronIntegral(k, sa, sb, sc, sn);
            energy -= R1 * integrals->GetOneElectronIntegral(sn, sd) / (En - Ed + delta);
        }

        it_n.Next();
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

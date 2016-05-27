#include "BruecknerSigmaCalculator.h"
#include "Include.h"
#ifdef AMBIT_USE_MPI
#include <mpi.h>
#endif

BruecknerSigmaCalculator::BruecknerSigmaCalculator(pOrbitalManagerConst orbitals, pSpinorOperatorConst one_body, pHartreeY two_body, const std::string& fermi_orbitals):
    MBPTCalculator(orbitals, fermi_orbitals), hf(one_body), hartreeY(two_body), core(orbitals->core), excited(orbitals->excited)
{}

void BruecknerSigmaCalculator::GetSecondOrderSigma(int kappa, SigmaPotential& sigma)
{
    if(DebugOptions.LogMBPT())
        *outstream << "\nkappa = " << kappa << std::endl;

    CalculateCorrelation1and3(kappa, sigma);
    CalculateCorrelation2(kappa, sigma);
    CalculateCorrelation4(kappa, sigma);
}

void BruecknerSigmaCalculator::CalculateCorrelation1and3(int kappa, SigmaPotential& sigma)
{
    const bool debug = DebugOptions.LogMBPT();

    int external_twoJ = 2*abs(kappa) - 1;
    int external_L = (kappa > 0)? kappa: (-kappa-1);

    if(debug)
        *outstream << "Cor 1+3:  ";
    double spacing = 1./double(core->size() * excited->size());
    double count = 0.;

    const double ValenceEnergy = ValenceEnergies.find(kappa)->second;
    MathConstant* constants = MathConstant::Instance()->Instance();
    unsigned int sigma_size = sigma.size();

    unsigned int k1, k1max;

#ifdef AMBIT_USE_MPI
    SigmaPotential new_sigma(sigma);
    new_sigma.clear();
    int proc = 0;  // Processor index
#endif

    // Firstly, get the loop 24
    auto it_n = core->begin();
    while(it_n != core->end())
    {
        const Orbital& sn = *it_n->second;

        auto it_alpha = excited->begin();
        while(it_alpha != excited->end())
        {
            const Orbital& salpha = *(it_alpha->second);

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

            k1max = kmax(it_n->first, it_alpha->first);

            while(k1 <= k1max)
            {
                #ifdef AMBIT_USE_MPI
                if(proc == ProcessorRank)
                {
                #endif

                double C_nalpha = constants->Electron3j(sn.TwoJ(), salpha.TwoJ(), k1);
                hartreeY->SetParameters(k1, it_n->second, it_alpha->second);

                if(C_nalpha && !hartreeY->isZero())
                {
                    C_nalpha = C_nalpha * C_nalpha * it_n->first.MaxNumElectrons() * it_alpha->first.MaxNumElectrons()
                                                / (2. * k1 + 1.);

                    // Correlation 1 has excited state beta
                    auto it_beta = excited->begin();
                    while(it_beta != excited->end())
                    {
                        const Orbital& sbeta = *(it_beta->second);

                        double coeff;
                        if(InQSpace(OrbitalInfo(sn), OrbitalInfo(salpha), OrbitalInfo(sbeta)) && (external_L + sbeta.L() + k1)%2 == 0)
                            coeff = constants->Electron3j(external_twoJ, sbeta.TwoJ(), k1);
                        else
                            coeff = 0.;

                        if(coeff)
                        {
                            coeff = coeff * coeff * C_nalpha * it_beta->first.MaxNumElectrons();
                            coeff = coeff/(ValenceEnergy + sn.Energy() - sbeta.Energy() - salpha.Energy() + delta);

                            // R1 = R_k1 (a n, beta alpha)
                            // R2 = R_k1 (b n, beta alpha)
                            SpinorFunction Ybeta = hartreeY->ApplyTo(sbeta, kappa);
                            Ybeta.resize(sigma_size);

                            #ifdef AMBIT_USE_MPI
                            new_sigma.AddToSigma(Ybeta, Ybeta, coeff);
                            #else
                            sigma.AddToSigma(Ybeta, Ybeta, coeff);
                            #endif
                        }

                        it_beta++;
                    }

                    // Correlation 3 has core state m
                    auto it_m = core->begin();
                    while(it_m != core->end())
                    {
                        const Orbital& sm = *(it_m->second);

                        double coeff;
                        if(InQSpace(OrbitalInfo(sn), OrbitalInfo(salpha), OrbitalInfo(sm)) && (external_L + sm.L() + k1)%2 == 0)
                            coeff =  constants->Electron3j(external_twoJ, sm.TwoJ(), k1);
                        else
                            coeff = 0.;

                        if(coeff)
                        {
                            coeff = coeff * coeff * C_nalpha * it_m->first.MaxNumElectrons();
                            coeff = coeff/(ValenceEnergy + salpha.Energy() - sn.Energy() - sm.Energy() - delta);

                            // R1 = R_k1 (a alpha, m n)
                            // R2 = R_k1 (b alpha, m n)
                            SpinorFunction Ym = hartreeY->ApplyTo(sm, kappa, true);
                            Ym.resize(sigma_size);

                            #ifdef AMBIT_USE_MPI
                            new_sigma.AddToSigma(Ym, Ym, coeff);
                            #else
                            sigma.AddToSigma(Ym, Ym, coeff);
                            #endif
                        }
                        it_m++;
                    }
                } // C_nalpha

                #ifdef AMBIT_USE_MPI
                }

                if(++proc == NumProcessors)
                    proc = 0;
                #endif

                k1 += 2;
            }
            it_alpha++;
        }
        it_n++;
    }

#ifdef AMBIT_USE_MPI
    SigmaMatrix reduced = SigmaMatrix::Zero(sigma.matrix_size, sigma.matrix_size);
    MPI_Allreduce(new_sigma.ff.data(), reduced.data(), sigma.matrix_size * sigma.matrix_size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    sigma.ff += reduced;

    if(sigma.use_fg)
    {   MPI_Allreduce(new_sigma.fg.data(), reduced.data(), sigma.matrix_size * sigma.matrix_size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        sigma.fg += reduced;
        MPI_Allreduce(new_sigma.gf.data(), reduced.data(), sigma.matrix_size * sigma.matrix_size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        sigma.gf += reduced;
    }
    if(sigma.use_gg)
    {   MPI_Allreduce(new_sigma.gg.data(), reduced.data(), sigma.matrix_size * sigma.matrix_size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        sigma.gg += reduced;
    }
#endif
}

void BruecknerSigmaCalculator::CalculateCorrelation2(int kappa, SigmaPotential& sigma)
{
    const bool debug = DebugOptions.LogMBPT();

    int external_twoJ = 2*abs(kappa) - 1;
    int external_L = (kappa > 0)? kappa: (-kappa-1);
    double external_J = double(external_twoJ)/2.;

    if(debug)
        *outstream << "Cor 2:    ";
    double spacing = 1./double(core->size() * excited->size());
    double count = 0.;

    const double ValenceEnergy = ValenceEnergies.find(kappa)->second;
    MathConstant* constants = MathConstant::Instance();
    unsigned int sigma_size = sigma.size();

    unsigned int k1, k1max;
    unsigned int k2, k2max;

    pHartreeY hartreeY1(hartreeY->Clone());
    pHartreeY hartreeY2(hartreeY->Clone());

#ifdef AMBIT_USE_MPI
    SigmaPotential new_sigma(sigma);
    new_sigma.clear();
    int proc = 0;  // Processor index
#endif

    auto it_n = core->begin();
    while(it_n != core->end())
    {
        const Orbital& sn = *(it_n->second);

        auto it_alpha = excited->begin();
        while(it_alpha != excited->end())
        {
            const Orbital& salpha = *(it_alpha->second);

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
                #ifdef AMBIT_USE_MPI
                if(proc == ProcessorRank)
                {
                #endif

                double C_nalpha = constants->Electron3j(sn.TwoJ(), salpha.TwoJ(), k1);
                hartreeY1->SetParameters(k1, it_n->second, it_alpha->second);

                if(C_nalpha && !hartreeY1->isZero())
                {
                    C_nalpha = C_nalpha * it_n->first.MaxNumElectrons() * it_alpha->first.MaxNumElectrons();

                    auto it_beta = excited->begin();
                    while(it_beta != excited->end())
                    {
                        const Orbital& sbeta = *(it_beta->second);

                        double C_abeta;
                        if(InQSpace(OrbitalInfo(sn), OrbitalInfo(salpha), OrbitalInfo(sbeta)) && (external_L + sbeta.L() + k1)%2 == 0)
                            C_abeta = constants->Electron3j(external_twoJ, sbeta.TwoJ(), k1);
                        else
                            C_abeta = 0.;

                        if(C_abeta && (external_L + salpha.L() + sn.L() + sbeta.L())%2 == 0)
                        {
                            C_abeta = C_abeta * it_beta->first.MaxNumElectrons();
                            C_abeta = C_abeta/(ValenceEnergy + sn.Energy() - sbeta.Energy() - salpha.Energy() + delta);

                            // R1 = R_k1 (a n, beta alpha)
                            SpinorFunction Y1 = hartreeY1->ApplyTo(sbeta, kappa);
                            Y1.resize(sigma_size);

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
                                = C_abeta * C_nalpha * constants->Electron3j(external_twoJ, salpha.TwoJ(), k2)
                                * constants->Electron3j(sbeta.TwoJ(), sn.TwoJ(), k2)
                                * constants->Wigner6j(external_J, sbeta.J(), k1, sn.J(), salpha.J(), k2);
                                // Note: The 6j symbol is given incorrectly in Berengut et al. PRA 73, 012504 (2006)

                                if(coeff)
                                {
                                    hartreeY2->SetParameters(k2, it_n->second, it_beta->second);

                                    // R2 = R_k2 (beta alpha, n b) = R_k2 (b n, alpha beta)
                                    SpinorFunction Y2 = hartreeY2->ApplyTo(salpha, kappa);
                                    if(Y2.size())
                                    {
                                        Y2.resize(sigma_size);

                                        #ifdef AMBIT_USE_MPI
                                        new_sigma.AddToSigma(Y1, Y2, coeff);
                                        #else
                                        sigma.AddToSigma(Y1, Y2, coeff);
                                        #endif
                                    }
                                }
                                k2 += 2;
                            }
                        }
                        it_beta++;
                    }
                } // C_nalpha

                #ifdef AMBIT_USE_MPI
                }

                if(++proc == NumProcessors)
                    proc = 0;
                #endif

                k1 += 2;
            }
            it_alpha++;
        }
        it_n++;
    }

#ifdef AMBIT_USE_MPI
    SigmaMatrix reduced = SigmaMatrix::Zero(sigma.matrix_size, sigma.matrix_size);
    MPI_Allreduce(new_sigma.ff.data(), reduced.data(), sigma.matrix_size * sigma.matrix_size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    sigma.ff += reduced;

    if(sigma.use_fg)
    {   MPI_Allreduce(new_sigma.fg.data(), reduced.data(), sigma.matrix_size * sigma.matrix_size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        sigma.fg += reduced;
        MPI_Allreduce(new_sigma.gf.data(), reduced.data(), sigma.matrix_size * sigma.matrix_size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        sigma.gf += reduced;
    }
    if(sigma.use_gg)
    {   MPI_Allreduce(new_sigma.gg.data(), reduced.data(), sigma.matrix_size * sigma.matrix_size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        sigma.gg += reduced;
    }
#endif
}

void BruecknerSigmaCalculator::CalculateCorrelation4(int kappa, SigmaPotential& sigma)
{
    const bool debug = DebugOptions.LogMBPT();

    int external_twoJ = 2*abs(kappa) - 1;
    int external_L = (kappa > 0)? kappa: (-kappa-1);
    double external_J = double(external_twoJ)/2.;

    if(debug)
        *outstream << "Cor 4:    ";
    double spacing = 1./double(core->size() * excited->size());
    double count = 0.;

    const double ValenceEnergy = ValenceEnergies.find(kappa)->second;
    MathConstant* constants = MathConstant::Instance();
    unsigned int sigma_size = sigma.size();

    unsigned int k1, k1max;
    unsigned int k2, k2max;

    pHartreeY hartreeY1(hartreeY->Clone());
    pHartreeY hartreeY2(hartreeY->Clone());

#ifdef AMBIT_USE_MPI
    SigmaPotential new_sigma(sigma);
    new_sigma.clear();
    int proc = 0;  // Processor index
#endif

    auto it_n = core->begin();
    while(it_n != core->end())
    {
        const OrbitalInfo& info_n = it_n->first;
        const Orbital& sn = *(it_n->second);

        auto it_alpha = excited->begin();
        while(it_alpha != excited->end())
        {
            const OrbitalInfo& info_alpha = it_alpha->first;
            const Orbital& salpha = *(it_alpha->second);

            if(debug)
            {   count += spacing;
                if(count >= 0.02)
                {   *logstream << ".";
                    count -= 0.02;
                }
            }

            k1 = kmin(info_n, info_alpha);
            k1max = kmax(info_n, info_alpha);

            while(k1 <= k1max)
            {
                #ifdef AMBIT_USE_MPI
                if(proc == ProcessorRank)
                {
                #endif

                double C_nalpha = constants->Electron3j(sn.TwoJ(), salpha.TwoJ(), k1);
                hartreeY1->SetParameters(k1, it_alpha->second, it_n->second);

                if(C_nalpha && !hartreeY1->isZero())
                {
                    C_nalpha = C_nalpha * info_n.MaxNumElectrons() * info_alpha.MaxNumElectrons();

                    auto it_m = core->begin();
                    while(it_m != core->end())
                    {
                        const OrbitalInfo& info_m = it_m->first;
                        const Orbital& sm = *(it_m->second);

                        double C_am;
                        if(InQSpace(OrbitalInfo(sn), OrbitalInfo(salpha), OrbitalInfo(sm)) && (external_L + sm.L() + k1)%2 == 0)
                            C_am = constants->Electron3j(external_twoJ, sm.TwoJ(), k1);
                        else
                            C_am = 0.;

                        if(C_am && (external_L + sn.L() + sm.L() + salpha.L())%2 == 0)
                        {
                            C_am = C_am * info_m.MaxNumElectrons();
                            C_am = C_am/(ValenceEnergy + salpha.Energy() - sn.Energy() - sm.Energy() - delta);

                            // R1 = R_k1 (a alpha, m n)
                            SpinorFunction Y1 = hartreeY1->ApplyTo(sm, kappa);
                            Y1.resize(sigma_size);

                            // Note: we can rely on HartreeY to check triangle and parity for (k2, n, b)
                            k2 = kmin(info_m, info_alpha);
                            k2max = kmax(info_m, info_alpha);

                            // Sign
                            if((k1 + k2)%2)
                                C_am = -C_am;

                            while(k2 <= k2max)
                            {
                                double coeff
                                = C_am * C_nalpha * constants->Electron3j(external_twoJ, sn.TwoJ(), k2)
                                    * constants->Electron3j(sm.TwoJ(), salpha.TwoJ(), k2)
                                    * constants->Wigner6j(external_J, sm.J(), k1, salpha.J(), sn.J(), k2);

                                if(coeff)
                                {
                                    hartreeY2->SetParameters(k2, it_alpha->second, it_m->second);

                                    // R2 = R_k2 (m n, alpha b) = R_k2 (b alpha, n m)
                                    SpinorFunction Y2 = hartreeY2->ApplyTo(sn, kappa);
                                    if(Y2.size())
                                    {
                                        Y2.resize(sigma_size);

                                        #ifdef AMBIT_USE_MPI
                                        new_sigma.AddToSigma(Y1, Y2, coeff);
                                        #else
                                        sigma.AddToSigma(Y1, Y2, coeff);
                                        #endif
                                    }
                                }
                                k2 += 2;
                            }
                        }
                        it_m++;
                    }
                } // C_nalpha

                #ifdef AMBIT_USE_MPI
                }

                if(++proc == NumProcessors)
                    proc = 0;
                #endif

                k1 += 2;
            }
            it_alpha++;
        }
        it_n++;
    }

#ifdef AMBIT_USE_MPI
    SigmaMatrix reduced = SigmaMatrix::Zero(sigma.matrix_size, sigma.matrix_size);
    MPI_Allreduce(new_sigma.ff.data(), reduced.data(), sigma.matrix_size * sigma.matrix_size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    sigma.ff += reduced;

    if(sigma.use_fg)
    {   MPI_Allreduce(new_sigma.fg.data(), reduced.data(), sigma.matrix_size * sigma.matrix_size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        sigma.fg += reduced;
        MPI_Allreduce(new_sigma.gf.data(), reduced.data(), sigma.matrix_size * sigma.matrix_size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        sigma.gf += reduced;
    }
    if(sigma.use_gg)
    {   MPI_Allreduce(new_sigma.gg.data(), reduced.data(), sigma.matrix_size * sigma.matrix_size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        sigma.gg += reduced;
    }
#endif
}

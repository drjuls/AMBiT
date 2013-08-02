#ifdef _MPI
#include <mpi.h>
#endif
#include "Include.h"
#include "CIIntegralsMBPT.h"
#include "HartreeFock/StateIntegrator.h"
#include "Universal/CoulombIntegrator.h"
#include "Universal/PhysicalConstant.h"
#include <time.h>

inline void swap(unsigned int& i1, unsigned int& i2)
{   unsigned int temp = i1;
    i1 = i2;
    i2 = temp;
}

CIIntegralsMBPT::~CIIntegralsMBPT()
{
    std::map<int, SigmaPotential*>::iterator it = Sigma1.begin();
    while(it != Sigma1.end())
    {   delete it->second;
        it++;
    }
    Sigma1.clear();
}

/** Calculate number of elements that will be stored. */
unsigned int CIIntegralsMBPT::GetStorageSize() const
{
    if(PT)
    {   unsigned int slater_storage = PT->GetStorageSize(&states);
        *outstream << "Num Slater Integrals stored for MBPT: " << slater_storage << std::endl;
    }
    if(ValencePT)
    {   unsigned int slater_storage = ValencePT->GetStorageSize(&states);
        *outstream << "Num Slater Integrals stored for MBPT: " << slater_storage << std::endl;
    }


    unsigned int num_states = states.NumStates();
    unsigned int size1 = 0;
    unsigned int size2 = 0;
    unsigned int size3 = 0;

    unsigned int i1, i2, i3, i4;
    unsigned int k, kmax;

    ConstStateIterator it_1 = states.GetConstStateIterator();
    ConstStateIterator it_2 = states.GetConstStateIterator();
    ConstStateIterator it_3 = states.GetConstStateIterator();
    ConstStateIterator it_4 = states.GetConstStateIterator();

    // One electron integrals
    it_1.First();
    while(!it_1.AtEnd())
    {
        it_2 = it_1;
        while(!it_2.AtEnd())
        {
            const Orbital* si = it_1.GetState();
            const Orbital* sj = it_2.GetState();

            // Calculate any remaining one electron integrals
            if(si->Kappa() == sj->Kappa())
                size1++;

            it_2.Next();
        }
        it_1.Next();
    }
    *outstream << "Num one-electron integrals: " << size1 << std::endl;

    // Two-electron integrals
    it_2.First(); i2 = 0;
    while(!it_2.AtEnd())
    {   const Orbital* s2 = it_2.GetState();

        it_4.First(); i4 = 0;
        while(!it_4.AtEnd())
        {   const Orbital* s4 = it_4.GetState();

            // Limits on k
            k = abs(int(s2->L()) - int(s4->L()));
            if(fabs(s2->J() - s4->J()) > double(k))
                k += 2;

            kmax = s2->L() + s4->L();
            if(s2->J() + s4->J() < double(kmax))
                kmax -= 2;

            while(k <= kmax)
            {
                // s1 is the smallest
                it_1.First(); i1 = 0;
                while((i1 <= i2) && (i1 <= i4) && (it_1.GetState()->RequiredPQN() <= max_pqn_1))
                {   const Orbital* s1 = it_1.GetState();
                    
                    it_3 = it_1; i3 = i1;
                    unsigned int i3_limit;
                    if(i1 == i2)
                        i3_limit = i4;
                    else
                        i3_limit = num_states;
                    while((i3 <= i3_limit) && !it_3.AtEnd())
                    {   const Orbital* s3 = it_3.GetState();
                       
                        // Check conditions:
                        //     if(i1 == i4)  i2 <= i3
                        //     if(i1 == i3)  i2 <= i4
                        //     Triangle(s1, s3, k)
                        if(   !((i1 == i4) && (i2 > i3))
                           && !((i1 == i3) && (i2 > i4))
                           && ((s1->L() + s3->L() + k)%2 == 0)
                           && (double(k) >= fabs(s1->J() - s3->J()))
                           && (k <= s1->L() + s3->L())
                           && (double(k) <= s1->J() + s3->J()))
                        {
                            // Check max_pqn conditions - these are slightly redefined
                            //    because of the reduced symmetry in integrals.
                            // (one of remaining states s2, s3 or s4) < max_pqn_2
                            // (at least two  of remaining states) < max_pqn_3
                            unsigned int smallest_pqn = s2->RequiredPQN();
                            unsigned int second_smallest_pqn = s3->RequiredPQN();
                            if(smallest_pqn > second_smallest_pqn)
                                swap(smallest_pqn, second_smallest_pqn);
                            if(s4->RequiredPQN() < smallest_pqn)
                            {   second_smallest_pqn = smallest_pqn;
                                smallest_pqn = s4->RequiredPQN();
                            }
                            else if(s4->RequiredPQN() < second_smallest_pqn)
                                second_smallest_pqn = s4->RequiredPQN();

                            if((smallest_pqn <= max_pqn_2) &&
                               (second_smallest_pqn <= max_pqn_3))
                            {
                                size2++;
                            }
                        }
                        it_3.Next(); i3++;
                    }
                    it_1.Next(); i1++;
                }
                k += 2;
            }
            it_4.Next(); i4++;
        }
        it_2.Next(); i2++;
    }
    *outstream << "Num two-electron integrals: " << size2 << std::endl;

    // Count box diagrams with wrong parity
    if(include_extra_box)
    {
        it_2.First(); i2 = 0;
        while(!it_2.AtEnd())
        {   const Orbital* s2 = it_2.GetState();

            it_4.First(); i4 = 0;
            while(!it_4.AtEnd())
            {   const Orbital* s4 = it_4.GetState();

                // Limits on k; start at wrong parity
                k = (unsigned int)fabs(s2->J() - s4->J());
                if((k + s2->L() + s4->L())%2 == 0)
                    k++;

                kmax = (unsigned int)(s2->J() + s4->J());

                while(k <= kmax)
                {
                    // s1 is the smallest
                    it_1.First(); i1 = 0;
                    while((i1 <= i2) && (i1 <= i4) && (it_1.GetState()->RequiredPQN() <= box_max_pqn_1))
                    {   const Orbital* s1 = it_1.GetState();
                        
                        it_3 = it_1; i3 = i1;
                        unsigned int i3_limit;
                        if(i1 == i2)
                            i3_limit = i4;
                        else
                            i3_limit = NumStates;

                        while((i3 <= i3_limit) && !it_3.AtEnd())
                        {   const Orbital* s3 = it_3.GetState();
                           
                            // Check conditions:
                            //     if(i1 == i4)  i2 <= i3
                            //     if(i1 == i3)  i2 <= i4
                            //     Triangle(s1.J, s3.J, k)
                            if(   !((i1 == i4) && (i2 > i3))
                            && !((i1 == i3) && (i2 > i4))
                            && ((s1->L() + s3->L() + k)%2 == 1)
                            && (double(k) >= fabs(s1->J() - s3->J()))
                            && (double(k) <= s1->J() + s3->J()))
                            {
                                // Check max_pqn conditions - these are slightly redefined
                                //    because of the reduced symmetry in integrals.
                                // (one of remaining states s2, s3 or s4) < max_pqn_2
                                // (at least two  of remaining states) < max_pqn_3
                                unsigned int smallest_pqn = s2->RequiredPQN();
                                unsigned int second_smallest_pqn = s3->RequiredPQN();
                                if(smallest_pqn > second_smallest_pqn)
                                    swap(smallest_pqn, second_smallest_pqn);
                                if(s4->RequiredPQN() < smallest_pqn)
                                {   second_smallest_pqn = smallest_pqn;
                                    smallest_pqn = s4->RequiredPQN();
                                }
                                else if(s4->RequiredPQN() < second_smallest_pqn)
                                    second_smallest_pqn = s4->RequiredPQN();

                                if((smallest_pqn <= box_max_pqn_2) &&
                                (second_smallest_pqn <= box_max_pqn_3))
                                {
                                    size3++;
                                }
                            }
                            it_3.Next(); i3++;
                        }
                        it_1.Next(); i1++;
                    }
                    k += 2;
                }
                it_4.Next(); i4++;
            }
            it_2.Next(); i2++;
        }
        *outstream << "Num extra box integrals: " << size3 << std::endl;
    }

    return (size1 + size2 + size3);
}

void CIIntegralsMBPT::Update()
{
    Update("");
}

void CIIntegralsMBPT::Update(const std::string& sigma_id)
{
    Clear();

    // Update the valence energies and Slater integrals for perturbation theory
    if((include_sigma1 || include_mbpt1 || include_mbpt2 || include_extra_box) && PT)
        PT->UpdateIntegrals(&states);
    if((include_valence_mbpt1 || include_valence_mbpt2 || include_valence_extra_box) && ValencePT)
        ValencePT->UpdateIntegrals(&states);

    UpdateOneElectronIntegrals(sigma_id);
    UpdateTwoElectronIntegrals();
    if(include_extra_box || include_valence_extra_box)
        UpdateTwoElectronBoxDiagrams();
}

void CIIntegralsMBPT::UpdateOneElectronIntegrals(const std::string& sigma_id)
{
    StateIntegrator SI(states.GetLattice());
    unsigned int i, j;
    double alphasquared = PhysicalConstant::Instance()->GetAlphaSquared();

    ConstStateIterator it_i = states.GetConstStateIterator();
    ConstStateIterator it_j = states.GetConstStateIterator();

    // If fp is NULL, calculate one electron integrals, otherwise read them in.
    std::string file1 = read_id + ".one.int";
    FILE* fp = NULL;

    if(!read_id.empty())
    {   fp = fopen(file1.c_str(), "rb");
        if(fp)
        {   ReadOneElectronIntegrals(fp);
            fclose(fp);
        }
    }

    // Calculate sigma1 potentials
    if(include_sigma1)
    {
        std::map<int, SigmaPotential*>::iterator it_sigma = Sigma1.begin();
        while(it_sigma != Sigma1.end())
        {   delete it_sigma->second;
            it_sigma++;
        }
        Sigma1.clear();

        unsigned int max_l = 0;
        it_i.First();
        while(!it_i.AtEnd())
        {   max_l = mmax(it_i.GetState()->L(), max_l);
            it_i.Next();
        }

        for(int kappa = - (int)max_l - 1; kappa <= (int)max_l; kappa++)
            if(kappa != 0)
            {
                std::string sigma_file;
                if(!sigma_id.empty())
                    sigma_file = sigma_id + ".";
                else
                    sigma_file = read_id + ".";

                if(kappa > 0)
                    sigma_file = sigma_file + MathConstant::Instance()->GetSpectroscopicNotation(kappa);
                else
                    sigma_file = sigma_file + MathConstant::Instance()->GetSpectroscopicNotation(-kappa-1);

                if(kappa < -1)
                    sigma_file = sigma_file + '+';

                sigma_file = sigma_file + ".matrix";

                SigmaPotential* pot = new SigmaPotential(states.GetLattice(), sigma_file);

                // If (pot->Size() == 0), then it didn't already exist on disk, so make it.
                if(!pot->Size() && PT)
                {
                    PT->GetSecondOrderSigma(kappa, pot);
                    pot->Store();
                }

                Sigma1.insert(std::pair<int, SigmaPotential*>(kappa, pot));
            }
    }

#ifdef _MPI
    // If MPI, then only calculate our integrals (these will presumably be stored).
    int count = 0;
#endif

    // Want to save progress every hour or so.
    time_t start, gap, now;
    time(&start);
    gap = 47 * 60;  // 47 minutes

    // Get single particle integrals
    it_i.First(); i = 0;
    while(!it_i.AtEnd())
    {
        it_j = it_i; j = i;
        while(!it_j.AtEnd())
        {
            const Orbital* si = it_i.GetState();
            const Orbital* sj = it_j.GetState();

            // calculate one electron integrals
            if(si->Kappa() == sj->Kappa())
            {
                unsigned int key = i * NumStates + j;

                // Check that this integral doesn't already exist
                // (it may have been read in)
                if(OneElectronIntegrals.find(key) == OneElectronIntegrals.end())
                {
                  #ifdef _MPI
                    // MPI: Check if this is our integral.
                    if((!PT && !ValencePT) || (count == ProcessorRank))
                    {
                  #endif
                    double integral = SI.HamiltonianMatrixElement(*si, *sj, *states.GetCore());
                    if(include_sigma1)
                    {   SigmaPotential* pot = Sigma1[si->Kappa()];
                        integral += pot->GetMatrixElement(si->f, sj->f);
                    }
                    else if(include_mbpt1 && PT)
                    {   integral += PT->GetOneElectronDiagrams(si, sj);
                    }
                    if(include_mbpt1_subtraction && PT)
                    {   integral += PT->GetOneElectronSubtraction(si, sj);
                    }
                    if(include_valence_mbpt1 && ValencePT)
                    {   integral += ValencePT->GetOneElectronValence(si, sj);
                    }
                    OneElectronIntegrals.insert(std::pair<unsigned int, double>(key, integral));
                  #ifdef _MPI
                    }
                    count++;
                    if(count == NumProcessors)
                        count = 0;
                  #endif

                    if((include_mbpt1 && PT) || (include_valence_mbpt1 && ValencePT))
                    {   time(&now);
                        if(now - start > gap)
                        {   WriteOneElectronIntegrals();
                            start = now;
                        }
                    }
                }
            }

            // i.pqn <= j.pqn, so calculate using derivative of i instead of j
            SMSIntegrals.insert(std::pair<unsigned int, double>(i* NumStates + j,
                -SI.IsotopeShiftIntegral(*it_j.GetState(), *it_i.GetState())));

            // Overlap integrals
            const SingleParticleWavefunction& p1 = *it_i.GetState();
            const SingleParticleWavefunction& p2 = *it_j.GetState();
            if(p1.L() == p2.L())
            {
                double overlap = 0.;
                const double* dR = states.GetLattice()->dR();
                for(unsigned int x=0; x<mmin(p1.Size(), p2.Size()); x++)
                    overlap += (p1.f[x] * p2.f[x] + alphasquared * p1.g[x] * p2.g[x]) * dR[x];

                OverlapIntegrals.insert(std::pair<unsigned int, double>(i * NumStates + j, overlap));
            }

            it_j.Next(); j++;
        }
        it_i.Next(); i++;
    }

    if((include_mbpt1 && PT) || (include_valence_mbpt1 && ValencePT))
    {
        WriteOneElectronIntegrals();
    #ifdef _MPI
        MPI::COMM_WORLD.Barrier();
    #endif
    }
}

void CIIntegralsMBPT::UpdateTwoElectronIntegrals()
{
    // Read stored two electron integrals.
    std::string file1 = read_id + ".two.int";
    FILE* fp = NULL;

    if(!read_id.empty())
    {   fp = fopen(file1.c_str(), "rb");
        if(fp)
        {   ReadTwoElectronIntegrals(fp);
            fclose(fp);
        }
    }

#ifdef _MPI
    // If MPI, then only calculate our integrals (these will presumably be stored).
    int count = 0;
#endif

    // Want to save progress every hour or so.
    time_t start, gap, now;
    time(&start);
    gap = 47 * 60;  // 47 minutes

    // Calculate any remaining two electron integrals.
    CoulombIntegrator CI(states.GetLattice());
    std::vector<double> density(states.GetCore()->GetHFPotential().size());
    std::vector<double> Pot24(states.GetCore()->GetHFPotential().size());
    const double* dR = states.GetLattice()->dR();

    unsigned int i1, i2, i3, i4;
    unsigned int k, kmax;
    unsigned int p;  // just a counter
    double alphasquared = PhysicalConstant::Instance()->GetAlphaSquared();

    ConstStateIterator it_1 = states.GetConstStateIterator();
    ConstStateIterator it_2 = states.GetConstStateIterator();
    ConstStateIterator it_3 = states.GetConstStateIterator();
    ConstStateIterator it_4 = states.GetConstStateIterator();

    it_2.First(); i2 = 0;
    while(!it_2.AtEnd())
    {   const Orbital* s2 = it_2.GetState();

        it_4.First(); i4 = 0;
        while(!it_4.AtEnd())
        {   const Orbital* s4 = it_4.GetState();

            // Limits on k
            k = abs(int(s2->L()) - int(s4->L()));
            if(fabs(s2->J() - s4->J()) > double(k))
                k += 2;

            kmax = s2->L() + s4->L();
            if(s2->J() + s4->J() < double(kmax))
                kmax -= 2;

            // Get density24
            if(k <= kmax)
            {   for(p=0; p < mmin(s2->Size(), s4->Size()); p++)
                {   density[p] = s2->f[p] * s4->f[p] + alphasquared * s2->g[p] * s4->g[p];
                }
            }

            while(k <= kmax)
            {
                // Get Pot24
                CI.FastCoulombIntegrate(density, Pot24, k, mmin(s2->Size(), s4->Size()));

                // s1 is the smallest
                it_1.First(); i1 = 0;
                while((i1 <= i2) && (i1 <= i4) && (it_1.GetState()->RequiredPQN() <= max_pqn_1))
                {   const Orbital* s1 = it_1.GetState();
                    
                    it_3 = it_1; i3 = i1;
                    unsigned int i3_limit;
                    if(i1 == i2)
                        i3_limit = i4;
                    else
                        i3_limit = NumStates;

                    while((i3 <= i3_limit) && !it_3.AtEnd())
                    {   const Orbital* s3 = it_3.GetState();
                       
                        // Check conditions:
                        //     if(i1 == i4)  i2 <= i3
                        //     if(i1 == i3)  i2 <= i4
                        //     Triangle(s1.J, s3.J, k)
                        if(   !((i1 == i4) && (i2 > i3))
                           && !((i1 == i3) && (i2 > i4))
                           && ((s1->L() + s3->L() + k)%2 == 0)
                           && (double(k) >= fabs(s1->J() - s3->J()))
                           && (k <= s1->L() + s3->L())
                           && (double(k) <= s1->J() + s3->J()))
                        {
                            // Check max_pqn conditions - these are slightly redefined
                            //    because of the reduced symmetry in integrals.
                            // (one of remaining states s2, s3 or s4) < max_pqn_2
                            // (at least two  of remaining states) < max_pqn_3
                            unsigned int smallest_pqn = s2->RequiredPQN();
                            unsigned int second_smallest_pqn = s3->RequiredPQN();
                            if(smallest_pqn > second_smallest_pqn)
                                swap(smallest_pqn, second_smallest_pqn);
                            if(s4->RequiredPQN() < smallest_pqn)
                            {   second_smallest_pqn = smallest_pqn;
                                smallest_pqn = s4->RequiredPQN();
                            }
                            else if(s4->RequiredPQN() < second_smallest_pqn)
                                second_smallest_pqn = s4->RequiredPQN();

                            if((smallest_pqn <= max_pqn_2) &&
                               (second_smallest_pqn <= max_pqn_3))
                            {
                                unsigned int key = k  * NumStates*NumStates*NumStates*NumStates +
                                                   i1 * NumStates*NumStates*NumStates +
                                                   i2 * NumStates*NumStates +
                                                   i3 * NumStates +
                                                   i4;

                                // Check that this integral doesn't already exist
                                // (it may have been read in)
                                if(TwoElectronIntegrals.find(key) == TwoElectronIntegrals.end())
                                {
                                  #ifdef _MPI
                                    // MPI: Check if this is our integral.
                                    if((!PT && !ValencePT) || (count == ProcessorRank))
                                    {
                                  #endif
                                    double radial = 0.;
                                    unsigned int limit = mmin(s1->Size(), s3->Size());
                                    limit = mmin(limit, Pot24.size());
                                    for(p=0; p<limit; p++)
                                    {
                                        radial += (s1->f[p] * s3->f[p] + alphasquared * s1->g[p] * s3->g[p])
                                                * Pot24[p] * dR[p];
                                    }

                                    // Include the SMS directly into the integral
                                    if(include_valence_sms && (k == 1))
                                    {   double SMS = GetNuclearInverseMass();
                                        if(SMS)
                                        {   // i1 <= i3 always, but i2 is not necessarily less than i4.
                                            if(i2 <= i4)
                                                SMS = SMS * SMSIntegrals.find(i2*NumStates + i4)->second;
                                            else
                                                SMS = - SMS * SMSIntegrals.find(i4*NumStates + i2)->second;

                                            radial = radial - SMS * SMSIntegrals.find(i1*NumStates + i3)->second;
                                        }
                                    }

                                    if(include_mbpt2 && PT)
                                    {   radial += PT->GetTwoElectronDiagrams(k, s1, s2, s3, s4);
                                        if(include_mbpt2_subtraction)
                                            radial += PT->GetTwoElectronSubtraction(k, s1, s2, s3, s4);
                                    }
                                    if(include_valence_mbpt2 && ValencePT)
                                    {   radial += ValencePT->GetTwoElectronValence(s1, s2, s3, s4, k);
                                    }

                                    TwoElectronIntegrals.insert(std::pair<unsigned int, double>(key, radial));
                                  #ifdef _MPI
                                    }
                                    count++;
                                    if(count == NumProcessors)
                                        count = 0;
                                  #endif

                                    if((include_mbpt2 && PT) || (include_valence_mbpt2 && ValencePT))
                                    {   time(&now);
                                        if(now - start > gap)
                                        {   WriteTwoElectronIntegrals();
                                            start = now;
                                        }
                                    }
                                }
                            }
                        }
                        it_3.Next(); i3++;
                    }
                    it_1.Next(); i1++;
                }
                k += 2;
            }
            it_4.Next(); i4++;
        }
        it_2.Next(); i2++;
    }

    if((include_mbpt2 && PT) || (include_valence_mbpt2 && ValencePT))
    {
        WriteTwoElectronIntegrals();
    #ifdef _MPI
        MPI::COMM_WORLD.Barrier();
    #endif
    }
}

void CIIntegralsMBPT::UpdateTwoElectronBoxDiagrams()
{
    if(!(include_extra_box && PT) && !(include_valence_extra_box && ValencePT))
        return;

#ifdef _MPI
    // If MPI, then only calculate our integrals (these will presumably be stored).
    int count = 0;
#endif

    // Want to save progress every hour or so.
    time_t start, gap, now;
    time(&start);
    gap = 47 * 60;  // 47 minutes

    unsigned int i1, i2, i3, i4;
    unsigned int k, kmax;

    ConstStateIterator it_1 = states.GetConstStateIterator();
    ConstStateIterator it_2 = states.GetConstStateIterator();
    ConstStateIterator it_3 = states.GetConstStateIterator();
    ConstStateIterator it_4 = states.GetConstStateIterator();

    it_2.First(); i2 = 0;
    while(!it_2.AtEnd())
    {   const Orbital* s2 = it_2.GetState();

        it_4.First(); i4 = 0;
        while(!it_4.AtEnd())
        {   const Orbital* s4 = it_4.GetState();

            // Limits on k; start at wrong parity
            k = (unsigned int)fabs(s2->J() - s4->J());
            if((k + s2->L() + s4->L())%2 == 0)
                k++;

            kmax = (unsigned int)(s2->J() + s4->J());

            while(k <= kmax)
            {
                // s1 is the smallest
                it_1.First(); i1 = 0;
                while((i1 <= i2) && (i1 <= i4) && (it_1.GetState()->RequiredPQN() <= box_max_pqn_1))
                {   const Orbital* s1 = it_1.GetState();
                    
                    it_3 = it_1; i3 = i1;
                    unsigned int i3_limit;
                    if(i1 == i2)
                        i3_limit = i4;
                    else
                        i3_limit = NumStates;

                    while((i3 <= i3_limit) && !it_3.AtEnd())
                    {   const Orbital* s3 = it_3.GetState();
                       
                        // Check conditions:
                        //     if(i1 == i4)  i2 <= i3
                        //     if(i1 == i3)  i2 <= i4
                        //     Triangle(s1.J, s3.J, k)
                        if(   !((i1 == i4) && (i2 > i3))
                           && !((i1 == i3) && (i2 > i4))
                           && ((s1->L() + s3->L() + k)%2 == 1)
                           && (double(k) >= fabs(s1->J() - s3->J()))
                           && (double(k) <= s1->J() + s3->J()))
                        {
                            // Check max_pqn conditions - these are slightly redefined
                            //    because of the reduced symmetry in integrals.
                            // (one of remaining states s2, s3 or s4) < max_pqn_2
                            // (at least two  of remaining states) < max_pqn_3
                            unsigned int smallest_pqn = s2->RequiredPQN();
                            unsigned int second_smallest_pqn = s3->RequiredPQN();
                            if(smallest_pqn > second_smallest_pqn)
                                swap(smallest_pqn, second_smallest_pqn);
                            if(s4->RequiredPQN() < smallest_pqn)
                            {   second_smallest_pqn = smallest_pqn;
                                smallest_pqn = s4->RequiredPQN();
                            }
                            else if(s4->RequiredPQN() < second_smallest_pqn)
                                second_smallest_pqn = s4->RequiredPQN();

                            if((smallest_pqn <= box_max_pqn_2) &&
                               (second_smallest_pqn <= box_max_pqn_3))
                            {
                                unsigned int key = k  * NumStates*NumStates*NumStates*NumStates +
                                                   i1 * NumStates*NumStates*NumStates +
                                                   i2 * NumStates*NumStates +
                                                   i3 * NumStates +
                                                   i4;

                                // Check that this integral doesn't already exist
                                // (it may have been read in)
                                if(TwoElectronIntegrals.find(key) == TwoElectronIntegrals.end())
                                {
                                  #ifdef _MPI
                                    // MPI: Check if this is our integral.
                                    if(count == ProcessorRank)
                                    {
                                  #endif
                                        double radial = 0;
                                        if(include_extra_box && PT)
                                            radial += PT->GetTwoElectronBoxDiagrams(k, s1, s2, s3, s4);
                                        if(include_valence_extra_box && ValencePT)
                                            radial += ValencePT->GetTwoElectronBoxValence(s1, s2, s3, s4, k);

                                        TwoElectronIntegrals.insert(std::pair<unsigned int, double>(key, radial));
                                  #ifdef _MPI
                                    }
                                    count++;
                                    if(count == NumProcessors)
                                        count = 0;
                                  #endif

                                    time(&now);
                                    if(now - start > gap)
                                    {   WriteTwoElectronIntegrals();
                                        start = now;
                                    }
                                }
                            }
                        }
                        it_3.Next(); i3++;
                    }
                    it_1.Next(); i1++;
                }
                k += 2;
            }
            it_4.Next(); i4++;
        }
        it_2.Next(); i2++;
    }

    WriteTwoElectronIntegrals();
#ifdef _MPI
    MPI::COMM_WORLD.Barrier();
#endif
}

/** GetTwoElectronIntegral(k, i, j, l, m) = R_k(ij, lm): i->l, j->m */
double CIIntegralsMBPT::GetTwoElectronIntegral(unsigned int k, const OrbitalInfo& s1, const OrbitalInfo& s2, const OrbitalInfo& s3, const OrbitalInfo& s4) const
{
    unsigned int i1 = state_index.find(s1)->second;
    unsigned int i2 = state_index.find(s2)->second;
    unsigned int i3 = state_index.find(s3)->second;
    unsigned int i4 = state_index.find(s4)->second;

    TwoElectronIntegralOrdering(i1, i2, i3, i4);

    unsigned int key = k  * NumStates*NumStates*NumStates*NumStates +
                       i1 * NumStates*NumStates*NumStates +
                       i2 * NumStates*NumStates +
                       i3 * NumStates +
                       i4;

    double radial = 0.;
    if(TwoElectronIntegrals.find(key) != TwoElectronIntegrals.end())
    {
        radial = TwoElectronIntegrals.find(key)->second;
    }
    else
    {   // Check triangle and parity conditions on k
        if(((k + s1.L() + s3.L())%2 == 1) ||
             (double(k) < fabs(s1.J() - s3.J())) ||
             (double(k) > s1.J() + s3.J()) ||
           ((k + s2.L() + s4.L())%2 == 1) ||
             (double(k) < fabs(s2.J() - s4.J())) ||
             (double(k) > s2.J() + s4.J()))
            return 0.;

        const SingleParticleWavefunction* s_1 = states.GetState(reverse_state_index.find(i1)->second);
        const SingleParticleWavefunction* s_2 = states.GetState(reverse_state_index.find(i2)->second);
        const SingleParticleWavefunction* s_3 = states.GetState(reverse_state_index.find(i3)->second);
        const SingleParticleWavefunction* s_4 = states.GetState(reverse_state_index.find(i4)->second);

        unsigned int p;
        CoulombIntegrator CI(states.GetLattice());
        const double* R = states.GetLattice()->R();
        const double* dR = states.GetLattice()->dR();
        const double core_pol = states.GetCore()->GetPolarisability();
        const double core_rad = states.GetCore()->GetClosedShellRadius();
        double alphasquared = PhysicalConstant::Instance()->GetAlphaSquared();

        // Get density24
        std::vector<double> density(mmin(s_2->Size(), s_4->Size()));
        for(p=0; p<density.size(); p++)
        {
            density[p] = s_2->f[p] * s_4->f[p] + alphasquared * s_2->g[p] * s_4->g[p];
        }
        density.resize(states.GetCore()->GetHFPotential().size());

        // Get Pot24
        std::vector<double> Pot24(density.size());
        CI.FastCoulombIntegrate(density, Pot24, k);

        unsigned int limit = mmin(s_1->Size(), s_3->Size());
        limit = mmin(limit, Pot24.size());
        for(p=0; p<limit; p++)
        {
            radial += (s_1->f[p] * s_3->f[p] + alphasquared * s_1->g[p] * s_3->g[p])
                        * Pot24[p] * dR[p];
        }

        if(core_pol && k == 1)
        {
            double R1 = 0.;
            double R2 = 0.;
            for(p=0; p<limit; p++)
            {
                double r2 = R[p]*R[p] + core_rad*core_rad;
                R1 += (s_1->f[p] * s_3->f[p] + alphasquared * s_1->g[p] * s_3->g[p])/r2 * dR[p];
                R2 += density[p]/r2 * dR[p];
            }

            radial -= core_pol * R1 * R2;
        }

        if(include_valence_sms && (k == 1))
        {   double SMS = GetNuclearInverseMass();
            if(SMS)
            {   // i1 <= i3 always, but i2 is not necessarily less than i4.
                if(i2 <= i4)
                    SMS = SMS * SMSIntegrals.find(i2*NumStates + i4)->second;
                else
                    SMS = - SMS * SMSIntegrals.find(i4*NumStates + i2)->second;

                radial = radial - SMS * SMSIntegrals.find(i1*NumStates + i3)->second;
            }
        }
    }
    return radial;
}

void CIIntegralsMBPT::ReadMultipleOneElectronIntegrals(const std::string& name, unsigned int num_files)
{
    FILE* fp;
    for(unsigned int i = 0; i < num_files; i++)
    {
        std::stringstream ss;
        ss << i;
        std::string filename = name;
        if(num_files > 1)
            filename = filename + '_' + ss.str();
        filename += ".one.int";

        fp = fopen(filename.c_str(), "rb");
        if(fp)
        {   ReadOneElectronIntegrals(fp);
            fclose(fp);
        }
        else
        {   *errstream << "Missing file: " << filename << std::endl;
        }
    }
}

void CIIntegralsMBPT::ReadMultipleTwoElectronIntegrals(const std::string& name, unsigned int num_files)
{
    FILE* fp;
    for(unsigned int i = 0; i < num_files; i++)
    {
        std::stringstream ss;
        ss << i;
        std::string filename = name;
        if(num_files > 1)
            filename = filename + '_' + ss.str();
        filename += ".two.int";

        fp = fopen(filename.c_str(), "rb");
        if(fp)
        {   ReadTwoElectronIntegrals(fp);
            fclose(fp);
        }
        else
        {   *errstream << "Missing file: " << filename << std::endl;
        }
    }
}

void CIIntegralsMBPT::WriteSigmaPotentials() const
{
    std::map<int, SigmaPotential*>::const_iterator it = Sigma1.begin();

    while(it != Sigma1.end())
    {   it->second->Store();
        it++;
    }
}

bool CIIntegralsMBPT::TwoElectronIntegralOrdering(unsigned int& i1, unsigned int& i2, unsigned int& i3, unsigned int& i4) const
{
    // Ordering of indices:
    // i1 is smallest && (if i1 == i2, then (i3 <= i4))
    //                && (if i1 == i3, then (i2 <= i4))
    //                && (if i1 == i4, then (i2 <= i3))

    // Assert one of i1, i3 is smallest
    if(mmin(i1, i3) > mmin(i2, i4))
    {   swap(i1, i2);
        swap(i3, i4);
    }
    // Assert i1 <= i3
    if(i1 > i3)
    {   swap(i1, i3);
        swap(i2, i4);
    }

    if((i1 == i2) && (i4 < i3))
        swap(i3, i4);
    if((i1 == i3) && (i4 < i2))
        swap(i2, i4);
    if((i1 == i4) && (i3 < i2))
        swap(i2, i3);

    return true;
}

void CIIntegralsMBPT::AddSMSToTwoElectronIntegrals(const std::string& name)
{
    const double SMS_scaling = GetNuclearInverseMass();
    if(!include_valence_sms || !SMS_scaling)
        return;

    unsigned int stored_key;
    double radial;
    unsigned int k, i1, i2, i3, i4;

    Clear();
    UpdateOneElectronIntegrals("");
    ReadMultipleTwoElectronIntegrals(name, 1);

    std::map<unsigned int, double>::iterator it = TwoElectronIntegrals.begin();

    while(it != TwoElectronIntegrals.end())
    {
        stored_key = it->first;
        radial = it->second;

        // Get k and stored states
        i4 = stored_key%NumStates;
            stored_key = stored_key/NumStates;
        i3 = stored_key%NumStates;
            stored_key = stored_key/NumStates;
        i2 = stored_key%NumStates;
            stored_key = stored_key/NumStates;
        i1 = stored_key%NumStates;
            stored_key = stored_key/NumStates;
        k = stored_key%NumStates;

        if(k == 1)
        {
            // Get states
            OrbitalInfo& s1 = reverse_state_index.find(i1)->second;
            OrbitalInfo& s2 = reverse_state_index.find(i2)->second;
            OrbitalInfo& s3 = reverse_state_index.find(i3)->second;
            OrbitalInfo& s4 = reverse_state_index.find(i4)->second;

            // Check for correct parity (as opposed to box diagrams)
            if(((s1.L() + s3.L() + k)%2 == 0) && ((s2.L() + s4.L() + k)%2 == 0))
            {
                double SMS = SMS_scaling;

                // i1 <= i3 always, but i2 is not necessarily less than i4.
                if(i2 <= i4)
                    SMS = SMS * SMSIntegrals.find(i2*NumStates + i4)->second;
                else
                    SMS = - SMS * SMSIntegrals.find(i4*NumStates + i2)->second;

                radial = radial - SMS * SMSIntegrals.find(i1*NumStates + i3)->second;
                it->second = radial;
            }
        }

        it++;
    }

    WriteTwoElectronIntegrals();
}

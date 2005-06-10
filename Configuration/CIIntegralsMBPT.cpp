#include "Include.h"
#include "CIIntegralsMBPT.h"
#include "HartreeFock/StateIntegrator.h"
#include "Universal/CoulombIntegrator.h"

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
    unsigned int size1 = 0;
    unsigned int size2 = 0;

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
            const DiscreteState* si = it_1.GetState();
            const DiscreteState* sj = it_2.GetState();

            // Calculate any remaining one electron integrals
            if(si->Kappa() == sj->Kappa())
                size1++;

            it_2.Next();
        }
        it_1.Next();
    }
    *logstream << "Num one-electron integrals: " << size1 << std::endl;

    // Two-electron integrals
    // s1 is the smallest
    it_1.First(); i1 = 0;
    while(!it_1.AtEnd() && (it_1.GetState()->RequiredPQN() <= max_pqn_1))
    {   const DiscreteState* s1 = it_1.GetState();

        it_2 = it_1; i2 = i1;
        while(!it_2.AtEnd())
        {   const DiscreteState* s2 = it_2.GetState();
            
            it_3 = it_1; i3 = i1;
            while(!it_3.AtEnd())
            {   const DiscreteState* s3 = it_3.GetState();

                if(i1 == i2)
                {   it_4 = it_3; i4 = i3;
                }
                else if(i1 == i3)
                {   it_4 = it_2; i4 = i2;
                }
                else
                {   it_4 = it_1; i4 = i1;
                }

                while(!it_4.AtEnd())
                {   const DiscreteState* s4 = it_4.GetState();

                    // Skip if(i1 == i4 && i2 > i3)
                    if(!((i1 == i4) && (i2 > i3))
                       && ((s1->L() + s3->L())%2 == (s2->L() + s4->L())%2))
                    {
                        // Limits on k
                        k = abs(int(s2->L()) - int(s4->L()));
                        k = mmax(k, (unsigned int)abs(int(s1->L()) - int(s3->L())));
                        if(fabs(s2->J() - s4->J()) > double(k))
                            k += 2;
                        if(fabs(s1->J() - s3->J()) > double(k))
                            k += 2;

                        kmax = s2->L() + s4->L();
                        kmax = mmin(kmax, s1->L() + s3->L());
                        if(s2->J() + s4->J() < double(kmax))
                            kmax -= 2;
                        if(s1->J() + s3->J() < double(kmax))
                            kmax -= 2;

                        while(k <= kmax)
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
                            k += 2;
                        }
                    }
                    it_4.Next(); i4++;
                }
                it_3.Next(); i3++;
            }
            it_2.Next(); i2++;
        }
        it_1.Next(); i1++;
    }

    *logstream << "Num two-electron integrals: " << size2 << std::endl;
    return (size1 + size2);
}

void CIIntegralsMBPT::Update()
{
    Update("");
}

void CIIntegralsMBPT::Update(const std::string& sigma_id)
{
    NumStates = states.NumStates();
    state_index.clear();
    reverse_state_index.clear();
    OneElectronIntegrals.clear();
    SMSIntegrals.clear();
    TwoElectronIntegrals.clear();
    OverlapIntegrals.clear();

    UpdateStateIndexes();
    UpdateOneElectronIntegrals(sigma_id);
    UpdateTwoElectronIntegrals();
}

void CIIntegralsMBPT::UpdateOneElectronIntegrals(const std::string& sigma_id)
{
    StateIntegrator SI(*states.GetLattice());
    unsigned int i, j;

    ConstStateIterator it_i = states.GetConstStateIterator();
    ConstStateIterator it_j = states.GetConstStateIterator();

    // If fp is NULL, calculate one electron integrals, otherwise read them in.
    std::string file1 = id + ".one.int";
    FILE* fp = NULL;

    if(!id.empty())
    {   fp = fopen(file1.c_str(), "rb");
        if(fp)
        {   ReadOneElectronIntegrals(fp);
            fclose(fp);
        }
    }

    // Update the valence energies for perturbation theory
    if(include_sigma1 || include_mbpt1)
        SetValenceEnergies();

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
                    sigma_file = id + ".";

                if(kappa > 0)
                    sigma_file = sigma_file + Constant::SpectroscopicNotation[kappa];
                else
                    sigma_file = sigma_file + Constant::SpectroscopicNotation[-kappa-1];

                if(kappa < -1)
                    sigma_file = sigma_file + '+';

                sigma_file = sigma_file + ".matrix";

                SigmaPotential* pot = new SigmaPotential(states.GetLattice(), sigma_file);

                // If (pot->Size() == 0), then it didn't already exist on disk, so make it.
                if(!pot->Size() && PT)
                {
                    double valence_energy = ValenceEnergies[kappa];

                    PT->UseBrillouinWignerPT(valence_energy);
                    PT->GetSecondOrderSigma(kappa, pot);
                }

                Sigma1.insert(std::pair<int, SigmaPotential*>(kappa, pot));
            }
    }

#ifdef _MPI
    // If MPI, then only calculate our integrals (these will presumably be stored).
    int count = 0;
#endif

    // Get single particle integrals
    it_i.First(); i = 0;
    while(!it_i.AtEnd())
    {
        it_j = it_i; j = i;
        while(!it_j.AtEnd())
        {
            const DiscreteState* si = it_i.GetState();
            const DiscreteState* sj = it_j.GetState();

            // calculate one electron integrals
            if(si->Kappa() == sj->Kappa())
            {
              #ifdef _MPI
                // MPI: Check if this is our integral.
                if(count == ProcessorRank)
                {
              #endif
                unsigned int key = i * NumStates + j;

                // Check that this integral doesn't already exist
                // (it may have been read in)
                if(OneElectronIntegrals.find(key) == OneElectronIntegrals.end())
                {
                    double integral = SI.HamiltonianMatrixElement(*si, *sj, *states.GetCore());
                    if(include_sigma1)
                    {   SigmaPotential* pot = Sigma1[si->Kappa()];
                        integral += pot->GetMatrixElement(si->f, sj->f);
                    }
                    else if(include_mbpt1 && PT)
                    {   PT->UseBrillouinWignerPT(ValenceEnergies[si->Kappa()]);
                        integral += PT->GetSecondOrderSigma(si, sj);
                    }
                    if(include_mbpt1_subtraction && PT)
                    {   integral += PT->GetSigmaSubtraction(si, sj);
                    }
                    OneElectronIntegrals.insert(std::pair<unsigned int, double>(key, integral));
                }
              #ifdef _MPI
                }
                count++;
                if(count == NumProcessors)
                    count = 0;
              #endif
            }

            // i.pqn <= j.pqn, so calculate using derivative of i instead of j
            SMSIntegrals.insert(std::pair<unsigned int, double>(i* NumStates + j,
                -SI.IsotopeShiftIntegral(*it_j.GetState(), *it_i.GetState())));

            // Overlap integrals
            const State& p1 = *it_i.GetState();
            const State& p2 = *it_j.GetState();
            if(p1.L() == p2.L())
            {
                double overlap = 0.;
                const double* dR = states.GetLattice()->dR();
                for(unsigned int x=0; x<mmin(p1.Size(), p2.Size()); x++)
                    overlap += (p1.f[x] * p2.f[x] + Constant::AlphaSquared * p1.g[x] * p2.g[x]) * dR[x];

                OverlapIntegrals.insert(std::pair<unsigned int, double>(i * NumStates + j, overlap));
            }

            it_j.Next(); j++;
        }
        it_i.Next(); i++;
    }
}

void CIIntegralsMBPT::UpdateTwoElectronIntegrals()
{
    // Read stored two electron integrals.
    std::string file1 = id + ".two.int";
    FILE* fp = NULL;

    if(!id.empty())
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

    // Calculate any remaining two electron integrals.
    CoulombIntegrator CI(*states.GetLattice());
    const double* dR = states.GetLattice()->dR();
    const double* R = states.GetLattice()->R();
    const double core_pol = states.GetCore()->GetPolarisability();
    const double core_rad = states.GetCore()->GetClosedShellRadius();

    unsigned int i1, i2, i3, i4;
    unsigned int k, kmax;
    unsigned int p;  // just a counter

    ConstStateIterator it_1 = states.GetConstStateIterator();
    ConstStateIterator it_2 = states.GetConstStateIterator();
    ConstStateIterator it_3 = states.GetConstStateIterator();
    ConstStateIterator it_4 = states.GetConstStateIterator();

    // s1 is the smallest
    it_1.First(); i1 = 0;
    while(!it_1.AtEnd() && (it_1.GetState()->RequiredPQN() <= max_pqn_1))
    {   const DiscreteState* s1 = it_1.GetState();

        it_2 = it_1; i2 = i1;
        while(!it_2.AtEnd())
        {   const DiscreteState* s2 = it_2.GetState();
            
            it_3 = it_1; i3 = i1;
            while(!it_3.AtEnd())
            {   const DiscreteState* s3 = it_3.GetState();

                if(i1 == i2)
                {   it_4 = it_3; i4 = i3;
                }
                else if(i1 == i3)
                {   it_4 = it_2; i4 = i2;
                }
                else
                {   it_4 = it_1; i4 = i1;
                }

                while(!it_4.AtEnd())
                {   const DiscreteState* s4 = it_4.GetState();

                    // Skip if(i1 == i4 && i2 > i3)
                    if(!((i1 == i4) && (i2 > i3))
                       && ((s1->L() + s3->L())%2 == (s2->L() + s4->L())%2))
                    {
                        // Limits on k
                        k = abs(int(s2->L()) - int(s4->L()));
                        k = mmax(k, (unsigned int)abs(int(s1->L()) - int(s3->L())));
                        if(fabs(s2->J() - s4->J()) > double(k))
                            k += 2;
                        if(fabs(s1->J() - s3->J()) > double(k))
                            k += 2;

                        kmax = s2->L() + s4->L();
                        kmax = mmin(kmax, s1->L() + s3->L());
                        if(s2->J() + s4->J() < double(kmax))
                            kmax -= 2;
                        if(s1->J() + s3->J() < double(kmax))
                            kmax -= 2;

                        // Get density24
                        std::vector<double> density(mmin(s2->Size(), s4->Size()));
                        if(k <= kmax)
                        {   for(p=0; p<density.size(); p++)
                            {   density[p] = s2->f[p] * s4->f[p] + Constant::AlphaSquared * s2->g[p] * s4->g[p];
                            }
                            density.resize(states.GetCore()->GetHFPotential().size());
                        }

                        while(k <= kmax)
                        {
                            // Get Pot24
                            std::vector<double> Pot24(density.size());
                            CI.FastCoulombIntegrate(density, Pot24, k);

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
                              #ifdef _MPI
                                // MPI: Check if this is our integral.
                                if(count == ProcessorRank)
                                {
                              #endif
                                unsigned int key = k  * NumStates*NumStates*NumStates*NumStates +
                                                   i1 * NumStates*NumStates*NumStates +
                                                   i2 * NumStates*NumStates +
                                                   i3 * NumStates +
                                                   i4;

                                // Check that this integral doesn't already exist
                                // (it may have been read in)
                                if(TwoElectronIntegrals.find(key) == TwoElectronIntegrals.end())
                                {
                                    double radial = 0.;
                                    unsigned int limit = mmin(s1->Size(), s3->Size());
                                    limit = mmin(limit, Pot24.size());
                                    for(p=0; p<limit; p++)
                                    {
                                        radial += (s1->f[p] * s3->f[p] + Constant::AlphaSquared * s1->g[p] * s3->g[p])
                                                * Pot24[p] * dR[p];
                                    }

                                    if(include_mbpt2 && PT)
                                    {   PT->SetTwoParticleEnergy(ValenceEnergies[s1->Kappa()] + ValenceEnergies[s2->Kappa()]);
                                        radial += PT->GetTwoElectronDiagrams(s1, s2, s3, s4, k);
                                    }

                                    TwoElectronIntegrals.insert(std::pair<unsigned int, double>(key, radial));
                                }
                              #ifdef _MPI
                                }
                                count++;
                                if(count == NumProcessors)
                                    count = 0;
                              #endif
                            }
                            k += 2;
                        }
                    }
                    it_4.Next(); i4++;
                }
                it_3.Next(); i3++;
            }
            it_2.Next(); i2++;
        }
        it_1.Next(); i1++;
    }
}

/** GetTwoElectronIntegral(k, i, j, l, m) = R_k(ij, lm): i->l, j->m */
double CIIntegralsMBPT::GetTwoElectronIntegral(unsigned int k, const StateInfo& s1, const StateInfo& s2, const StateInfo& s3, const StateInfo& s4) const
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
    {   // Check triangle conditions on k
        if((int(k) < abs(int(s1.L()) - int(s3.L()))) ||
             (double(k) < fabs(s1.J() - s3.J())) ||
           (k > s1.L() + s3.L()) ||
             (double(k) > s1.J() + s3.J()) ||
           (int(k) < abs(int(s2.L()) - int(s4.L()))) ||
             (double(k) < fabs(s2.J() - s4.J())) ||
           (k > s2.L() + s4.L()) ||
             (double(k) > s2.J() + s4.J()))
            return 0.;

        const State* s_1 = states.GetState(reverse_state_index.find(i1)->second);
        const State* s_2 = states.GetState(reverse_state_index.find(i2)->second);
        const State* s_3 = states.GetState(reverse_state_index.find(i3)->second);
        const State* s_4 = states.GetState(reverse_state_index.find(i4)->second);

        unsigned int p;
        CoulombIntegrator CI(*states.GetLattice());
        const double* R = states.GetLattice()->R();
        const double* dR = states.GetLattice()->dR();
        const double core_pol = states.GetCore()->GetPolarisability();
        const double core_rad = states.GetCore()->GetClosedShellRadius();

        // Get density24
        std::vector<double> density(mmin(s_2->Size(), s_4->Size()));
        for(p=0; p<density.size(); p++)
        {
            density[p] = s_2->f[p] * s_4->f[p] + Constant::AlphaSquared * s_2->g[p] * s_4->g[p];
        }
        density.resize(states.GetCore()->GetHFPotential().size());

        // Get Pot24
        std::vector<double> Pot24(density.size());
        CI.FastCoulombIntegrate(density, Pot24, k);

        unsigned int limit = mmin(s_1->Size(), s_3->Size());
        limit = mmin(limit, Pot24.size());
        for(p=0; p<limit; p++)
        {
            radial += (s_1->f[p] * s_3->f[p] + Constant::AlphaSquared * s_1->g[p] * s_3->g[p])
                        * Pot24[p] * dR[p];
        }

        if(core_pol && k == 1)
        {
            double R1 = 0.;
            double R2 = 0.;
            for(p=0; p<limit; p++)
            {
                double r2 = R[p]*R[p] + core_rad*core_rad;
                R1 += (s_1->f[p] * s_3->f[p] + Constant::AlphaSquared * s_1->g[p] * s_3->g[p])/r2 * dR[p];
                R2 += density[p]/r2 * dR[p];
            }

            radial -= core_pol * R1 * R2;
        }
    }

    if(include_valence_sms && (k == 1))
    {   double SMS = GetNuclearInverseMass();
        if(SMS)
        {   // i1 <= i3 always, but i2 is not necessarily less than i4.
            if(i2 <= i4)
                SMS = SMS * SMSIntegrals.find(i2*NumStates + i4)->second;
            else
                SMS = SMS * SMSIntegrals.find(i4*NumStates + i2)->second;

            radial = radial - SMS * SMSIntegrals.find(i1*NumStates + i3)->second;
        }
    }

    return radial;
}

void CIIntegralsMBPT::ReadMultipleOneElectronIntegrals(const std::string& name, unsigned int num_files)
{
    UpdateStateIndexes();

    FILE* fp;
    for(unsigned int i = 0; i < num_files; i++)
    {
        std::stringstream ss;
        ss << i;
        std::string filename = name + '_' + ss.str() + ".one.int";

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
    UpdateStateIndexes();

    FILE* fp;
    for(unsigned int i = 0; i < num_files; i++)
    {
        std::stringstream ss;
        ss << i;
        std::string filename = name + '_' + ss.str() + ".two.int";

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

void CIIntegralsMBPT::SetValenceEnergies()
{
    ConstStateIterator it_i = states.GetConstStateIterator();
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

bool CIIntegralsMBPT::TwoElectronIntegralOrdering(unsigned int& i1, unsigned int& i2, unsigned int& i3, unsigned int& i4) const
{
    // Ordering of indices:
    // i1 is smallest && (if i1 == i2, then (i3 <= i4))
    //                && (if i1 == i3, then (i2 <= i4))
    //                && (if i1 == i4, then (i2 <= i3))
    unsigned int min = i1;
    unsigned int minpos = 1;
    if(i2 < min)
    {   min = i2;
        minpos = 2;
    }
    if(i3 < min)
    {   min = i3;
        minpos = 3;
    }
    if(i4 < min)
    {   min = i4;
        minpos = 4;
    }

    switch(minpos)
    {   case 1:
            break;
        case 2:
            swap(i2, i1);
            swap(i4, i3);
            break;
        case 3:
            swap(i3, i1);
            swap(i4, i2);
            break;
        case 4:
            swap(i4, i1);
            swap(i3, i2);
            break;
    }

    if((i1 == i2) && (i4 < i3))
        swap(i3, i4);
    if((i1 == i3) && (i4 < i2))
        swap(i2, i4);
    if((i1 == i4) && (i3 < i2))
        swap(i2, i3);

    return true;
}

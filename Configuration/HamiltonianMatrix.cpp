#include "Include.h"
#include "HamiltonianMatrix.h"
#include "SmallMatrix.h"
#include "SymMatrix.h"
#include "HartreeFock/StateIntegrator.h"
#include "Universal/CoulombIntegrator.h"
#include "HartreeFock/State.h"
#include "Universal/Eigensolver.h"
#include "ConfigGenerator.h"

#define SMALL_MATRIX_LIM 2000

HamiltonianMatrix::HamiltonianMatrix(const ExcitedStates& excited_states, const RelativisticConfigList& rconfigs):
    states(excited_states), configs(rconfigs), include_sms_v2(false), NumSolutions(0)
{
    // Set up matrix
    N = 0;
    RelativisticConfigList::const_iterator it = configs.begin();
    while(it != configs.end())
    {   N += it->GetJCoefficients().size();
        it++;
    }

    std::cerr << " " << N << std::flush;
    UpdateIntegrals();

    if(N <= SMALL_MATRIX_LIM)
        M = new SmallMatrix(N);
    else
        M = new SymMatrix(N);

    std::cout << " Number of J-configurations = " << N << std::endl;
}

void HamiltonianMatrix::UpdateIntegrals()
{
    // One electron integrals
    StateIntegrator SI(*states.GetLattice());
    CoulombIntegrator CI(*states.GetLattice());

    NumStates = states.NumStates();
    state_index.clear();
    OneElectronIntegrals.clear();
    SMSIntegrals.clear();
    TwoElectronIntegrals.clear();

    unsigned int i, j;

    ConstStateIterator it_i = states.GetConstStateIterator();
    ConstStateIterator it_j = states.GetConstStateIterator();

    it_i.First(); i = 0;
    while(!it_i.AtEnd())
    {   
        // Use it_i to build index of states
        state_index.insert(std::pair<StateInfo, unsigned int>(StateInfo(it_i.GetState()), i));
        it_j = it_i; j = i;
        while(!it_j.AtEnd())
        {
            OneElectronIntegrals.insert(std::pair<unsigned int, double>(i* NumStates + j,
                SI.HamiltonianMatrixElement(*it_i.GetState(), *it_j.GetState(), *states.GetCore())));

            // SMS integrals too
            SMSIntegrals.insert(std::pair<unsigned int, double>(i* NumStates + j,
                CI.IsotopeShiftIntegral(*it_i.GetState(), *it_j.GetState())));

            it_j.Next(); j++;
        }
        it_i.Next(); i++;
    }

    // Two electron integrals
    const double* dR = states.GetLattice()->dR();

    unsigned int i1, i2, i3, i4;
    unsigned int k, kmax;
    unsigned int p;  // just a counter

    ConstStateIterator it_1 = states.GetConstStateIterator();
    ConstStateIterator it_2 = states.GetConstStateIterator();
    ConstStateIterator it_3 = states.GetConstStateIterator();
    ConstStateIterator it_4 = states.GetConstStateIterator();

    // Get 2 -> 4
    it_2.First(); i2 = 0;
    while(!it_2.AtEnd())
    {
        const State* s_2 = it_2.GetState();
        it_4 = it_2; i4 = i2;
        while(!it_4.AtEnd())
        {
            const State* s_4 = it_4.GetState();

            // Limits on k
            k = abs(int(s_2->L()) - int(s_4->L()));
            if(fabs(s_2->J() - s_4->J()) > double(k))
                k += 2;

            kmax = s_2->L() + s_4->L();
            if(s_2->J() + s_4->J() < double(kmax))
                kmax -= 2;

            // Get density24
            std::vector<double> density(mmin(s_2->Size(), s_4->Size()));
            if(k <= kmax)
            {   
                for(p=0; p<density.size(); p++)
                {
                    density[p] = s_2->f[p] * s_4->f[p] + Constant::AlphaSquared * s_2->g[p] * s_4->g[p];
                }
                density.resize(states.GetCore()->GetHFPotential().size());
            }

            while(k <= kmax)
            {
                // Get Pot24
                std::vector<double> Pot24(density.size());
                CI.FastCoulombIntegrate(density, Pot24, k);

                // s1 is the smallest
                it_1.First(); i1 = 0;
                while(i1 <= i2)
                {
                    const State* s_1 = it_1.GetState();

                    it_3 = it_1; i3 = i1;
                    unsigned int i3_limit;
                    if(i1 == i2)
                        i3_limit = i4;
                    else
                        i3_limit = NumStates;
                    while((i3 <= i3_limit) && !it_3.AtEnd())
                    {
                        const State* s_3 = it_3.GetState();

                        if((int(k) >= abs(int(s_1->L()) - int(s_3->L()))) &&
                           (double(k) >= fabs(s_1->J() - s_3->J())) &&
                           (k <= s_1->L() + s_3->L()) &&
                           (double(k) <= s_1->J() + s_3->J()))
                        {
                            unsigned int key = k  * NumStates*NumStates*NumStates*NumStates +
                                               i1 * NumStates*NumStates*NumStates +
                                               i2 * NumStates*NumStates + 
                                               i3 * NumStates +
                                               i4;

                            double radial = 0.;
                            unsigned int limit = mmin(s_1->Size(), s_3->Size());
                            for(p=0; p<limit; p++)
                            {
                                radial += (s_1->f[p] * s_3->f[p] + Constant::AlphaSquared * s_1->g[p] * s_3->g[p])
                                          * Pot24[p] * dR[p];
                            }

                            TwoElectronIntegrals.insert(std::pair<unsigned int, double>(key, radial));
                        }

                        it_3.Next(); i3++;
                    }
                    it_1.Next(); i1++;
                }

                k+=2;
            }
            it_4.Next(); i4++;
        }
        it_2.Next(); i2++;
    }
}

void HamiltonianMatrix::GenerateMatrix()
{
    unsigned int i, j;

    M->Clear();
    M->WriteMode(true);

    // Loop through relativistic configurations
    RelativisticConfigList::const_iterator list_it = configs.begin();
    i=0;
    while(list_it != configs.end())
    {
        const ProjectionSet& proj_i = list_it->GetProjections();
        const std::vector< std::vector<double> >& coefficients_i = list_it->GetJCoefficients();

        RelativisticConfigList::const_iterator list_jt = list_it;
        j = i;
        
        while(list_jt != configs.end())
        {
            // Iterate over projections
            const ProjectionSet& proj_j = list_jt->GetProjections();
            const std::vector< std::vector<double> >& coefficients_j = list_jt->GetJCoefficients();

            ProjectionSet::const_iterator proj_it = proj_i.begin();
            ProjectionSet::const_iterator proj_jt = proj_j.begin();
            unsigned int pi = 0, pj;

            while(proj_it != proj_i.end())
            {
                proj_jt = proj_j.begin();
                pj = 0;
                while(proj_jt != proj_j.end())
                {
                    double operatorH = GetProjectionH(*proj_it, *proj_jt);

                    if(fabs(operatorH) > 1.e-16)
                    // Loop through JStates of the relativistic configurations and update M
                    for(unsigned int jstate_i = 0; jstate_i < coefficients_i.size(); jstate_i++)
                    {
                        unsigned int jstate_j_start = 0;
                        if(j == i)
                            jstate_j_start = jstate_i;

                        for(unsigned int jstate_j = jstate_j_start; jstate_j < coefficients_j.size(); jstate_j++)
                        {
                            double matrix_element = coefficients_i[jstate_i][pi] * coefficients_j[jstate_j][pj];
                            matrix_element = matrix_element * operatorH;

                            M->At(i + jstate_i, j + jstate_j) += matrix_element;
                        }
                    }
                    proj_jt++; pj++;
                }
                proj_it++; pi++;
            }

            list_jt++; j += coefficients_j.size();
        }

        list_it++; i += coefficients_i.size();
    }
}

double HamiltonianMatrix::GetProjectionH(const Projection& first, const Projection& second) const
{
    unsigned int diff[4];
    int numdiff = Projection::GetProjectionDifferences(first, second, diff);

    int sign;
    if(numdiff >= 0)
        sign = 1;
    else
        sign = -1;

    double value = 0.;

    if(numdiff == 0)
    {
        // Sum <i|f|i>
        for(unsigned int i=0; i<first.Size(); i++)
        {
            value += GetOneElectronIntegral(first[i].GetStateInfo(), first[i].GetStateInfo());
            
            // Sum(i < j) <ij|g|ij> - <ij|g|ji>
            for(unsigned int j=i+1; j<first.Size(); j++)
            {
                value += CoulombMatrixElement(first[i], first[j], first[i], first[j])
                         - CoulombMatrixElement(first[i], first[j], first[j], first[i]);
            }
        }
    }
    else if(abs(numdiff) == 1)
    {
        const ElectronInfo& f1 = first[diff[0]];
        const ElectronInfo& s1 = second[diff[1]];

        // a->b
        // <a|f|b>
        if((f1.M() == s1.M()) && (f1.Kappa() == s1.Kappa()))
            value = GetOneElectronIntegral(f1.GetStateInfo(), s1.GetStateInfo()) * sign;

        // Sum(e) <ae|g|be> - <ae|g|eb>
        for(unsigned int i=0; i<first.Size(); i++)
        {
            if(i != diff[0])
            {   const ElectronInfo& e(first[i]);
                value += sign * (CoulombMatrixElement(f1, e, s1, e) 
                                - CoulombMatrixElement(f1, e, e, s1));
            }
        }
    }
    else if(abs(numdiff) == 2)
    {
        const ElectronInfo& f1 = first[diff[0]];
        const ElectronInfo& s1 = second[diff[1]];
        const ElectronInfo& f2 = first[diff[2]];
        const ElectronInfo& s2 = second[diff[3]];

        // a->b, c->d
        // <ac|g|bd> - <ac|g|db>
        value = sign * (CoulombMatrixElement(f1, f2, s1, s2)
                        - CoulombMatrixElement(f1, f2, s2, s1));
    }

    return value;
}

double HamiltonianMatrix::CoulombMatrixElement(const ElectronInfo& e1, const ElectronInfo& e2, const ElectronInfo& e3, const ElectronInfo& e4) const
{
    if((e1.L() + e2.L() + e3.L() + e4.L())%2)
        return 0.;

    int two_q = e1.TwoM() - e3.TwoM();
    if(two_q != - e2.TwoM() + e4.TwoM())
        return 0.;

    unsigned int k = mmax(abs(int(e1.L()) - int(e3.L())), abs(int(e2.L()) - int(e4.L())));
    if((fabs(e1.J() - e3.J()) > double(k)) || (fabs(e2.J() - e4.J()) > double(k)))
        k += 2;

    unsigned int kmax = mmin(e1.L() + e3.L(), e2.L() + e4.L());
    if((e1.J() + e3.J() < double(kmax)) || (e2.J() + e4.J() < double(kmax)))
        kmax -= 2;

    if(k > kmax)
        return 0.;

    double q = double(two_q)/2.;

    const double* dR = states.GetLattice()->dR();

    double total = 0.;

    while(k <= kmax)
    {
        double coeff = 0.;
        if(fabs(q) <= k)
            coeff = Constant::Electron3j(e1.TwoJ(), e3.TwoJ(), k, -e1.TwoM(), e3.TwoM()) *
                    Constant::Electron3j(e2.TwoJ(), e4.TwoJ(), k, -e2.TwoM(), e4.TwoM());
            
        if(coeff)
            coeff = coeff * Constant::Electron3j(e1.TwoJ(), e3.TwoJ(), k, 1, -1) *
                            Constant::Electron3j(e2.TwoJ(), e4.TwoJ(), k, 1, -1);

        if(coeff)
        {
            if(int(q - e1.M() - e2.M() + 1.)%2)
                coeff = - coeff;

            coeff = coeff * sqrt(double(e1.MaxNumElectrons() * e2.MaxNumElectrons() *
                                        e3.MaxNumElectrons() * e4.MaxNumElectrons()));

            double radial = GetTwoElectronIntegral(k, e1.GetStateInfo(), e2.GetStateInfo(), e3.GetStateInfo(), e4.GetStateInfo());

            total += coeff * radial;
        }

        k = k+2;
    }

    return total;
}

double HamiltonianMatrix::GetProjectionSMS(const Projection& first, const Projection& second) const
{
    unsigned int diff[4];
    int numdiff = Projection::GetProjectionDifferences(first, second, diff);

    int sign;
    if(numdiff >= 0)
        sign = 1;
    else
        sign = -1;

    double value = 0.;

    if(numdiff == 0)
    {
        // Sum <i|f|i>
        for(unsigned int i=0; i<first.Size(); i++)
        {
//            value += CoreSMS(p1[i], p1[i]);
            
            // Sum(i < j) <ij|g|ij> - <ij|g|ji>
            for(unsigned int j=i+1; j<first.Size(); j++)
            {
                value += SMSMatrixElement(first[i], first[j], first[i], first[j])
                         - SMSMatrixElement(first[i], first[j], first[j], first[i]);
            }
        }
    }
    else if(abs(numdiff) == 1)
    {
        const ElectronInfo& f1 = first[diff[0]];
        const ElectronInfo& s1 = second[diff[1]];
        // a->b
        // <a|f|b>
//        if((p1[0].M() == p2[0].M()) && (p1[0].Kappa() == p2[0].Kappa()))
//            value = CoreSMS(p1[0], p2[0]) * sign;

        // Sum(e) <ae|g|be> - <ae|g|eb>
        for(unsigned int i=0; i<first.Size(); i++)
        {   
            if(i != diff[0])
            {   const ElectronInfo& e(first[i]);
                value += sign * (SMSMatrixElement(f1, e, s1, e) 
                                - SMSMatrixElement(f1, e, e, s1));
            }
        }
    }
    else if(abs(numdiff) == 2)
    {
        const ElectronInfo& f1 = first[diff[0]];
        const ElectronInfo& s1 = second[diff[1]];
        const ElectronInfo& f2 = first[diff[2]];
        const ElectronInfo& s2 = second[diff[3]];

        // a->b, c->d
        // <ac|g|bd> - <ac|g|db>
        value = sign * (SMSMatrixElement(f1, f2, s1, s2) 
                        - SMSMatrixElement(f1, f2, s2, s1));
    }

    return value;
}

double HamiltonianMatrix::CoreSMS(const ElectronInfo& e1, const ElectronInfo& e2) const
{
    if((!states.GetCore()->GetNuclearInverseMass()) || (e1.Kappa() != e2.Kappa()))
        return 0.;

    double E = 0.;

    const State& state1 = *states.GetState(e1.GetStateInfo());
    const State& state2 = *states.GetState(e2.GetStateInfo());

    const Core* core = states.GetCore();

    // Sum over all core states
    CoulombIntegrator I(*states.GetLattice());
    ConstDiscreteStateIterator cs = core->GetConstDiscreteStateIterator();
    while(!cs.AtEnd())
    {
        const DiscreteState& other = *(cs.GetState());

        // k == 1
        if((other.L() + state1.L())%2 == 1)
        {
            const unsigned int k = 1;

            double coefficient = Constant::Wigner3j(k, state1.J(), other.J(), 0., .5, -.5);
            coefficient = -(2. * abs(other.Kappa())) * coefficient * coefficient;

            E = E + coefficient * core->GetNuclearInverseMass() * I.IsotopeShiftIntegral(state1, other) * I.IsotopeShiftIntegral(other, state2);
        }
        cs.Next();
    }

    return E;
}

double HamiltonianMatrix::SMSMatrixElement(const ElectronInfo& e1, const ElectronInfo& e2, const ElectronInfo& e3, const ElectronInfo& e4) const
{
    if((e1.L() + e2.L() + e3.L() + e4.L())%2)
        return 0.;

    int two_q = e1.TwoM() - e3.TwoM();
    if(two_q != - e2.TwoM() + e4.TwoM())
        return 0.;

    unsigned int k = mmax(abs(int(e1.L()) - int(e3.L())), abs(int(e2.L()) - int(e4.L())));
    if((fabs(e1.J() - e3.J()) > double(k)) || (fabs(e2.J() - e4.J()) > double(k)))
        k += 2;

    unsigned int kmax = mmin(e1.L() + e3.L(), e2.L() + e4.L());
    if((e1.J() + e3.J() < double(kmax)) || (e2.J() + e4.J() < double(kmax)))
        kmax -= 2;

    if((k > kmax) || (k != 1))
        return 0.;

    double q = double(two_q)/2.;

    const double* dR = states.GetLattice()->dR();

    const State* s1 = states.GetState(e1.GetStateInfo());
    const State* s2 = states.GetState(e2.GetStateInfo());
    const State* s3 = states.GetState(e3.GetStateInfo());
    const State* s4 = states.GetState(e4.GetStateInfo());

    CoulombIntegrator I(*states.GetLattice());

    double total = 0.;

    // k == 1 only
    double coeff = 0.;
    if(fabs(q) <= k)
        coeff = Constant::Electron3j(e1.TwoJ(), e3.TwoJ(), k, -e1.TwoM(), e3.TwoM()) *
                Constant::Electron3j(e2.TwoJ(), e4.TwoJ(), k, -e2.TwoM(), e4.TwoM());
        
    if(coeff)
        coeff = coeff * Constant::Electron3j(e1.TwoJ(), e3.TwoJ(), k, 1, -1) *
                        Constant::Electron3j(e2.TwoJ(), e4.TwoJ(), k, 1, -1);

    if(coeff)
    {
        if(int(q - e1.M() - e2.M() + 1)%2)
            coeff = - coeff;

        coeff = coeff * sqrt(double(e1.MaxNumElectrons() * e2.MaxNumElectrons() *
                                    e3.MaxNumElectrons() * e4.MaxNumElectrons()));

        // Specific Mass Shift
        double SMS = states.GetCore()->GetNuclearInverseMass();
        if(SMS)
            SMS = SMS * I.IsotopeShiftIntegral(*s1, *s3) * I.IsotopeShiftIntegral(*s2, *s4);

        total = - coeff * SMS;
    }

    return total;
}

void HamiltonianMatrix::PollMatrix()
{
    M->WriteMode(false);

    unsigned int i, j, count;
    double value;
    unsigned int range[10];
    for(i=0; i<10; i++)
        range[i] = 0;

    for(i=0; i<N; i++)
        for(j=0; j<N; j++)
        {
            value = fabs(M->At(i, j));
            if(value >= 1.)
                range[9]++;
            else
            {   count = 0;
                while(value > 1.e-16)
                {   value = value/100.;
                    count++;
                }
                range[count]++;
            }
        }

    for(i=0; i<10; i++)
        printf("%d %d %.2f\n", i, range[i], double(range[i])/double(N*N)*100.);
}

void HamiltonianMatrix::SolveMatrix(unsigned int num_solutions, unsigned int two_j, bool gFactors)
{
    M->WriteMode(false);

    std::cout << "\nFinding solutions" << std::endl;

    if(NumSolutions)
    {   delete[] V;
        delete[] E;
    }
    NumSolutions = mmin(num_solutions, N);

    V = new double[NumSolutions * N];
    E = new double[NumSolutions];

    Eigensolver solver;
    solver.SolveLargeSymmetric(M, E, V, N, NumSolutions);

    // Calculate g-Factors
    double* g_factors;
    if(gFactors)
    {   g_factors = new double[NumSolutions];
        GetgFactors(two_j, g_factors);
    }

    unsigned int i, j;

    std::cout << "Solutions for J = " << double(two_j)/2. << ": " << std::endl;
    for(i=0; i<NumSolutions; i++)
    {
        unsigned int solution = i;

        printf("%d: %.7f    %.5f /cm\n", i, E[solution], E[solution]*Constant::HartreeEnergy_cm);

        // Get non-rel configuration percentages
        RelativisticConfigList::const_iterator list_it = configs.begin();
        std::map<Configuration, double> percentages;  // Map non-rel configurations to percentages

        j = 0;
        while(list_it != configs.end())
        {
            Configuration nrconfig(list_it->GetNonRelConfiguration());
            if(percentages.find(nrconfig) == percentages.end())
                percentages[nrconfig] = 0.;

            const std::vector< std::vector<double> >& coefficients = list_it->GetJCoefficients();
            for(unsigned int Jstate = 0; Jstate < coefficients.size(); Jstate++)
            {
                double coeff = V[solution*N + j];
                coeff = coeff * coeff * 100;

                percentages[nrconfig] += coeff;
                j++;
            }

            list_it++;
        }

        std::map<Configuration, double>::const_iterator it = percentages.begin();
        while(it != percentages.end())
        {
            if(it->second > 3.)
                printf("    %s\t%.1f%%\n", it->first.Name().c_str(), it->second);
            it++;
        }

        if(gFactors)
            printf("    g-factor = %.4f\n", g_factors[solution]);

        printf("\n");
    }
    
    if(gFactors)
        delete[] g_factors;
}

void HamiltonianMatrix::GetEigenvalues() const
{
    double* total = new double[NumSolutions];
    double* coeff = new double[NumSolutions];
    unsigned int solution;
    for(solution = 0; solution < NumSolutions; solution++)
        total[solution] = 0.;

    // Iterate over different relativistic configurations
    unsigned int i=0, j;
    RelativisticConfigList::const_iterator list_it = configs.begin();
    while(list_it != configs.end())
    {
        const ProjectionSet& proj_i = list_it->GetProjections();
        const std::vector< std::vector<double> >& coefficients_i = list_it->GetJCoefficients();

        RelativisticConfigList::const_iterator list_jt = configs.begin();
        j = 0;
        while(list_jt != configs.end())
        {
            const ProjectionSet& proj_j = list_jt->GetProjections();
            const std::vector< std::vector<double> >& coefficients_j = list_jt->GetJCoefficients();

            // Iterate over projections
            ProjectionSet::const_iterator pi_it = proj_i.begin();
            unsigned int pi = 0;
            while(pi_it != proj_i.end())
            {
                ProjectionSet::const_iterator pj_it = proj_j.begin();
                unsigned int pj = 0;
                while(pj_it != proj_j.end())
                {
                    double matrix_element = GetProjectionSMS(*pi_it, *pj_it);

                    if(matrix_element)
                    {   // Summation over jstates
                        for(solution = 0; solution < NumSolutions; solution++)
                            coeff[solution] = 0.;

                        for(unsigned int jstate_i = 0; jstate_i < coefficients_i.size(); jstate_i++)
                        {
                            for(unsigned int jstate_j = 0; jstate_j < coefficients_j.size(); jstate_j++)
                            {
                                for(solution = 0; solution < NumSolutions; solution++)
                                {
                                    coeff[solution] += coefficients_i[jstate_i][pi] * coefficients_j[jstate_j][pj]
                                                    * V[solution*N + i + jstate_i]
                                                    * V[solution*N + j + jstate_j];
                                }
                            }
                        }

                        // If the relativistic configs are different, count twice
                        if(i != j)
                        {   for(solution = 0; solution < NumSolutions; solution++)
                                coeff[solution] = coeff[solution] * 2.;
                        }

                        for(solution = 0; solution < NumSolutions; solution++)
                            total[solution] += coeff[solution] * matrix_element;
                    }

                    pj_it++; pj++;
                }

                pi_it++; pi++;
            }

            if(list_jt == list_it)
                list_jt = configs.end();
            else
                list_jt++; j+=coefficients_j.size();
        }

        list_it++; i+=coefficients_i.size();
    }

    for(solution = 0; solution < NumSolutions; solution++)
        printf("%d: %.7f    %.5f /cm\n", solution, total[solution], total[solution]*Constant::HartreeEnergy_cm);

    delete[] total;
    delete[] coeff;
}

void HamiltonianMatrix::GetgFactors(unsigned int two_j, double* g_factors) const
{
    // Case where J=0
    if(two_j == 0)
    {   for(unsigned int solution = 0; solution < NumSolutions; solution++)
            g_factors[solution] = 0.;
        return;
    }

    double* total = new double[NumSolutions];
    double* coeff = new double[NumSolutions];
    unsigned int solution;
    for(solution = 0; solution < NumSolutions; solution++)
        total[solution] = 0.;

    unsigned int diff[4];   // Storage for projection differences.
    unsigned int num_electrons = configs.front().NumParticles();

    // Iterate over different relativistic configurations
    unsigned int i=0, j;
    RelativisticConfigList::const_iterator list_it = configs.begin();
    while(list_it != configs.end())
    {
        const ProjectionSet& proj_i = list_it->GetProjections();
        const std::vector< std::vector<double> >& coefficients_i = list_it->GetJCoefficients();

        RelativisticConfigList::const_iterator list_jt = configs.begin();
        j = 0;
        while(list_jt != configs.end())
        {
            const ProjectionSet& proj_j = list_jt->GetProjections();
            const std::vector< std::vector<double> >& coefficients_j = list_jt->GetJCoefficients();

            // Iterate over projections
            ProjectionSet::const_iterator pi_it = proj_i.begin();
            unsigned int pi = 0;
            while(pi_it != proj_i.end())
            {
                ProjectionSet::const_iterator pj_it = proj_j.begin();
                unsigned int pj = 0;
                while(pj_it != proj_j.end())
                {
                    // <pi| Jz | pj>
                    double matrix_element = 0.;
                    int num_diff = Projection::GetProjectionDifferences(*pi_it, *pj_it, diff);
                    if(num_diff == 0)
                    {
                        for(unsigned int i=0; i<num_electrons; i++)
                            matrix_element += GetSz((*pi_it)[i]);
                    }
                    else if(abs(num_diff) == 1)
                    {
                        matrix_element = GetSz((*pi_it)[diff[0]], (*pj_it)[diff[1]]);
                        if(num_diff == -1)
                            matrix_element = -matrix_element;
                    }

                    // coefficients
                    if(matrix_element)
                    {
                        // Summation over jstates
                        for(solution = 0; solution < NumSolutions; solution++)
                            coeff[solution] = 0.;

                        for(unsigned int jstate_i = 0; jstate_i < coefficients_i.size(); jstate_i++)
                        {
                            for(unsigned int jstate_j = 0; jstate_j < coefficients_j.size(); jstate_j++)
                            {
                                for(solution = 0; solution < NumSolutions; solution++)
                                {
                                    coeff[solution] += coefficients_i[jstate_i][pi] * coefficients_j[jstate_j][pj]
                                                    * V[solution*N + i + jstate_i]
                                                    * V[solution*N + j + jstate_j];
                                }
                            }
                        }

                        // If the relativistic configs are different, count twice
                        if(i != j)
                        {   for(solution = 0; solution < NumSolutions; solution++)
                                coeff[solution] = coeff[solution] * 2.;
                        }

                        for(solution = 0; solution < NumSolutions; solution++)
                            total[solution] += coeff[solution] * matrix_element;
                    }

                    pj_it++; pj++;
                }

                pi_it++; pi++;
            }

            if(list_jt == list_it)
                list_jt = configs.end();
            else
                list_jt++; j+=coefficients_j.size();
        }

        list_it++; i+=coefficients_i.size();
    }

    double J = two_j/2.;
    for(solution = 0; solution < NumSolutions; solution++)
        g_factors[solution] = total[solution]/J + 1.;

    delete[] coeff;
    delete[] total;
}

double HamiltonianMatrix::GetSz(const ElectronInfo& e) const
{
    if(e.Kappa() > 0)   // J = L - 1/2
        return - e.M()/(2.*e.L() + 1.);
    else
        return e.M()/(2.*e.L() + 1.);
}

double HamiltonianMatrix::GetSz(const ElectronInfo& e1, const ElectronInfo& e2) const
{
    if((e1.L() == e2.L()) &&
       (e1.TwoM() == e2.TwoM()))
    {
        double overlap = 0.;
        const State* p1 = states.GetState(e1.GetStateInfo());
        const State* p2 = states.GetState(e2.GetStateInfo());
        const double* dR = states.GetLattice()->dR();
        for(int i=0; i<p1->Size(); i++)
            overlap += (p1->f[i] * p2->f[i] + Constant::AlphaSquared * p1->g[i] * p2->g[i]) * dR[i];

        double Lplushalf = double(e1.L()) + 0.5;
        double M = e1.M();
        return -sqrt((Lplushalf + M)*(Lplushalf - M))/(2.*Lplushalf) * overlap;
    }
    else
        return 0.;
}

inline void swap(unsigned int& i1, unsigned int& i2)
{   unsigned int temp = i1;
    i1 = i2;
    i2 = temp;
}

double HamiltonianMatrix::GetOneElectronIntegral(const StateInfo& s1, const StateInfo& s2) const
{
    unsigned int i1 = state_index.find(s1)->second;
    unsigned int i2 = state_index.find(s2)->second;

    if(i1 <= i2)
        return OneElectronIntegrals.find(i1 * NumStates + i2)->second;
    else
        return OneElectronIntegrals.find(i2 * NumStates + i1)->second;
}

double HamiltonianMatrix::GetTwoElectronIntegral(unsigned int k, const StateInfo& s1, const StateInfo& s2, const StateInfo& s3, const StateInfo& s4) const
{
    unsigned int i1 = state_index.find(s1)->second;
    unsigned int i2 = state_index.find(s2)->second;
    unsigned int i3 = state_index.find(s3)->second;
    unsigned int i4 = state_index.find(s4)->second;

    bool sms_sign = true;

    if(i3 < i1)
    {   swap(i3, i1);
        sms_sign = !sms_sign;
    }
    if(i4 < i2)
    {   swap(i4, i2);
        sms_sign = !sms_sign;
    }
    if(i2 < i1)
    {   swap(i2, i1);
        swap(i3, i4);
    }
    if((i1 == i2) && (i4 < i3))
        swap(i3, i4);

    unsigned int key = k  * NumStates*NumStates*NumStates*NumStates +
                       i1 * NumStates*NumStates*NumStates + 
                       i2 * NumStates*NumStates + 
                       i3 * NumStates +
                       i4;

    double radial = 0.;
    if(TwoElectronIntegrals.find(key) != TwoElectronIntegrals.end())
    {   
        radial = TwoElectronIntegrals.find(key)->second;
        
        if(include_sms_v2 && (k == 1))
        {   double SMS = states.GetCore()->GetNuclearInverseMass();
            if(SMS)
            {   SMS = SMS * SMSIntegrals.find(i1*NumStates + i3)->second * SMSIntegrals.find(i2*NumStates + i4)->second;
                if(!sms_sign)
                    SMS = -SMS;
                radial = radial - SMS;
            }
        }
    }
    else
        printf("Shit\n");

    return radial;
}

/*
double ConfigGenerator::gFactor(const RelativisticConfigList& rlist, const double* jstate_coefficients, double J)
{
    double C = Constant::Wigner3j(J, 1., J, -J, 0., J);
    if(C == 0.)
        return 0.;

    // Get list of coefficients of projections
    unsigned int num_projections = 0;
    RelativisticConfigList::const_iterator it = rlist.begin();
    while(it != rlist.end())
    {   num_projections += it->GetProjections().size();
        it++;
    }

    std::vector<double> projection_coefficients(num_projections);
    std::vector<Projection> projections(num_projections);
    unsigned int jstate = 0;
    unsigned int proj_num = 0;
    it = rlist.begin();
    while(it != rlist.end())
    {
        const std::vector< std::vector<double> >& coefficients = it->GetJCoefficients();
        const ProjectionSet& proj = it->GetProjections();
        ProjectionSet::const_iterator proj_it = proj.begin();
        unsigned int i=0;   // index of projection, tied to proj_it
        while(proj_it != proj.end())
        {
            projections[proj_num] = *proj_it;

            for(unsigned int j=0; j<coefficients.size(); j++)
            {
                projection_coefficients[proj_num] += jstate_coefficients[jstate + j] * coefficients[j][i];
            }

            proj_it++;
            i++;
            proj_num++;
        }
        jstate += coefficients.size();
        it++;
    }

    unsigned int i, j;

    // Make ro matrix
    unsigned int e_size = ElectronSet.size();
    double* ro = new double[e_size * e_size];
    for(i=0; i<e_size * e_size; i++)
        ro[i] = 0.;

    // Loop through all projections to get electron coefficients
    for(i=0; i<num_projections; i++)
    {
        for(j=i; j<num_projections; j++)
        {
            // Get differences between projections and update ro.
            // Make copies because they get munged by GetProjectionDifferences.
            Projection proj_i(projections[i]);
            Projection proj_j(projections[j]);

            unsigned int diff[4];
            int numdiff = Projection::GetProjectionDifferences(proj_i, proj_j, diff);
            if(numdiff == 0)
            {   // All electrons
                for(unsigned int e=0; e<proj_i.Size(); e++)
                {
                    double value = projection_coefficients[i] * projection_coefficients[j];
                    if(i != j)
                        value = value * 2.;

                    unsigned int electron_index = ElectronSet[proj_i[e]];
                    ro[electron_index * e_size + electron_index] += value;
                }
            }
            else if(abs(numdiff) == 1)
            {
                double value = projection_coefficients[i] * projection_coefficients[j];
                if(numdiff < 0)
                    value = -value;

                ro[ElectronSet[proj_i[diff[0]]] * e_size + ElectronSet[proj_j[diff[1]]]] += value;
                ro[ElectronSet[proj_j[diff[1]]] * e_size + ElectronSet[proj_i[diff[0]]]] += value;
            }
        }
    }

    // Reduce ro matrix
    double g = 0.;

    std::map<ElectronInfo, unsigned int>::const_iterator eit, ejt;
    eit = ElectronSet.begin();
    while(eit != ElectronSet.end())
    {
        ejt = ElectronSet.begin();
        while(ejt != ElectronSet.end())
        {
            const ElectronInfo& i_info = eit->first;
            const ElectronInfo& j_info = ejt->first;

            if((i_info.PQN() == j_info.PQN()) &&
               (i_info.L() == j_info.L()) &&
               (i_info.TwoM() == j_info.TwoM()))
            {
                double x = ro[eit->second * e_size + ejt->second];
                double value = Constant::Wigner3j(j_info.J(), 1., i_info.J(), -i_info.M(), 0., i_info.M());
                if(x && value)
                {
                    value = value * pow(-1., (j_info.TwoJ() - j_info.TwoM())/2.);

                    if(i_info.TwoJ() == j_info.TwoJ())
                        value = value * (i_info.J() + 0.5) * 
                                sqrt((2.*i_info.J() + 1.) * (4.*i_info.J() - 2.*i_info.L() + 1.) / (2.*i_info.L() + 1.));
                    else
                        value = value * pow(-1., j_info.J() + 0.5 + j_info.L()) *
                                sqrt(2. * j_info.L() * (j_info.L() + 1.) / (2. * j_info.L() + 1.));

                    g += value * x / C;
                }
            }

            ejt++;
        }
        eit++;
    }

    delete[] ro;
    return g/sqrt(J * (J + 1.) * (2.*J + 1.));
}
*/
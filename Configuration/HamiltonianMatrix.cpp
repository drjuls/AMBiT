#include "Include.h"
#include "HamiltonianMatrix.h"
#include "Universal/SmallMatrix.h"
#include "Universal/SymMatrix.h"
#include "HartreeFock/State.h"
#include "Universal/Eigensolver.h"
#include "Universal/Constant.h"

#define SMALL_MATRIX_LIM 2000

// Include this define for the box diagrams of "wrong" parity.
//#define INCLUDE_EXTRA_BOX_DIAGRAMS

// Include this define to only include sigma3 when both states are leading configurations
// (instead of just one).
#define SIGMA3_AND

HamiltonianMatrix::HamiltonianMatrix(const CIIntegrals& coulomb_integrals, const RelativisticConfigList& rconfigs):
    integrals(coulomb_integrals), configs(rconfigs), NumSolutions(0), M(NULL), include_sigma3(false)
{
    // Set up matrix
    N = 0;
    RelativisticConfigList::const_iterator it = configs.begin();
    while(it != configs.end())
    {   N += it->NumJStates();
        it++;
    }

    *logstream << " " << N << " " << std::flush;
    *outstream << " Number of J-configurations = " << N << std::endl;
}

void HamiltonianMatrix::GenerateMatrix()
{
    unsigned int i, j;

    if(M == NULL)
    {   if(N <= SMALL_MATRIX_LIM)
            M = new SmallMatrix(N);
        else
            M = new SymMatrix(N);
    }
    else
        M->Clear();

    M->WriteMode(true);

    // Set valence energies
    if(include_sigma3)
        sigma3calc->UseBrillouinWignerPT();

    // Loop through relativistic configurations
    RelativisticConfigList::const_iterator list_it = configs.begin();
    i=0;
    while(list_it != configs.end())
    {
        const ProjectionSet& proj_i = list_it->GetProjections();
        unsigned int proj_i_size = proj_i.size();
        unsigned int num_states_i = list_it->NumJStates();
        const double* coefficients_i = list_it->GetJCoefficients();

        RelativisticConfigList::const_iterator list_jt = list_it;
        j = i;
        
        while(list_jt != configs.end())
        {
            // Iterate over projections
            const ProjectionSet& proj_j = list_jt->GetProjections();
            unsigned int proj_j_size = proj_j.size();
            unsigned int num_states_j = list_jt->NumJStates();
            const double* coefficients_j = list_jt->GetJCoefficients();

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
                    for(unsigned int jstate_i = 0; jstate_i < num_states_i; jstate_i++)
                    {
                        unsigned int jstate_j_start = 0;
                        if(j == i)
                            jstate_j_start = jstate_i;

                        for(unsigned int jstate_j = jstate_j_start; jstate_j < num_states_j; jstate_j++)
                        {
                            double matrix_element = coefficients_i[jstate_i*proj_i_size + pi]
                                                   * coefficients_j[jstate_j*proj_j_size + pj];
                            matrix_element = matrix_element * operatorH;

                            M->At(i + jstate_i, j + jstate_j) += matrix_element;
                        }
                    }
                    proj_jt++; pj++;
                }
                proj_it++; pi++;
            }

            list_jt++; j += num_states_j;
        }

        list_it++; i += num_states_i;
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
            value += integrals.GetOneElectronIntegral(first[i], first[i]);
            
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
            value = integrals.GetOneElectronIntegral(f1, s1) * sign;

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

    if(include_sigma3)
        value += GetSigma3(first, second);

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

    double q = double(two_q)/2.;

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

            double radial = integrals.GetTwoElectronIntegral(k, e1, e2, e3, e4);

            total += coeff * radial;
        }

        k = k+2;
    }

#ifdef INCLUDE_EXTRA_BOX_DIAGRAMS
    // Include the box diagrams with "wrong" parity.
    k = (unsigned int)mmax(fabs(e1.J() - e3.J()), fabs(e2.J() - e4.J()));
    if((k + e1.L() + e3.L())%2 == 0)
        k++;

    kmax = (unsigned int)mmin(e1.J() + e3.J(), e2.J() + e4.J());

    while(k <= kmax)
    {
        double radial = integrals.GetTwoElectronIntegral(k, e1, e2, e3, e4);

        if(radial)
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

                total += coeff * radial;
            }
        }

        k = k+2;
    }
#endif

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
        double SMS = integrals.GetNuclearInverseMass();
        if(SMS)
            SMS = SMS * integrals.GetSMSIntegral(e1, e3) * integrals.GetSMSIntegral(e2, e4);

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
        *outstream << i << " " << range[i] << " " << double(range[i])/double(N*N)*100. << std::endl;
}

void HamiltonianMatrix::SolveMatrix(unsigned int num_solutions, unsigned int two_j, bool gFactors)
{
    M->WriteMode(false);

    *outstream << "\nFinding solutions" << std::endl;

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

    *outstream << "Solutions for J = " << double(two_j)/2. << ": " << std::endl;
    for(i=0; i<NumSolutions; i++)
    {
        unsigned int solution = i;

        *outstream << i << ": " << std::setprecision(8) << E[solution] << "    "
            << std::setprecision(12) << E[solution]*Constant::HartreeEnergy_cm << " /cm" << std::endl;

        // Get non-rel configuration percentages
        RelativisticConfigList::const_iterator list_it = configs.begin();
        std::map<Configuration, double> percentages;  // Map non-rel configurations to percentages

        j = 0;
        while(list_it != configs.end())
        {
            Configuration nrconfig(list_it->GetNonRelConfiguration());
            if(percentages.find(nrconfig) == percentages.end())
                percentages[nrconfig] = 0.;

            for(unsigned int Jstate = 0; Jstate < list_it->NumJStates(); Jstate++)
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
            if(it->second > 1.)
                *outstream << std::setw(20) << it->first.Name() << "  "<< std::setprecision(2)
                    << it->second << "%" << std::endl;
            it++;
        }

        if(gFactors)
            *outstream << "    g-factor = " << std::setprecision(5) << g_factors[solution] << std::endl;

        *outstream << std::endl;
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
        unsigned int proj_i_size = proj_i.size();
        unsigned int num_states_i = list_it->NumJStates();
        const double* coefficients_i = list_it->GetJCoefficients();

        RelativisticConfigList::const_iterator list_jt = list_it;
        j = i;
        while(list_jt != configs.end())
        {
            const ProjectionSet& proj_j = list_jt->GetProjections();
            unsigned int proj_j_size = proj_j.size();
            unsigned int num_states_j = list_jt->NumJStates();
            const double* coefficients_j = list_jt->GetJCoefficients();

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

                        for(unsigned int jstate_i = 0; jstate_i < num_states_i; jstate_i++)
                        {
                            for(unsigned int jstate_j = 0; jstate_j < num_states_j; jstate_j++)
                            {
                                for(solution = 0; solution < NumSolutions; solution++)
                                {
                                    coeff[solution] += coefficients_i[jstate_i*proj_i_size + pi]
                                                    * coefficients_j[jstate_j*proj_j_size + pj]
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

            list_jt++; j+=num_states_j;
        }

        list_it++; i+=num_states_i;
    }

    for(solution = 0; solution < NumSolutions; solution++)
        *outstream << solution << ": " 
            << std::setprecision(8) << total[solution] << "    "
            << total[solution]*Constant::HartreeEnergy_cm
            << " /cm" << std::endl;

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
        unsigned int proj_i_size = proj_i.size();
        unsigned int num_states_i = list_it->NumJStates();
        const double* coefficients_i = list_it->GetJCoefficients();

        RelativisticConfigList::const_iterator list_jt = list_it;
        j = i;
        while(list_jt != configs.end())
        {
            const ProjectionSet& proj_j = list_jt->GetProjections();
            unsigned int proj_j_size = proj_j.size();
            unsigned int num_states_j = list_jt->NumJStates();
            const double* coefficients_j = list_jt->GetJCoefficients();

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

                        for(unsigned int jstate_i = 0; jstate_i < num_states_i; jstate_i++)
                        {
                            for(unsigned int jstate_j = 0; jstate_j < num_states_j; jstate_j++)
                            {
                                for(solution = 0; solution < NumSolutions; solution++)
                                {
                                    coeff[solution] += coefficients_i[jstate_i*proj_i_size + pi]
                                                    * coefficients_j[jstate_j*proj_j_size + pj]
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

            list_jt++; j+=num_states_j;
        }

        list_it++; i+=num_states_i;
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
        double overlap = integrals.GetOverlapIntegral(e1, e2);

        double Lplushalf = double(e1.L()) + 0.5;
        double M = e1.M();
        return -sqrt((Lplushalf + M)*(Lplushalf - M))/(2.*Lplushalf) * overlap;
    }
    else
        return 0.;
}

void HamiltonianMatrix::AddLeadingConfigurations(const ConfigList& list)
{
    ConfigList::const_iterator it = list.begin();
    while(it != list.end())
    {   leading_configs.insert(*it);
        it++;
    }
}

void HamiltonianMatrix::AddLeadingConfigurations(const Configuration& config)
{
    leading_configs.insert(config);
}

double HamiltonianMatrix::GetSigma3(const Projection& first, const Projection& second) const
{
#ifdef SIGMA3_AND
    // Check that first AND second are leading configurations
    if((leading_configs.find(first.GetNonRelConfiguration()) == leading_configs.end()) ||
       (leading_configs.find(second.GetNonRelConfiguration()) == leading_configs.end()))
        return 0.;
#else
    // Check that first OR second is a leading configuration
    if((leading_configs.find(first.GetNonRelConfiguration()) == leading_configs.end()) &&
       (leading_configs.find(second.GetNonRelConfiguration()) == leading_configs.end()))
        return 0.;
#endif

    unsigned int diff[6];
    int numdiff = Projection::GetProjectionDifferences3(first, second, diff);

    int sign;
    if(numdiff >= 0)
        sign = 1;
    else
        sign = -1;

    double value = 0.;

    if(numdiff == 0)
    {
        // Sum(i < j < k) Sigma3(ijk, ijk)
        for(unsigned int i=0; i<first.Size(); i++)
        {
            for(unsigned int j=i+1; j<first.Size(); j++)
            {
                for(unsigned int k=j+1; k<first.Size(); k++)
                {
                    value += Sigma3(first[i], first[j], first[k], first[i], first[j], first[k]);
                }
            }
        }
    }
    else if(abs(numdiff) == 1)
    {
        const ElectronInfo& f1 = first[diff[0]];
        const ElectronInfo& s1 = second[diff[1]];

        // Sum(i < j) Sigma3(aij, bij)
        for(unsigned int i=0; i<first.Size(); i++)
        {
            for(unsigned int j=i+1; j<first.Size(); j++)
            {
                if((i != diff[0]) && (j != diff[0]))
                {   
                    value += sign * Sigma3(f1, first[i], first[j], s1, first[i], first[j]);
                }
            }
        }
    }
    else if(abs(numdiff) == 2)
    {
        const ElectronInfo& f1 = first[diff[0]];
        const ElectronInfo& s1 = second[diff[1]];
        const ElectronInfo& f2 = first[diff[2]];
        const ElectronInfo& s2 = second[diff[3]];

        // Sum(i) Sigma3(abi, cdi)
        for(unsigned int i=0; i<first.Size(); i++)
        {
            if((i != diff[0]) && (i != diff[2]))
               value += sign * Sigma3(f1, f2, first[i], s1, s2, first[i]);
        }
    }
    else if(abs(numdiff) == 3)
    {
        const ElectronInfo& f1 = first[diff[0]];
        const ElectronInfo& s1 = second[diff[1]];
        const ElectronInfo& f2 = first[diff[2]];
        const ElectronInfo& s2 = second[diff[3]];
        const ElectronInfo& f3 = first[diff[4]];
        const ElectronInfo& s3 = second[diff[5]];

        // Sigma3(abc, def)
        value = sign * Sigma3(f1, f2, f3, s1, s2, s3);
    }

    return value;
}

double HamiltonianMatrix::Sigma3(const ElectronInfo& e1, const ElectronInfo& e2, const ElectronInfo& e3,
           const ElectronInfo& e4, const ElectronInfo& e5, const ElectronInfo& e6) const
{
    // Check momentum projections
    if(e1.TwoM() + e2.TwoM() + e3.TwoM() != e4.TwoM() + e5.TwoM() + e6.TwoM())
        return 0.;

    // Check parity
    if((e1.L() + e2.L() + e4.L())%2 != (e3.L() + e5.L() + e6.L())%2)
        return 0.;

    double value = 0.;

    // The sign changes for odd permutations
    value =   Sigma3LinePermutations(e1, e2, e3, e4, e5, e6)
            + Sigma3LinePermutations(e1, e2, e3, e5, e6, e4)
            + Sigma3LinePermutations(e1, e2, e3, e6, e4, e5)
            - Sigma3LinePermutations(e1, e2, e3, e5, e4, e6)
            - Sigma3LinePermutations(e1, e2, e3, e6, e5, e4)
            - Sigma3LinePermutations(e1, e2, e3, e4, e6, e5);

    return value;
}

/** This function does the line permutations, putting the pairs on different levels
    of the three-body interaction.
    */
inline double HamiltonianMatrix::Sigma3LinePermutations(const ElectronInfo& e1, const ElectronInfo& e2, const ElectronInfo& e3,
           const ElectronInfo& e4, const ElectronInfo& e5, const ElectronInfo& e6) const
{
    double value = 0.;

    // There are no sign changes, since there are the same number
    // of permutations on both sides
    value =   sigma3calc->GetSecondOrderSigma3(e1, e2, e3, e4, e5, e6)
            + sigma3calc->GetSecondOrderSigma3(e2, e3, e1, e5, e6, e4)
            + sigma3calc->GetSecondOrderSigma3(e3, e1, e2, e6, e4, e5)
            + sigma3calc->GetSecondOrderSigma3(e3, e2, e1, e6, e5, e4)
            + sigma3calc->GetSecondOrderSigma3(e2, e1, e3, e5, e4, e6)
            + sigma3calc->GetSecondOrderSigma3(e1, e3, e2, e4, e6, e5);

    return value;
}

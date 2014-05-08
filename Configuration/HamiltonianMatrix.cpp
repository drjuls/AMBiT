#include "Include.h"
#include "HamiltonianMatrix.h"
#include "Universal/SmallMatrix.h"
#include "Universal/SymMatrix.h"
#include "HartreeFock/SingleParticleWavefunction.h"
#include "Universal/Eigensolver.h"
#include "Universal/MathConstant.h"
//#include "ConfigFileGenerator.h"

#define SMALL_MATRIX_LIM 2000

// Include this define for the box diagrams of "wrong" parity.
//#define INCLUDE_EXTRA_BOX_DIAGRAMS

// Include this define to only include sigma3 when both states are leading configurations
// (instead of just one).
//#define SIGMA3_AND

class HamiltonianOperator
{
public:
    HamiltonianOperator(const CIIntegrals& ci_integrals): integrals(ci_integrals) {}

    inline double GetMatrixElement(const ElectronInfo& e1, const ElectronInfo& e2) const
    {
        if(e1.TwoM() == e2.TwoM() && e1.Kappa() == e2.Kappa())
            return integrals.GetOneElectronIntegral(e1, e2);
        else
            return 0.;
    }

    double GetMatrixElement(const ElectronInfo& e1, const ElectronInfo& e2, const ElectronInfo& e3, const ElectronInfo& e4) const;

protected:
    const CIIntegrals& integrals;
};

HamiltonianMatrix::HamiltonianMatrix(const CIIntegrals* coulomb_integrals, pRelativisticConfigListConst relconfigs):
integrals(*coulomb_integrals), configs(relconfigs), M(nullptr), Sz(*coulomb_integrals)//, hamiltonian(coulomb_integrals)//, include_sigma3(false)
{
    // Set up matrix
    N = configs->NumCSFs();

    *logstream << " " << N << " " << std::flush;
    *outstream << " Number of CSFs = " << N << std::endl;
}

HamiltonianMatrix::~HamiltonianMatrix()
{
    if(M)
        delete M;
}

void HamiltonianMatrix::GenerateMatrix()
{
    if(M == nullptr)
    {   if(N <= SMALL_MATRIX_LIM)
            M = new SmallMatrix(N);
        else
            M = new SymMatrix(N);
    }
    else
        M->Clear();

    M->WriteMode(true);

    HamiltonianOperator H_operator(integrals);
    ManyBodyOperator<HamiltonianOperator*, HamiltonianOperator*> H(&H_operator, &H_operator);

    // Loop through projections
    auto proj_it = configs->projection_begin();
    while(proj_it != configs->projection_end())
    {
        auto proj_jt = proj_it;

        while(proj_jt != configs->projection_end())
        {
            double operatorH = H.GetMatrixElement(*proj_it, *proj_jt);

            if(fabs(operatorH) > 1.e-15)
            {
                for(auto coeff_i = proj_it.CSF_begin(); coeff_i != proj_it.CSF_end(); coeff_i++)
                {
                    RelativisticConfigList::const_CSF_iterator start_j = proj_jt.CSF_begin();

                    if(proj_it == proj_jt)
                        start_j = coeff_i;

                    for(auto coeff_j = start_j; coeff_j != proj_jt.CSF_end(); coeff_j++)
                    {
                        // See notes for an explanation
                        int i = coeff_i.index();
                        int j = coeff_j.index();

                        if(i < j)
                            M->At(i, j) += operatorH * (*coeff_i) * (*coeff_j);
                        else if(i > j)
                            M->At(j, i) += operatorH * (*coeff_i) * (*coeff_j);
                        else if(proj_it == proj_jt)
                            M->At(i, j) += operatorH * (*coeff_i) * (*coeff_j);
                        else
                            M->At(i, j) += 2. * operatorH * (*coeff_i) * (*coeff_j);
                    }
                }
            }
            proj_jt++;
        }
        proj_it++;
    }
}

//double HamiltonianMatrix::GetProjectionH(const Projection& first, const Projection& second) const
//{
//    // Create mutable vectors of pointers to projection elements
//    std::vector<const ElectronInfo*> left, right;
//    for(auto& it: first)
//        left.push_back(&it);
//    for(auto& it: second)
//        right.push_back(&it);
//
//    int numdiff = GetProjectionDifferences(left, right);
//
//    int sign;
//    if(numdiff >= 0)
//        sign = 1;
//    else
//        sign = -1;
//
//    double value = 0.;
//
//    if(numdiff == 0)
//    {
//        // Sum <i|f|i>
//        for(unsigned int i=0; i<left.size(); i++)
//        {
//            value += integrals.GetOneElectronIntegral(*left[i], *left[i]);
//            
//            // Sum(i < j) <ij|g|ij> - <ij|g|ji>
//            for(unsigned int j=i+1; j<first.size(); j++)
//            {
//                value += CoulombMatrixElement(*left[i], *left[j], *left[i], *left[j])
//                         - CoulombMatrixElement(*left[i], *left[j], *left[j], *left[i]);
//            }
//        }
//    }
//    else if(abs(numdiff) == 1)
//    {
//        const ElectronInfo& f1 = *left[0];
//        const ElectronInfo& s1 = *right[0];
//
//        // a->b
//        // <a|f|b>
//        if((f1.M() == s1.M()) && (f1.Kappa() == s1.Kappa()))
//            value = integrals.GetOneElectronIntegral(f1, s1) * sign;
//
//        // Sum(e) <ae|g|be> - <ae|g|eb>
//        for(unsigned int i = 1; i < left.size(); i++)
//        {
//            const ElectronInfo& e(*left[i]);
//            value += sign * (CoulombMatrixElement(f1, e, s1, e)
//                            - CoulombMatrixElement(f1, e, e, s1));
//        }
//    }
//    else if(abs(numdiff) == 2)
//    {
//        const ElectronInfo& f1 = *left[0];
//        const ElectronInfo& s1 = *right[0];
//        const ElectronInfo& f2 = *left[1];
//        const ElectronInfo& s2 = *right[1];
//
//        // a->b, c->d
//        // <ac|g|bd> - <ac|g|db>
//        value = sign * (CoulombMatrixElement(f1, f2, s1, s2)
//                        - CoulombMatrixElement(f1, f2, s2, s1));
//    }
//
////    if(include_sigma3)
////        value += GetSigma3(first, second);
//
//    return value;
//}

double HamiltonianOperator::GetMatrixElement(const ElectronInfo& e1, const ElectronInfo& e2, const ElectronInfo& e3, const ElectronInfo& e4) const
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

    MathConstant* constants = MathConstant::Instance();

    while(k <= kmax)
    {
        double coeff = 0.;
        if(fabs(q) <= k)
            coeff = constants->Electron3j(e1.TwoJ(), e3.TwoJ(), k, -e1.TwoM(), e3.TwoM()) *
                    constants->Electron3j(e2.TwoJ(), e4.TwoJ(), k, -e2.TwoM(), e4.TwoM());
            
        if(coeff)
            coeff = coeff * constants->Electron3j(e1.TwoJ(), e3.TwoJ(), k, 1, -1) *
                            constants->Electron3j(e2.TwoJ(), e4.TwoJ(), k, 1, -1);

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
                coeff = constants->Electron3j(e1.TwoJ(), e3.TwoJ(), k, -e1.TwoM(), e3.TwoM()) *
                        constants->Electron3j(e2.TwoJ(), e4.TwoJ(), k, -e2.TwoM(), e4.TwoM());
                
            if(coeff)
                coeff = coeff * constants->Electron3j(e1.TwoJ(), e3.TwoJ(), k, 1, -1) *
                                constants->Electron3j(e2.TwoJ(), e4.TwoJ(), k, 1, -1);

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

/*
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
        for(unsigned int i=0; i<first.size(); i++)
        {
            // Sum(i < j) <ij|g|ij> - <ij|g|ji>
            for(unsigned int j=i+1; j<first.size(); j++)
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
        for(unsigned int i=0; i<first.size(); i++)
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

    MathConstant* constants = MathConstant::Instance();

    // k == 1 only
    double coeff = 0.;
    if(fabs(q) <= k)
        coeff = constants->Electron3j(e1.TwoJ(), e3.TwoJ(), k, -e1.TwoM(), e3.TwoM()) *
                constants->Electron3j(e2.TwoJ(), e4.TwoJ(), k, -e2.TwoM(), e4.TwoM());
        
    if(coeff)
        coeff = coeff * constants->Electron3j(e1.TwoJ(), e3.TwoJ(), k, 1, -1) *
                        constants->Electron3j(e2.TwoJ(), e4.TwoJ(), k, 1, -1);

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
*/

void HamiltonianMatrix::WriteToFile(const std::string& filename)
{
    return;
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

void HamiltonianMatrix::SolveMatrix(const Symmetry& sym, unsigned int num_solutions, pLevelMap levels, bool get_gFactors)
{
    if(N == 0)
    {   *outstream << "\nNo solutions" << std::endl;
        return;
    }

    M->WriteMode(false);

    *outstream << "\nFinding solutions" << std::endl;

    unsigned int NumSolutions = mmin(num_solutions, N);
    
    double* V = new double[NumSolutions * N];
    double* E = new double[NumSolutions];

    Eigensolver solver;
    solver.SolveLargeSymmetric(M, E, V, N, NumSolutions);

    for(unsigned int i = 0; i < NumSolutions; i++)
    {
        pLevel level(new Level(E[i], (V + N * i), configs, N));
        (*levels)[LevelID(sym, i)] = level;
    }

    delete[] E;
    delete[] V;

    // Calculate g-Factors
    if(get_gFactors && sym.GetTwoJ())
        GetgFactors(sym, levels);
}

/*
void HamiltonianMatrix::GetEigenvalues(const Eigenstates& eigenstates) const
{
    unsigned int NumSolutions = eigenstates.GetNumEigenvectors();
    const double* V = eigenstates.GetEigenvectors();

    double* total = new double[NumSolutions];
    double* coeff = new double[NumSolutions];
    unsigned int solution;
    for(solution = 0; solution < NumSolutions; solution++)
        total[solution] = 0.;

    // Iterate over different relativistic configurations
    unsigned int i=0, j;
    RelativisticConfigList::const_iterator list_it = configs->begin();
    while(list_it != configs->end())
    {
        const ProjectionSet& proj_i = list_it->GetProjections();
        unsigned int proj_i_size = proj_i.size();
        unsigned int num_states_i = list_it->NumJStates();
        const double* coefficients_i = list_it->GetJCoefficients();

        RelativisticConfigList::const_iterator list_jt = list_it;
        j = i;
        while(list_jt != configs->end())
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
            << total[solution]*MathConstant::Instance()->HartreeEnergyInInvCm()
            << " /cm" << std::endl;

    delete[] total;
    delete[] coeff;
}
*/

void HamiltonianMatrix::GetgFactors(const Symmetry& sym, pLevelMap levels) const
{
    unsigned int NumSolutions = levels->size(sym);
    if(NumSolutions == 0)
        return;

    const SzOperator* pSz = &Sz;
    ManyBodyOperator<const SzOperator*> many_body_Sz(pSz);

    std::vector<double> total_Sz = many_body_Sz.GetMatrixElement(levels->begin(sym), levels->end(sym));

    unsigned int solution = 0;
    auto it = levels->begin(sym);
    while(solution < NumSolutions && it != levels->end(sym))
    {
        it->second->SetgFactor(total_Sz[solution]/sym.GetJ() + 1.);
        solution++;
        it++;
    }
}

//double HamiltonianMatrix::GetSz(const ElectronInfo& e) const
//{
//    if(e.Kappa() > 0)   // J = L - 1/2
//        return -e.M()/(2.*e.L() + 1.);
//    else
//        return e.M()/(2.*e.L() + 1.);
//}

double HamiltonianMatrix::SzOperator::GetMatrixElement(const ElectronInfo& e1, const ElectronInfo& e2) const
{
    double val = 0.;

    if((e1.L() == e2.L()) &&
       (e1.TwoM() == e2.TwoM()))
    {
        double overlap = integrals.GetOverlapIntegral(e1, e2);

        if(e1.Kappa() == e2.Kappa())
        {   val = e1.M()/(2*e1.L() + 1.);

            if(e1.Kappa() > 0)
                val = -val;
        }
        else
        {   double Lplushalf = double(e1.L()) + 0.5;
            double M = e1.M();
            val = -sqrt((Lplushalf + M)*(Lplushalf - M))/(2.*Lplushalf);
        }
        val = val * overlap;
    }

    return val;
}

//double HamiltonianMatrix::GetSigma3(const Projection& first, const Projection& second) const
//{
//    const std::set<Configuration>* leading_configs = confgen->GetLeadingConfigs();
//
//#ifdef SIGMA3_AND
//    // Check that first AND second are leading configurations
//    if((leading_configs->find(first.GetNonRelConfiguration()) == leading_configs->end()) ||
//       (leading_configs->find(second.GetNonRelConfiguration()) == leading_configs->end()))
//        return 0.;
//#else
//    // Check that first OR second is a leading configuration
//    if((leading_configs->find(first.GetNonRelConfiguration()) == leading_configs->end()) &&
//       (leading_configs->find(second.GetNonRelConfiguration()) == leading_configs->end()))
//        return 0.;
//#endif
//
//    unsigned int diff[6];
//    int numdiff = Projection::GetProjectionDifferences3(first, second, diff);
//
//    int sign;
//    if(numdiff >= 0)
//        sign = 1;
//    else
//        sign = -1;
//
//    double value = 0.;
//
//    if(numdiff == 0)
//    {
//        // Sum(i < j < k) Sigma3(ijk, ijk)
//        for(unsigned int i=0; i<first.size(); i++)
//        {
//            for(unsigned int j=i+1; j<first.size(); j++)
//            {
//                for(unsigned int k=j+1; k<first.size(); k++)
//                {
//                    value += Sigma3(first[i], first[j], first[k], first[i], first[j], first[k]);
//                }
//            }
//        }
//    }
//    else if(abs(numdiff) == 1)
//    {
//        const ElectronInfo& f1 = first[diff[0]];
//        const ElectronInfo& s1 = second[diff[1]];
//
//        // Sum(i < j) Sigma3(aij, bij)
//        for(unsigned int i=0; i<first.size(); i++)
//        {
//            for(unsigned int j=i+1; j<first.size(); j++)
//            {
//                if((i != diff[0]) && (j != diff[0]))
//                {   
//                    value += sign * Sigma3(f1, first[i], first[j], s1, first[i], first[j]);
//                }
//            }
//        }
//    }
//    else if(abs(numdiff) == 2)
//    {
//        const ElectronInfo& f1 = first[diff[0]];
//        const ElectronInfo& s1 = second[diff[1]];
//        const ElectronInfo& f2 = first[diff[2]];
//        const ElectronInfo& s2 = second[diff[3]];
//
//        // Sum(i) Sigma3(abi, cdi)
//        for(unsigned int i=0; i<first.size(); i++)
//        {
//            if((i != diff[0]) && (i != diff[2]))
//               value += sign * Sigma3(f1, f2, first[i], s1, s2, first[i]);
//        }
//    }
//    else if(abs(numdiff) == 3)
//    {
//        const ElectronInfo& f1 = first[diff[0]];
//        const ElectronInfo& s1 = second[diff[1]];
//        const ElectronInfo& f2 = first[diff[2]];
//        const ElectronInfo& s2 = second[diff[3]];
//        const ElectronInfo& f3 = first[diff[4]];
//        const ElectronInfo& s3 = second[diff[5]];
//
//        // Sigma3(abc, def)
//        value = sign * Sigma3(f1, f2, f3, s1, s2, s3);
//    }
//
//    return value;
//}
//
//double HamiltonianMatrix::Sigma3(const ElectronInfo& e1, const ElectronInfo& e2, const ElectronInfo& e3,
//           const ElectronInfo& e4, const ElectronInfo& e5, const ElectronInfo& e6) const
//{
//    // Check momentum projections
//    if(e1.TwoM() + e2.TwoM() + e3.TwoM() != e4.TwoM() + e5.TwoM() + e6.TwoM())
//        return 0.;
//
//    // Check parity
//    if((e1.L() + e2.L() + e3.L() + e4.L() + e5.L() + e6.L())%2)
//        return 0.;
//
//    double value = 0.;
//    
//    // The sign changes for odd permutations
//    value =   Sigma3LinePermutations(e1, e2, e3, e4, e5, e6)
//            + Sigma3LinePermutations(e1, e2, e3, e5, e6, e4)
//            + Sigma3LinePermutations(e1, e2, e3, e6, e4, e5)
//            - Sigma3LinePermutations(e1, e2, e3, e5, e4, e6)
//            - Sigma3LinePermutations(e1, e2, e3, e6, e5, e4)
//            - Sigma3LinePermutations(e1, e2, e3, e4, e6, e5);
//
//    return value;
//}
//
///** This function does the line permutations, putting the pairs on different levels
//    of the three-body interaction.
//    */
//inline double HamiltonianMatrix::Sigma3LinePermutations(const ElectronInfo& e1, const ElectronInfo& e2, const ElectronInfo& e3,
//           const ElectronInfo& e4, const ElectronInfo& e5, const ElectronInfo& e6) const
//{
//    double value = 0.;
//
//    // There are no sign changes, since there are the same number
//    // of permutations on both sides
//
//    value =   sigma3calc->GetSecondOrderSigma3(e1, e2, e3, e4, e5, e6)
//            + sigma3calc->GetSecondOrderSigma3(e2, e3, e1, e5, e6, e4)
//            + sigma3calc->GetSecondOrderSigma3(e3, e1, e2, e6, e4, e5)
//            + sigma3calc->GetSecondOrderSigma3(e3, e2, e1, e6, e5, e4)
//            + sigma3calc->GetSecondOrderSigma3(e2, e1, e3, e5, e4, e6)
//            + sigma3calc->GetSecondOrderSigma3(e1, e3, e2, e4, e6, e5);
//
//    return value;
//}

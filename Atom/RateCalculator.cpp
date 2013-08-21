#include "RateCalculator.h"
#include "Include.h"
#include "Universal/MathConstant.h"
#include "HartreeFock/StateIntegrator.h"
#include "Universal/CoulombIntegrator.h"
#include "Universal/Integrator.h"
#include "Universal/Enums.h"
#include "HartreeFock/ContinuumBuilder.h"

#include <complex>
#include <boost/math/special_functions/bessel.hpp>
#include <boost/math/special_functions/factorials.hpp>

RateCalculator::RateCalculator(ExcitedStates* basis):
    excited(basis)
{
    NumStates = excited->NumStates();

    ConstStateIterator it_i = excited->GetConstStateIterator();
    unsigned int i;

    // Iterate through states, assign in order
    it_i.First(); i = 0;
    while(!it_i.AtEnd())
    {
        state_index.insert(std::pair<OrbitalInfo, unsigned int>(OrbitalInfo(it_i.GetState()), i));
        it_i.Next(); i++;
    }
}

void RateCalculator::DielectronicRecombination(Atom* A)
{
    FILE* fp;
    if(ProcessorRank == 0)
        fp = fopen("dr.txt", "wt");

    double Ec_width = 0.01; // Beam energy width (eV)
    Ec_width = Ec_width/MathConstant::Instance()->HartreeEnergyIneV();

    Eigenstates* gs_eigenstates = A->GetEigenstates(Symmetry(0, even));

    // C4+:
    //double continuum_energy = -3952437.7/MathConstant::Instance()->HartreeEnergyInInvCm();
    //      gs_eigenstates->GetEigenvalues()[0] + 198310.6672/MathConstant::Instance()->HartreeEnergyInInvCm();

    // C2+, HF(bare), E(2s):
    //double continuum_energy = -2.36589786067;
    // C2+, HF(2s), E(GS + 386241.0):
    double continuum_energy = gs_eigenstates->GetEigenvalues()[0] + 386241.0/MathConstant::Instance()->HartreeEnergyInInvCm();

    double energy_limit = -2.2;//continuum_energy + 100.;//2.0/MathConstant::Instance()->HartreeEnergyIneV();
    *outstream << continuum_energy * MathConstant::Instance()->HartreeEnergyInInvCm() << " -> "
               << energy_limit * MathConstant::Instance()->HartreeEnergyInInvCm() << std::endl;

    double auger_rate, radiative_rate;

    for(int parity_loop = 0; parity_loop < 2; parity_loop++)
    {   Parity p = even;
        if(parity_loop)
            p = odd;

        for(unsigned int two_j = 0; two_j <= 10; two_j+=2)
        {   Symmetry sym(two_j, p);
            Eigenstates* ES = A->GetEigenstates(sym);
            if(ES)
            {
                unsigned int i = 0;
                unsigned int num_eigenvectors = ES->GetNumEigenvectors();
                while((i < num_eigenvectors) && (ES->GetEigenvalues()[i] < continuum_energy))
                    i++;

                while((i < num_eigenvectors) && (ES->GetEigenvalues()[i] < energy_limit))
                {
                    *outstream << "\n*** " << sym.GetString() << " ***" << std::endl;
                    ES->Print(i);
                    auger_rate = CalculateAugerRate(A, sym, i, continuum_energy);
                    *outstream << "  E = " << ES->GetEigenvalues()[i] * MathConstant::Instance()->HartreeEnergyInInvCm()
                               << "  W = " << 6.5855e-16*auger_rate/MathConstant::Instance()->HartreeEnergyIneV() << std::endl;

                    radiative_rate = CalculateAllDipoleStrengths(A, sym, i);
                    
                    double Ec = ES->GetEigenvalues()[i] - continuum_energy;
                    double cross_section = 6.68525e-15/(Ec*Ec_width) * (two_j + 1.)/4.
                                                                // Denominator is 2 * (stat. weight of N-electron target)
                                                                //      = 2 * [ 2(1/2) + 1] = 4
                                                                // for j = 1/2 ground state
                                           * auger_rate * radiative_rate
                                           /(auger_rate + radiative_rate);

                    if(ProcessorRank == 0)
                        fprintf(fp, "%12.6e\t%12.6e\t%s:%d\n", Ec * MathConstant::Instance()->HartreeEnergyIneV(), cross_section,
                            sym.GetString().c_str(), i);
                    i++;
                }
            }
        }
    }
    
    if(ProcessorRank == 0)
        fclose(fp);
}

void RateCalculator::AutoionisationRates(Atom* A)
{
    Eigenstates* gs_eigenstates = A->GetEigenstates(Symmetry(0, even));

    // C4+:
    //double continuum_energy = -3952437.7/MathConstant::Instance()->HartreeEnergyInInvCm();
    //      gs_eigenstates->GetEigenvalues()[0] + 198310.6672/MathConstant::Instance()->HartreeEnergyInInvCm();

    // C2+, HF(bare), E(2s):
    //double continuum_energy = -2.36589786067;
    // C2+, HF(2s), E(GS + 386241.0):
    double continuum_energy = gs_eigenstates->GetEigenvalues()[0] + 386241.0/MathConstant::Instance()->HartreeEnergyInInvCm();

    double energy_limit = continuum_energy + 10000.0/MathConstant::Instance()->HartreeEnergyInInvCm();
    *outstream << continuum_energy*MathConstant::Instance()->HartreeEnergyInInvCm() << " -> "
               << energy_limit*MathConstant::Instance()->HartreeEnergyInInvCm() << std::endl;

    double auger_rate;

    for(int parity_loop = 0; parity_loop < 2; parity_loop++)
    {   Parity p = even;
        if(parity_loop)
            p = odd;

        for(unsigned int two_j = 0; two_j <= 10; two_j+=2)
        {
            Symmetry sym(two_j, p);
            Eigenstates* ES = A->GetEigenstates(sym);
            if(ES)
            {
                unsigned int i = 0;
                unsigned int num_eigenvectors = ES->GetNumEigenvectors();
                while((i < num_eigenvectors) && (ES->GetEigenvalues()[i] < continuum_energy))
                    i++;

                while((i < num_eigenvectors) && (ES->GetEigenvalues()[i] < energy_limit))
                {
                    *outstream << "\nJ = " << double(two_j)/2. << ", P = ";
                    if(p == even)
                        *outstream << "even:" << std::endl;
                    else
                        *outstream << "odd:" << std::endl;
                    ES->Print(i);
                    auger_rate = CalculateAugerRate(A, sym, i, continuum_energy);
                    i++;
                }
            }
        }
    }
}

double RateCalculator::CalculateAllDipoleStrengths(Atom* A, Symmetry sym1, unsigned int solution1, bool print_rates, bool print_oscillator_strengths)
{
    Eigenstates* eigenstates1 = A->GetEigenstates(sym1);
    if(!eigenstates1 || !eigenstates1->Restore())
    {   *outstream << sym1.GetString() << " eigenstates restore failed" << std::endl;
        return 0.;
    }

    const double* V1 = eigenstates1->GetEigenvectors();
    const double* E1 = eigenstates1->GetEigenvalues();
    unsigned int N1 = eigenstates1->GetEigenvectorLength();
    RelativisticConfigList* configs1 = eigenstates1->GetRelConfigs();

    unsigned int diff[4];   // Storage for projection differences.

    unsigned int TwoJ_min, TwoJ_max;
    Parity other_P;

    if(sym1.GetTwoJ() == 0)
        TwoJ_min = 2;
    else if(sym1.GetTwoJ() == 1)
        TwoJ_min = 1;
    else
        TwoJ_min = sym1.GetTwoJ() - 2;

    TwoJ_max = sym1.GetTwoJ() + 2;

    if(sym1.GetParity() == even)
        other_P = odd;
    else
        other_P = even;

    double sum_over_all_J = 0.0;
    *outstream << std::setprecision(5);

    for(unsigned int TwoJ = TwoJ_min; TwoJ <= TwoJ_max; TwoJ += 2)
    {
        Symmetry sym2(TwoJ, other_P);
        Eigenstates* eigenstates2 = A->GetEigenstates(sym2);
        if(!eigenstates2 || !eigenstates2->Restore())
        {   *outstream << "WARNING: " << sym2.GetString() << "eigenstates restore failed" << std::endl;
            continue;
        }

        const double* V2 = eigenstates2->GetEigenvectors();
        const double* E2 = eigenstates2->GetEigenvalues();
        unsigned int N2 = eigenstates2->GetEigenvectorLength();
        RelativisticConfigList* configs2 = eigenstates2->GetRelConfigs();

        // Determine number of lower states
        unsigned int num_solutions2 = 0;
        while((num_solutions2 < eigenstates2->GetNumEigenvectors()) &&
              (E2[num_solutions2] < E1[solution1]))
        {
            num_solutions2++;
        }

        double* total = new double[num_solutions2];
        double* coeff = new double[num_solutions2];
        unsigned int solution2;
        memset(total, 0, num_solutions2*sizeof(double));

        // Iterate over different relativistic configurations
        unsigned int i=0, j;
        RelativisticConfigList::const_iterator list_it = configs1->begin();
        while(list_it != configs1->end())
        {
            const ProjectionSet& proj_i = list_it->GetProjections();
            unsigned int proj_i_size = proj_i.size();
            unsigned int num_states_i = list_it->NumJStates();
            const double* coefficients_i = list_it->GetJCoefficients();

            RelativisticConfigList::const_iterator list_jt = configs2->begin();
            j = 0;
            while(list_jt != configs2->end())
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
                        // <pi| r | pj>
                        double matrix_element = 0.;
                        int num_diff = Projection::GetProjectionDifferences(*pi_it, *pj_it, diff);

                        if(abs(num_diff) == 1)
                        {
                            matrix_element = GetE1MatrixElement((*pi_it)[diff[0]], (*pj_it)[diff[1]]);
                            if(num_diff == -1)
                                matrix_element = -matrix_element;
                        }

                        // coefficients
                        if(matrix_element)
                        {
                            memset(coeff, 0, num_solutions2*sizeof(double));

                            for(unsigned int jstate_i = 0; jstate_i < num_states_i; jstate_i++)
                            {
                                for(unsigned int jstate_j = 0; jstate_j < num_states_j; jstate_j++)
                                {
                                    for(solution2 = 0; solution2 < num_solutions2; solution2++)
                                    {
                                        coeff[solution2] += coefficients_i[jstate_i*proj_i_size + pi]
                                            * coefficients_j[jstate_j*proj_j_size + pj]
                                            * V1[solution1*N1 + i + jstate_i]
                                            * V2[solution2*N2 + j + jstate_j];
                                    }
                                }
                            }
                            
                            for(solution2 = 0; solution2 < num_solutions2; solution2++)
                                total[solution2] += coeff[solution2] * matrix_element;
                        }

                        pj_it++; pj++;
                    }

                    pi_it++; pi++;
                }

                list_jt++; j+=num_states_j;
            }

            list_it++; i+=num_states_i;
        }

        double wigner_coeff = MathConstant::Instance()->Electron3j(sym1.GetTwoJ(), TwoJ, 1, -sym1.GetTwoJ(), TwoJ);
        for(solution2 = 0; solution2 < num_solutions2; solution2++)
        {   // Get reduced matrix element squared
            total[solution2] = total[solution2]/wigner_coeff;
            total[solution2] = total[solution2] * total[solution2];
        }

        if(print_oscillator_strengths)
        {   *outstream << "(" << sym1.GetString() << ": " << solution1 << ") -> (" << sym2.GetString()
                       << ": i) reduced E1 amplitudes and weighted oscillator strengths (g.f)" << std::endl;

            for(solution2 = 0; solution2 < num_solutions2; solution2++)
            {   // Oscillator strength
                double deltaE = fabs(E1[solution1] - E2[solution2]);
                *outstream << "  " << solution2 << "  " << sqrt(total[solution2])
                           << "  " << deltaE * total[solution2]/3. << std::endl;
            }
        }

        if(print_rates)
            *outstream << "(" << sym1.GetString() << ": " << solution1 << ") -> (" << sym2.GetString()
                       << ": i) transition probabilities (A) (/sec)" << std::endl;

        double sum = 0.0;
        for(solution2 = 0; solution2 < num_solutions2; solution2++)
        {   // Radiative rate (per second)
            double sig = fabs(E1[solution1] - E2[solution2]) * MathConstant::Instance()->HartreeEnergyInInvCm();
            double rate = total[solution2] * sig * sig * sig * 2.0261e-6/(sym1.GetTwoJ() + 1);
            if(print_rates)
                *outstream << "  " << solution2 << "  " << rate << std::endl;
            sum += rate;
        }

        if(print_rates)
            *outstream << "  Sum = " << sum << std::endl;
        sum_over_all_J += sum;

        delete[] coeff;
        delete[] total;
    }
    
    if(print_rates)
        *outstream << "Radiative rate (sum over all J) = " << sum_over_all_J << std::endl;

    return sum_over_all_J;
}

double RateCalculator::CalculateDipoleStrength(Atom* A, Symmetry sym1, unsigned int solution1, Symmetry sym2, unsigned int solution2)
{
    Eigenstates* eigenstates1 = A->GetEigenstates(sym1);
    if(!eigenstates1 || !eigenstates1->Restore())
    {   *errstream << "eigenstates1 restore failed " << std::endl;
        return 0.0;
    }
    
    Eigenstates* eigenstates2 = A->GetEigenstates(sym2);
    if(!eigenstates2 || !eigenstates2->Restore())
    {   *errstream << "eigenstates2 restore failed" << std::endl;
        return 0.0;
    }
    
    const double* V1 = eigenstates1->GetEigenvectors();
    unsigned int N1 = eigenstates1->GetEigenvectorLength();

    const double* V2 = eigenstates2->GetEigenvectors();
    unsigned int N2 = eigenstates2->GetEigenvectorLength();

    RelativisticConfigList* configs1 = eigenstates1->GetRelConfigs();
    RelativisticConfigList* configs2 = eigenstates2->GetRelConfigs();

    unsigned int diff[4];   // Storage for projection differences.
    double total = 0.0;

    // Iterate over different relativistic configurations
    unsigned int i=0, j;
    RelativisticConfigList::const_iterator list_it = configs1->begin();
    while(list_it != configs1->end())
    {
        const ProjectionSet& proj_i = list_it->GetProjections();
        unsigned int proj_i_size = proj_i.size();
        unsigned int num_states_i = list_it->NumJStates();
        const double* coefficients_i = list_it->GetJCoefficients();

        RelativisticConfigList::const_iterator list_jt = configs2->begin();
        j = 0;
        while(list_jt != configs2->end())
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
                    // <pi| r | pj>
                    double matrix_element = 0.;
                    int num_diff = Projection::GetProjectionDifferences(*pi_it, *pj_it, diff);

                    if(abs(num_diff) == 1)
                    {
                        matrix_element = GetE1MatrixElement((*pi_it)[diff[0]], (*pj_it)[diff[1]]);
                        if(num_diff == -1)
                            matrix_element = -matrix_element;
                    }

                    // coefficients
                    if(matrix_element)
                    {
                        double coeff = 0.;
                        for(unsigned int jstate_i = 0; jstate_i < num_states_i; jstate_i++)
                        {
                            for(unsigned int jstate_j = 0; jstate_j < num_states_j; jstate_j++)
                            {
                                coeff += coefficients_i[jstate_i*proj_i_size + pi]
                                        * coefficients_j[jstate_j*proj_j_size + pj]
                                        * V1[solution1*N1 + i + jstate_i]
                                        * V2[solution2*N2 + j + jstate_j];
                                //*outstream << "Coefficient = " << coefficients_i[jstate_i*proj_i_size + pi] * coefficients_j[jstate_j*proj_j_size + pj] * V1[solution1*N1 + i + jstate_i] * V2[solution2*N2 + j + jstate_j] << " = " << coefficients_i[jstate_i*proj_i_size + pi] << " * " << coefficients_j[jstate_j*proj_j_size + pj] << " * " << V1[solution1*N1 + i + jstate_i] << " * " << V2[solution2*N2 + j + jstate_j] << std::endl;
                            }
                        }
                        
                        //*outstream << list_it->Name() << " to " << list_jt->Name() << ": " << matrix_element * coeff << " (" << matrix_element << " * " << coeff << ")" << std::endl;
                        total += matrix_element * coeff;
                    }

                    pj_it++; pj++;
                }

                pi_it++; pi++;
            }

            list_jt++; j+=num_states_j;
        }

        list_it++; i+=num_states_i;
    }

    *outstream << std::setprecision(5);

    //*outstream << "jcb simple sum = " << total << std::endl;
    //*outstream << "3j modifier = " << 1.0/MathConstant::Instance()->Electron3j(sym1.GetTwoJ(), sym2.GetTwoJ(), 1, -sym1.GetTwoJ(), sym2.GetTwoJ()) << std::endl;
    total = total/MathConstant::Instance()->Electron3j(sym1.GetTwoJ(), sym2.GetTwoJ(), 1, -sym1.GetTwoJ(), sym2.GetTwoJ());
    *outstream << "E1 reduced matrix element = " << total << std::endl;
    total = total * total;
    *outstream << "E1 reduced matrix element squared (S) = " << total << std::endl;
    
    double deltaE = fabs((eigenstates1->GetEigenvalues())[solution1] - (eigenstates2->GetEigenvalues())[solution2]);
    *outstream << "Weighted oscillator strength (g.f) = " << 2. * deltaE * total/3. << std::endl;

    double sig = deltaE * MathConstant::Instance()->HartreeEnergyInInvCm();
    *outstream << "sigma = " << sig <<std::endl;
    *outstream << "Weighted transition prob. (g.A) (/sec) = " << total * sig * sig * sig * 2.0261e-6 << std::endl;

    return total;
}

double RateCalculator::CalculateMultipoleStrength(MultipolarityType::Enum MType, unsigned int J, Atom* A, Symmetry sym1, unsigned int solution1, Symmetry sym2, unsigned int solution2, TransitionGaugeType::Enum aGaugeType, TransitionCalculationMethod::Enum aMethod)
{
    if(MType != MultipolarityType::M && MType != MultipolarityType::E)
    {
        *errstream << "RateCalculator::multipole_moment_integral: Unknown MultipolarityType" << (int) MType << "received" << std::endl;
        return 0.0;
    }
    
    Eigenstates* eigenstates1 = A->GetEigenstates(sym1);
    if(!eigenstates1 || !eigenstates1->Restore())
    {   *errstream << "eigenstates1 restore failed " << std::endl;
    return 0.0;
    }
    
    Eigenstates* eigenstates2 = A->GetEigenstates(sym2);
    if(!eigenstates2 || !eigenstates2->Restore())
    {   *errstream << "eigenstates2 restore failed" << std::endl;
    return 0.0;
    }
    
    const double* V1 = eigenstates1->GetEigenvectors();
    unsigned int N1 = eigenstates1->GetEigenvectorLength();
    
    const double* V2 = eigenstates2->GetEigenvectors();
    unsigned int N2 = eigenstates2->GetEigenvectorLength();
    
    RelativisticConfigList* configs1 = eigenstates1->GetRelConfigs();
    RelativisticConfigList* configs2 = eigenstates2->GetRelConfigs();
    
    unsigned int diff[4];   // Storage for projection differences.
    std::complex<double> total = 0.0;
    
    // Iterate over different relativistic configurations
    unsigned int i=0, j;
    RelativisticConfigList::const_iterator list_it = configs1->begin();
    while(list_it != configs1->end())
    {
        const ProjectionSet& proj_i = list_it->GetProjections();
        unsigned int proj_i_size = proj_i.size();
        unsigned int num_states_i = list_it->NumJStates();
        const double* coefficients_i = list_it->GetJCoefficients();
        
        RelativisticConfigList::const_iterator list_jt = configs2->begin();
        j = 0;
        while(list_jt != configs2->end())
        {
            //*outstream << "Relativisitic config from" << list_it->Name() << " to" << list_jt->Name() << std::endl;
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
                    // <pi| q | pj>
                    std::complex<double> matrix_element = 0.0;
                    int num_diff = Projection::GetProjectionDifferences(*pi_it, *pj_it, diff);
                    //*outstream << pi_it->Name() << " to " << pj_it->Name() << std::endl;

                    // If the configuration including projections are exactly the same
                    if(num_diff == 0)
                    {
                        // Determine if the operator acts on states with
                        // the same parity, if it does,
                        // apply the operator to all the electrons
                        // Warning: assumes all core states are closed shells!

                        if(MType == MultipolarityType::M)
                        {
                            for(unsigned int i=0; i<configs1->front().NumParticles(); i++)
                            {
                                matrix_element += GetMJMatrixElement(A, J, (*pi_it)[i], (*pi_it)[i]);
                            }
                        }
                        else if(MType == MultipolarityType::E)
                        {
                            for(unsigned int i=0; i<configs1->front().NumParticles(); i++)
                            {
                                matrix_element += GetEJMatrixElement(J, (*pi_it)[i], (*pi_it)[i]);
                            }
                        }
                        //*outstream << "Matrix element for 0 diff: " << matrix_element << std::endl;
                    }
                    
                    // If there is exactly one electron different between the configurations
                    // Act on the different electron with the operator
                    if(abs(num_diff) == 1)
                    {
                        if(MType == MultipolarityType::M)
                        {
                            matrix_element += GetMJMatrixElement(A, J, (*pi_it)[diff[0]], (*pj_it)[diff[1]]);
                        }
                        else if(MType == MultipolarityType::E)
                        {
                            if(aMethod == TransitionCalculationMethod::WalterJohnson)
                            {
                                matrix_element = GetEJMatrixElement(J, (*pi_it)[diff[0]], (*pj_it)[diff[1]], aGaugeType);
                            }
                            else if(aMethod == TransitionCalculationMethod::ZenonasRudzikas)
                            {
                                matrix_element += GetEJMatrixElementZenonas(J, (*pi_it)[diff[0]], (*pj_it)[diff[1]], aGaugeType);
                            }
                        }
                        
                        // Reverse the sign of the matrix element if the state I (X + e1) required
                        // an odd number of states to reach state J (X + e2)
                        // where X have exactly the same order and e1 and e2 are the different electrons
                        if(num_diff == -1)
                        {
                            matrix_element = -matrix_element;
                        }
                    }
                    
                    // The correct coefficient for the state was the amount of mixing calculated
                    // from previous calculations
                    if(std::real(matrix_element) != 0.0 || std::imag(matrix_element) != 0.0)
                    {
                        double coeff = 0.;
                        for(unsigned int jstate_i = 0; jstate_i < num_states_i; jstate_i++)
                        {
                            for(unsigned int jstate_j = 0; jstate_j < num_states_j; jstate_j++)
                            {                                                                                                                                   
                                coeff += coefficients_i[jstate_i*proj_i_size + pi] * coefficients_j[jstate_j*proj_j_size + pj] * V1[solution1*N1 + i + jstate_i] * V2[solution2*N2 + j + jstate_j];
                            }
                        }
                        
                        //if(fabs(matrix_element * coeff) > 0.05)
                        //{
                            //*outstream << list_it->Name() << " to " << list_jt->Name() << ": " << matrix_element * coeff << " (" << matrix_element << " * " << coeff << ")" << std::endl;
                            //*outstream << "Contribution = " << matrix_element * coeff  << " = " << matrix_element << " * " << coeff << std::endl;
                        //}
                        total += (matrix_element * coeff);
                    }
                    
                    pj_it++; 
                    pj++;
                }
                
                pi_it++;
                pi++;
            }
            
            list_jt++; j+=num_states_j;
        }
        
        list_it++; i+=num_states_i;
    }
    
    *outstream << std::setprecision(5);
    // Use Wigner-Eckart to convert from matrix element to reduced matrix element
    // Note that the projection of the state is always M = J
    total = total/MathConstant::Instance()->Electron3j(sym2.GetTwoJ(), sym1.GetTwoJ(), J, sym2.GetTwoJ(), -sym1.GetTwoJ())/pow(-1.0, (((double) sym1.GetTwoJ())/2.0) - (((double) sym1.GetTwoJ())/2.0));

    if(std::imag(total) == 0)
    {
        *outstream << MultipolarityType::Name(MType) << J << " reduced matrix element = " << double(std::real(total)) << std::endl;
    }
    else
    {
        *outstream << MultipolarityType::Name(MType) << J << " reduced matrix element = " << total << std::endl;
    }
    total = std::norm(total);
    *outstream << MultipolarityType::Name(MType) << J << " reduced matrix element squared (S) = " << std::real(total) << std::endl;
    
    double deltaE = fabs((eigenstates1->GetEigenvalues())[solution1] - (eigenstates2->GetEigenvalues())[solution2]);
    //*outstream << "Weighted oscillator strength (g.f) = " << 2. * deltaE * total/3. << std::endl;
    
    double sig = deltaE * MathConstant::Instance()->HartreeEnergyInInvCm();
    *outstream << "sigma = " << sig << std::endl;
    
    double initialJ = sym1.GetJ();
    double finalJ = sym2.GetJ();
    
    if((eigenstates1->GetEigenvalues())[solution1] < (eigenstates2->GetEigenvalues())[solution2] )
    {
        initialJ = sym2.GetJ();
        finalJ = sym1.GetJ();
    }
    
    // Einstein A coefficient
    double EinsteinA = ((2.0 * ((double) J)) + 2.0) * ((2.0 * ((double) J)) + 1.0) /((double) J);
    EinsteinA /= boost::math::double_factorial<double>((2.0 * ((double) J)) + 1.0);
    EinsteinA /= boost::math::double_factorial<double>((2.0 * ((double) J)) + 1.0);
    
    // Convert from wavenumber representation to wavelength representation
    // Convert wavelength from atomic units to meters (factor of Bohr radius)
    // Convert 1/m to 1/cm
    EinsteinA *= pow(2.0 * MathConstant::Instance()->Pi() * MathConstant::Instance()->BohrRadiusSI(), (2.0 * ((double) J)) + 1.0);
    EinsteinA *= pow(100.0, (2.0 * ((double) J)) + 1.0);
    
    EinsteinA /= ((2.0 * initialJ) + 1.0);
    
    // The magnetic operator is defined with additional factor of 2c
    // This must be removed in order to calculate the Einstein A coefficient
    if(MType == MultipolarityType::M)
    {
        EinsteinA *= PhysicalConstant::Instance()->GetAlphaSquared();
        EinsteinA /= 4.0;
    }

    // Sigma = 1/wavelength in cm-1
    EinsteinA *= pow(sig, ((2.0 * ((double) J))) + 1.0);
    EinsteinA *= std::real(total);
    
    // Convert from A.U. to SI
    EinsteinA *= MathConstant::Instance()->AtomicFrequencySI();
    
    // Equal to the Einstein A coefficient
    *outstream << "Probability per unit time for emission (/sec) = " << EinsteinA << std::endl;
    
    // Degeneracy of final state times Einstein A coefficient
    *outstream << "Weighted probability per unit time (/sec) = " << EinsteinA * ((2.0 * initialJ) + 1.0) << std::endl;
    
    return std::real(total);
}

double RateCalculator::GetE1MatrixElement(const ElectronInfo& e1, const ElectronInfo& e2)
{
    *outstream << "Electron from " << e1.Name() << " to " << e2.Name() << std::endl;
    double matrix_element = 0.0;

    // Check e1.L() + e2.L() + 1 is even
    if((e1.L() + e2.L())%2)
    {
        double coeff = MathConstant::Instance()->Electron3j(e1.TwoJ(), e2.TwoJ(), 1, -e1.TwoM(), e2.TwoM());
        
        if(coeff)
        {   coeff = coeff * MathConstant::Instance()->Electron3j(e1.TwoJ(), e2.TwoJ(), 1, 1, -1)
                          * sqrt(double(e1.MaxNumElectrons() * e2.MaxNumElectrons()));


            if((abs(e1.TwoM() + 1)/2)%2)
                coeff = -coeff;

            unsigned int key = state_index[e1] * NumStates + state_index[e2];
            double overlap = 0.;

            // Check that this integral doesn't already exist
            if(E1Integrals.find(key) != E1Integrals.end())
            {
                overlap = E1Integrals[key];
            }
            else
            {
                const Orbital& p1 = *excited->GetState(e1);
                const Orbital& p2 = *excited->GetState(e2);

                const double* R = excited->GetLattice()->R();
                const double* dR = excited->GetLattice()->dR();
                for(unsigned int x=0; x<mmin(p1.Size(), p2.Size()); x++)
                    overlap += (p1.f[x] * p2.f[x] + p1.g[x] * p2.g[x]) * R[x] * dR[x];

                E1Integrals[key] = overlap;
            }

            matrix_element = coeff * overlap;
        }
    }
    
    *outstream << "jcb me = " << matrix_element << std::endl;
    return matrix_element;
}

// Old M1 calculator code
double RateCalculator::GetM1MatrixElement(const ElectronInfo& e1, const ElectronInfo& e2)
{
    double coeff_j = delta_function(e1.J(), e2.J()) * delta_function(e1.L(), e2.L()) * delta_function(e1.PQN(), e2.PQN()) * pow(e1.J() * (e1.J() + 1.0) * ((2.0 * e1.J()) + 1.0), 0.5);
    coeff_j *= pow(-1.0, e1.J() - e1.M()) * MathConstant::Instance()->Electron3j(e2.J(), e1.J(), 1, e2.M(), -e1.M());
    
    double coeff_s = delta_function(e1.L(), e2.L()) * delta_function(e1.PQN(), e2.PQN()) * pow(-1.0, e1.L() + 0.5 + e1.J() + 1.0);
    coeff_s *= pow(((2.0 * e1.J()) + 1.0) * ((2.0 * e2.J()) + 1.0), 0.5);
    coeff_s *= MathConstant::Instance()->Wigner6j(e1.L(), 0.5, e1.J(), 1, e2.J(), 0.5);
    coeff_s *= pow(0.5 * 1.5 * 2.0, 0.5);
    coeff_s *= pow(-1.0, e1.J() - e1.M()) * MathConstant::Instance()->Electron3j(e2.J(), e1.J(), 1, e2.M(), -e1.M());

    /* Debug
    *outstream << pow(((2.0 * e1.J()) + 1.0) * ((2.0 * e2.J()) + 1.0), 0.5) << " * " << MathConstant::Wigner6j(e1.L(), 0.5, e1.J(), 1, e2.J(), 0.5) << " * " << pow(0.5 * 1.5 * 2.0, 0.5) << std::endl;
    *outstream << "coeff_j + s = " << coeff_j + coeff_s << " = " << coeff_j << " + " << coeff_s << std::endl;
    *outstream << "Electron from " << e1.Name() << " to " << e2.Name() << " with ";
    *outstream << "Reduced Matrix Element = " << (coeff_j + coeff_s) << std::endl;
    *outstream << "Matrix Element: " << (coeff_j + coeff_s) * pow(-1.0, e1.J() - e1.M()) * MathConstant::Instance()->Electron3j(e2.J(), e1.J(), J, e2.M(), -e1.M()) << std::endl;
    */
    
    return (coeff_j + coeff_s);
}

double RateCalculator::CalculateAugerRate(Atom* A, Symmetry sym1, unsigned int solution1, double continuum_energy)
{
    // Assume ion ground state is s_1/2 (for now)
    // Get energy of atom
    Eigenstates* atom_eigenstates = A->GetEigenstates(sym1);
    const double* V1 = atom_eigenstates->GetEigenvectors();
    const double* E1 = atom_eigenstates->GetEigenvalues();
    unsigned int N1 = atom_eigenstates->GetEigenvectorLength();
    RelativisticConfigList* configs1 = atom_eigenstates->GetRelConfigs();

    double cs_energy = E1[solution1] - continuum_energy;
    *outstream << "Continuum energy: " << std::setprecision(6) 
               << cs_energy * MathConstant::Instance()->HartreeEnergyIneV() << std::endl;
    if(cs_energy <= 0.0)
        return 0.0;

    unsigned int final_ion_pqn = 2;

    Core* core = A->GetCore();
    core->ToggleOpenShellCore();

    // Decide how to make the continuum wavefunctions
//    ContinuumBuilder cs_builder;
//    //    cs_builder.CreateNewLattice(1000, 1.e-5, 50.);
//    cs_builder.CopyLattice(core->GetLattice());
//    //    cs_builder.CopyCore(core, false);
//    cs_builder.CreateNewCore(6, 5);

    ContinuumBuilder cs_builder(A->GetCore());
    cs_builder.SetNormalisationType(Unitary);
//    Core* cs_core = cs_builder.GetCore();
//    cs_core->Ionise(OrbitalInfo(1, -1));

//    core->ToggleClosedShellCore();

    // Match parity
    Parity cs_parity;
    // Parity of ion + continuum state == Parity of excited state
    //if(gs.GetParity() == sym1.GetParity())
    //    cs_parity = even;
    //else
    //    cs_parity = odd;

    // Parity of s_1/2 ion is even
    cs_parity = sym1.GetParity();

    bool debug = DebugOptions.LogAugerRate();

    if(debug)
    {   *logstream << "\n\nCS Parity = ";
        if(cs_parity == even)
            *logstream << "even" << std::endl;
        else
            *logstream << "odd" << std::endl;
        *logstream << std::setprecision(8);
    }

    double total = 0.0;
 
    for(int cs_twoj = abs(int(sym1.GetTwoJ()) - 1); cs_twoj <= sym1.GetTwoJ() + 1; cs_twoj += 2)
    {
        double partial_total = 0.0;

        // kappa = (j + 1/2) * (-1)^(j + 1/2 + l)
        int cs_kappa = (cs_twoj + 1)/2;
        if(cs_kappa%2 == 1)
            cs_kappa = -cs_kappa;
        if(cs_parity == odd)
            cs_kappa = -cs_kappa;

        // Calculate continuum state
        ContinuumWave* cs = new ContinuumWave(cs_energy, cs_kappa);
        if(debug)
            *logstream << "Continuum SingleParticleWavefunction: " << cs->Name() << std::endl;

        unsigned int loops = cs_builder.CalculateContinuumWave(cs, core->GetLattice());
        //unsigned int loops = cs_builder.ReadContinuumWave(cs, core->GetLattice(), "fort.64", "fort.65");
        if(loops == 0)
        {   *outstream << "Failed to build continuum state: " << cs->Name() << std::endl;
            return 0.;
        }
/*
        // Orthogonalise
        ConstStateIterator ex_it = excited->GetConstStateIterator();
        while(!ex_it.AtEnd())
        {
            const Orbital* other = ex_it.GetState();
            if(other->Kappa() == cs->Kappa() && other->RequiredPQN() < 15)
            {
                double S = cs->Overlap(*other, excited->GetLattice());

                for(unsigned int i=0; i<mmin(other->Size(), cs->Size()); i++)
                {
                    cs->f[i] = cs->f[i] - S * other->f[i];
                    cs->g[i] = cs->g[i] - S * other->g[i];
                    cs->df[i] = cs->df[i] - S * other->df[i];
                    cs->dg[i] = cs->dg[i] - S * other->dg[i];
                }
            }
            ex_it.Next();
        }

        cs->Print("continuum.txt", core->GetLattice());
*/
        // Form ion+continuum wavefunction
        // For j = J - 1/2, ground state has projection +1/2
        // For j = J + 1/2, ground state can have projection +1/2 or -1/2
        int ion_TwoM = -1;
        if(cs_twoj + 1 == sym1.GetTwoJ())
        {   ion_TwoM = 1;
        }

        while(ion_TwoM <= 1)
        {
            // Here again is our ground state = 2s assumption
            Projection p_ion;
            p_ion.Add(ElectronInfo(final_ion_pqn, -1, ion_TwoM));
            
            double ion_coeff = 1.;
            if(cs_twoj == sym1.GetTwoJ() + 1)
            {   if(ion_TwoM == -1)
                    ion_coeff = - sqrt(double(sym1.GetTwoJ() + 1)/double(sym1.GetTwoJ() + 2));
                else // ion_TwoM == 1
                    ion_coeff = 1./sqrt(double(sym1.GetTwoJ() + 2));
            }

            ElectronInfo cs_electron_info(99, cs_kappa, sym1.GetTwoJ() - ion_TwoM);

            if(debug)
                *logstream << "CS twoM = " << cs_electron_info.TwoM() << "\n"
                           << "Ion twoM = " << ion_TwoM << std::endl;

            // Loop over configurations
            unsigned int i=0;
            RelativisticConfigList::const_iterator list_it = configs1->begin();
            while(list_it != configs1->end())
            {
                const ProjectionSet& proj_i = list_it->GetProjections();
                unsigned int proj_i_size = proj_i.size();
                unsigned int num_states_i = list_it->NumJStates();
                const double* coefficients_i = list_it->GetJCoefficients();
                if(debug)
                    *logstream << "\nproj_i_size = " << proj_i_size
                               << "\nnum_states_i = " << num_states_i << std::endl;

                // Iterate over projections
                ProjectionSet::const_iterator pi_it = proj_i.begin();
                unsigned int pi = 0;
                while(pi_it != proj_i.end())
                {
                    double matrix_element = GetProjectionH(p_ion, *pi_it, cs, cs_electron_info);
                    if(debug)
                        *logstream << "< " << p_ion.Name() << " " << cs_electron_info.Name() << " |H| "
                                   << pi_it->Name() << " >\n    = " << matrix_element << std::endl;

                    if(fabs(matrix_element) > 1.e-15)
                    {
                        double coeff = 0.0;

                        // Summation over jstates
                        for(unsigned int jstate_i = 0; jstate_i < num_states_i; jstate_i++)
                        {
                            coeff += coefficients_i[jstate_i*proj_i_size + pi]
                                     * ion_coeff
                                     * V1[solution1*N1 + i + jstate_i];

                            if(debug)
                                *logstream <<   "      * " << ion_coeff << "  (ion_coeff)"
                                    << "\n      * " << coefficients_i[jstate_i*proj_i_size + pi] << "  (coeff_i)"
                                    << "\n      * " << V1[solution1*N1 + i + jstate_i] << "  (V1)" << std::endl;
                        }

                        if(debug)
                            *logstream << "    = " << coeff*matrix_element << "\n" << std::endl;
                        partial_total += coeff * matrix_element;
                    }

                    pi_it++;
                    pi++;
                }

                list_it++; i+=num_states_i;
            }

            ion_TwoM += 2;
        }

        total += partial_total*partial_total;
        if(debug)
            *logstream << " Total += " << partial_total*partial_total << "\n" << std::endl;
    }

    switch(cs_builder.GetNormalisationType())
    {
        case Unitary:
            // Total using unitary normalisation of continuum wavefunctions.
            total = total * 4.0/sqrt(2.0 * cs_energy);
            break;

        case Cowan:
            // Total using Cowan normalisation of continuum wavefunctions.
            total = total * 2. * MathConstant::Instance()->Pi();
            break;
            
        default:
            *errstream << "Normalisation type not supported: " << cs_builder.GetNormalisationType();
            exit(1);
            break;
    }

    // Convert to inverse seconds (divide by hbar)
    total = total/2.418884e-17;

    *outstream << "Auger rate (/sec) = " << total << std::endl;

    return total;
}

/** Get the Hamiltonian matrix element <first + continuum | H | second> */
double RateCalculator::GetProjectionH(const Projection& first, const Projection& second, const ContinuumWave* cs, const ElectronInfo& cs_electron) const
{
    unsigned int diff[4];   // Storage for projection differences.
    int numdiff = Projection::GetProjectionDifferences(first, second, diff);

    int sign;
    if(numdiff >= 0)
        sign = 1;
    else
        sign = -1;

    double value = 0.;

    if(numdiff == 0)
    {   // Never happen unless forced
        const ElectronInfo& f1 = first[0];
        const ElectronInfo& s1 = second[1];
        const ElectronInfo& s2 = second[0];
        
        value = CoulombMatrixElement(cs_electron, f1, s2, s1, cs, sign)
                    - CoulombMatrixElement(cs_electron, f1, s1, s2, cs, -sign);
    }

    if(abs(numdiff) == 1)
    {
        if(diff[0] != first.Size())
        {   *errstream << "RateCalculator::GetProjectionH: diff doesn't include continuum state!" << std::endl;
            exit(1);
        }

        StateIntegrator SI(excited->GetLattice());
        const ElectronInfo& s1 = second[diff[1]];

        // <a|f|b>
        // Note: SI.HamiltonianMatrixElement uses derivative of second orbital, so put
        //       continuum first, bound second.
//        if((cs_electron.M() == s1.M()) && (cs_electron.Kappa() == s1.Kappa()))
//        {
//            const Orbital* sj = excited->GetState(s1);
//            value = SI.HamiltonianMatrixElement(*cs, *sj, *excited->GetCore()) * sign;
////            value = SI.HamiltonianMatrixElement(*sj, *cs, *excited->GetCore()) * sign;
//
//            if(DebugOptions.LogAugerRate())
//                *logstream << "\t < " << cs_electron.Name() << " | H(1) | " << s1.Name() << " >\n"
//                    << "\t     = " << sign << " * " << value*sign << " = " << value << std::endl;
//        }

        // Sum(e) <ae|g|be> - <ae|g|eb>
        for(unsigned int i=0; i<first.Size(); i++)
        {
            if(i != diff[0])
            {   const ElectronInfo& e(first[i]);
                value += sign * (CoulombMatrixElement(cs_electron, e, s1, e, cs, sign) 
                                - CoulombMatrixElement(cs_electron, e, e, s1, cs, -sign));
            }
        }
    }
    else if(abs(numdiff) == 2)
    {
        if(diff[2] != first.Size())
        {   *errstream << "RateCalculator::GetProjectionH: diff doesn't include continuum state!" << std::endl;
            exit(1);
        }

        const ElectronInfo& f1 = first[diff[0]];
        const ElectronInfo& s1 = second[diff[1]];
        const ElectronInfo& s2 = second[diff[3]];

        // a->b, c->d
        // <ac|g|bd> - <ac|g|db>
        value = sign * (CoulombMatrixElement(cs_electron, f1, s2, s1, cs, sign)
                        - CoulombMatrixElement(cs_electron, f1, s1, s2, cs, -sign));
    }

    return value;
}

/** Get the Coulomb matrix element < e1, e2 | 1/r | e3, e4 >.
    e1 is the continuum state.
 */
double RateCalculator::CoulombMatrixElement(const ElectronInfo& e1, const ElectronInfo& e2, const ElectronInfo& e3, const ElectronInfo& e4, const ContinuumWave* cs, int sign) const
{
    // Get two-body matrix element
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

    bool debug = DebugOptions.LogAugerRate();
    double total = 0.;

    // Prepare radial matrix element
    const SingleParticleWavefunction* s_1;
    const SingleParticleWavefunction* s_2 = excited->GetState(e2);
    const SingleParticleWavefunction* s_3 = excited->GetState(e3);
    const SingleParticleWavefunction* s_4 = excited->GetState(e4);
    
    if(!s_2 || !s_3 || !s_4)
    {   *errstream << "RateCalculator::CoulombMatrixElement couldn't find state." << std::endl;
        exit(1);
    }

    unsigned int p;
    CoulombIntegrator CI(excited->GetLattice());
    const double* R = excited->GetLattice()->R();
    const double* dR = excited->GetLattice()->dR();
    const double core_pol = excited->GetCore()->GetPolarisability();
    const double core_rad = excited->GetCore()->GetClosedShellRadius();

    // Get density24
    std::vector<double> density(mmin(s_2->Size(), s_4->Size()));
    for(p=0; p<density.size(); p++)
    {
        density[p] = s_2->f[p] * s_4->f[p] + s_2->g[p] * s_4->g[p];
    }
    density.resize(excited->GetCore()->GetHFPotential().size());

    // Orthogonalise continuum wave
    ContinuumWave* cs_orth = NULL;

    if(cs->Kappa() == s_3->Kappa())
    {
        cs_orth = new ContinuumWave(*cs);
        double S = cs_orth->Overlap(*s_3, excited->GetLattice());

        for(unsigned int i=0; i<mmin(cs_orth->Size(), s_3->Size()); i++)
        {
            cs_orth->f[i] = cs_orth->f[i] - S * s_3->f[i];
            cs_orth->g[i] = cs_orth->g[i] - S * s_3->g[i];
            cs_orth->df[i] = cs_orth->df[i] - S * s_3->df[i];
            cs_orth->dg[i] = cs_orth->dg[i] - S * s_3->dg[i];
        }
        
        s_1 = cs_orth;
    }
    else
        s_1 = cs;

    // Calculate diagram
    while(k <= kmax)
    {
        double coeff = 0.;
        if(fabs(q) <= k)
            coeff = MathConstant::Instance()->Electron3j(e1.TwoJ(), e3.TwoJ(), k, -e1.TwoM(), e3.TwoM()) *
                    MathConstant::Instance()->Electron3j(e2.TwoJ(), e4.TwoJ(), k, -e2.TwoM(), e4.TwoM());
            
        if(coeff)
            coeff = coeff * MathConstant::Instance()->Electron3j(e1.TwoJ(), e3.TwoJ(), k, 1, -1) *
                            MathConstant::Instance()->Electron3j(e2.TwoJ(), e4.TwoJ(), k, 1, -1);

        if(coeff)
        {
            if(int(q - e1.M() - e2.M() + 1.)%2)
                coeff = - coeff;

            coeff = coeff * sqrt(double(e1.MaxNumElectrons() * e2.MaxNumElectrons() *
                                        e3.MaxNumElectrons() * e4.MaxNumElectrons()));

            //integrals.GetTwoElectronIntegral(k, e1, e2, e3, e4)
            double radial = 0.;

            // Get Pot24
            std::vector<double> Pot24(density.size());
            CI.FastCoulombIntegrate(density, Pot24, k);

            unsigned int limit = mmin(s_1->Size(), s_3->Size());
            limit = mmin(limit, Pot24.size());
            for(p=0; p<limit; p++)
            {
                radial += (s_1->f[p] * s_3->f[p] + s_1->g[p] * s_3->g[p])
                            * Pot24[p] * dR[p];
            }
//            for(p=0; p<limit; p+=2)
//            {
//                radial += 4.* (s_1->f[p] * s_3->f[p] + MathConstant::AlphaSquared * s_1->g[p] * s_3->g[p])
//                            * Pot24[p] * dR[p];
//            }
//            for(p=1; p<limit; p+=2)
//            {
//                radial += 2.* (s_1->f[p] * s_3->f[p] + MathConstant::AlphaSquared * s_1->g[p] * s_3->g[p])
//                            * Pot24[p] * dR[p];
//            }
//            radial = radial/3.;

            if(core_pol && k == 1)
            {
                double R1 = 0.;
                double R2 = 0.;
                for(p=0; p<limit; p++)
                {
                    double r2 = R[p]*R[p] + core_rad*core_rad;
                    R1 += (s_1->f[p] * s_3->f[p] + s_1->g[p] * s_3->g[p])/r2 * dR[p];
                    R2 += density[p]/r2 * dR[p];
                }

                radial -= core_pol * R1 * R2;
            }

            if(debug)
                *logstream << "\t < " << e1.Name() << e2.Name() << " | R(k=" << k << ") | "
                    << e3.Name() << e4.Name() << " >\n"
                    << "\t     = " << sign*coeff << " * " << radial << " = " << sign*coeff * radial << std::endl;

            //radial += RateCalculator::SubtractionDiagram(cs, s_2, s_3, s_4, k);
            total += coeff * radial;
        }

        k = k+2;
    }

    if(cs_orth)
        delete cs_orth;

    return total;
}

/** Get the subtraction diagram.
 */
double RateCalculator::SubtractionDiagram(const ContinuumWave* cs, const SingleParticleWavefunction* sb, const SingleParticleWavefunction* sc, const SingleParticleWavefunction* sd, unsigned int k) const
{
    bool debug = DebugOptions.LogAugerRate();
    double total = 0.;

    if(!cs || !sb || !sc || !sd)
    {   *errstream << "RateCalculator::SubtractionDiagram couldn't find state." << std::endl;
        exit(1);
    }

    if(debug)
        *logstream << "\tSubtraction diagrams:" << std::endl;

    unsigned int p;
    StateIntegrator SI(excited->GetLattice());
    CoulombIntegrator CI(excited->GetLattice());
    const double* R = excited->GetLattice()->R();
    const double* dR = excited->GetLattice()->dR();
    const double core_pol = excited->GetCore()->GetPolarisability();
    const double core_rad = excited->GetCore()->GetClosedShellRadius();

    // Get density24
    std::vector<double> density(mmin(sb->Size(), sd->Size()));
    for(p=0; p<density.size(); p++)
    {
        density[p] = sb->f[p] * sd->f[p] + sb->g[p] * sd->g[p];
    }
    density.resize(excited->GetCore()->GetHFPotential().size());

    // Get Pot24
    std::vector<double> Pot24(density.size());
    CI.FastCoulombIntegrate(density, Pot24, k);

    // Loop over excited states
    ConstStateIterator it = excited->GetConstStateIterator();

    while(!it.AtEnd())
    {   const SingleParticleWavefunction* sa = it.GetState();

        if(cs->Kappa() == sa->Kappa())
        {
            double radial = 0.;

            unsigned int limit = mmin(sa->Size(), sc->Size());
            limit = mmin(limit, Pot24.size());
            for(p=0; p<limit; p++)
            {
                radial += (sa->f[p] * sc->f[p] + sa->g[p] * sc->g[p])
                            * Pot24[p] * dR[p];
            }

            if(debug)
                *logstream << "\t< " << cs->Name() << " |h| " << sa->Name() << "> = "
                           << cs->Overlap(*sa, excited->GetLattice())
//                           << SI.HamiltonianMatrixElement(*cs, *sa, *excited->GetCore())
                           << " * " << radial << " / " << (cs->Energy() - sa->Energy()) << std::endl;

            total -= radial
                      * cs->Overlap(*sa, excited->GetLattice());
//                      * SI.HamiltonianMatrixElement(*cs, *sa, *excited->GetCore())
//                      / (cs->Energy() - sa->Energy());
            
            if(debug)
                *logstream << "\t = " << total << std::endl;
        }

        it.Next();
    }

    return total;
}

double RateCalculator::GetEJMatrixElement(unsigned int J, const ElectronInfo& e1, const ElectronInfo& e2, TransitionGaugeType::Enum aGaugeType)
{
    //*outstream << "Electron from " << e1.Name() << " to " << e2.Name() << std::endl;
    
    double matrix_element = 0.0;
    double coeff = 0.0;
    
    // This is the reduced matrix element of the spherical tensor CJ
    /*
    coeff = MathConstant::Wigner3j(e1.L(), J, e2.L(), 0.0, 0.0, 0.0);
    coeff *= pow(-1.0, e1.L()) * pow((2.0 * e1.L()) + 1.0, 0.5) * pow((2.0 * e2.L()) + 1.0, 0.5);

    // This factor is (-1)^(L + S + J' + k) [J, J']^(1/2) {L  S  J }
    //                                                    {J' k  L'}
    // Because the operator depends only on L (kappa)/J
    coeff *= pow(-1.0, e1.L() + e2.J() + 0.5 + J) * pow((2.0 * e1.J()) + 1.0, 0.5) * pow((2.0 * e2.J()) + 1.0, 0.5) * MathConstant::Wigner6j(e1.L(), 0.5, e1.J(), e2.J(), J, e2.L());
    
    *outstream << "Coeff = " << coeff << " | Formula for same thing = " << SphericalTensorReducedMatrixElement(J, e1.Kappa(), e2.Kappa()) << std::endl;*/
    
    coeff = SphericalTensorReducedMatrixElement(J, e1.Kappa(), e2.Kappa());
    // Convert from reduced matrix element to full matrix element
    coeff *= pow(-1.0, e1.J() - e1.M()) * MathConstant::Instance()->Electron3j(e2.J(), e1.J(), J, e2.M(), -e1.M());
    
    // Don't bother computing the overlap if the angular part is zero
    if(coeff == 0)
    {
        return 0.0;
    }
    
    double overlap = 0.;
    
    const Orbital& p1 = *excited->GetState(e1);
    const Orbital& p2 = *excited->GetState(e2);
    const double* R = excited->GetLattice()->R();
    const double* dR = excited->GetLattice()->dR();
    const double alpha = PhysicalConstant::Instance()->GetAlpha();
    
    // k = omega/c, in atomic units
    // omega = Energy_1 - Energy_2 also in atomic units
    double k = (double) fabs((double) (p1.Energy() - p2.Energy())) * alpha;
    
    if(aGaugeType == TransitionGaugeType::Length)
    {
        if(k == 0)
        {
            // For the same state to itself, a special case is needed to compute the
            // overlap, because dependence of the spherical bessel function on k
            // exactly cancels the k dependence added by reintroducing dimensionality
            for(unsigned int x=0; x<mmin(p1.Size(), p2.Size()); x++)
            {
                overlap += ((p1.f[x] * p2.f[x]) + (p1.g[x] * p2.g[x])) * sph_bessel_small_limit(J, R[x]) * dR[x];
            }
            
            matrix_element = coeff * overlap;
            matrix_element = ((double) boost::math::double_factorial<double>((2.0 * ((double) J)) + 1.0)) * matrix_element;
        }
        else
        {
            for(unsigned int x=0; x<mmin(p1.Size(), p2.Size()); x++)
            {
                overlap += ((p1.f[x] * p2.f[x]) + (p1.g[x] * p2.g[x])) * boost::math::sph_bessel(J, R[x] * k) * dR[x];
                overlap += ((p1.f[x] * p2.g[x]) + (p1.g[x] * p2.f[x])) * alpha * boost::math::sph_bessel(J + 1, R[x] * k) * ((((double) e1.Kappa()) - ((double) e2.Kappa()))/(((double) J) + 1.0)) * dR[x];
                overlap += ((p1.f[x] * p2.g[x]) - (p1.g[x] * p2.f[x])) * alpha * boost::math::sph_bessel(J + 1, R[x] * k) * dR[x];
            }
            
            matrix_element = coeff * overlap;
            // This is the matrix operator with dimensions in atomic units
            matrix_element = ((double) boost::math::double_factorial<double>((2.0 * ((double) J)) + 1.0)/pow((double) k, J)) * matrix_element;
            //*outstream << "Overlap = " << overlap /pow((double) k, J) * (double) boost::math::double_factorial<double>((2.0 * ((double) J)) + 1.0) << std::endl;

            //*outstream << "Reduced matrix element = " << matrix_element/(pow(-1.0, e1.J() - e1.M()) * MathConstant::Instance()->Electron3j(e2.J(), e1.J(), J, e2.M(), -e1.M())) << std::endl;
        }
    }
    else if(aGaugeType == TransitionGaugeType::Velocity)
    {
        if(k == 0)
        {
            // For the same state to itself, a special case is needed to compute the
            // overlap, because dependence of the spherical bessel function on k
            // exactly cancels the k dependence added by reintroducing dimensionality
            for(unsigned int x=0; x<mmin(p1.Size(), p2.Size()); x++)
            {
                //overlap += ((p1.f[x] * p2.f[x]) + (p1.g[x] * p2.g[x])) * sph_bessel_small_limit(J, R[x]) * dR[x];
            }
            
            matrix_element = ((double) boost::math::double_factorial<double>((2.0 * ((double) J)) + 1.0)) * matrix_element;
        }
        else
        {
            for(unsigned int x=0; x<mmin(p1.Size(), p2.Size()); x++)
            {
                overlap += -1.0 * ((((double) e1.Kappa()) - ((double) e2.Kappa()))/(((double) J) + 1.0)) * (sph_bessel_prime(J, k * R[x]) + (boost::math::sph_bessel(J, k * R[x])/(k *R[x]))) * ((p1.f[x] * p2.g[x]) + (p1.g[x] * p2.f[x])) * alpha * dR[x];
                overlap += ((p1.f[x] * p2.g[x]) - (p1.g[x] * p2.f[x])) * J * (boost::math::sph_bessel(J, k * R[x])/(k * R[x])) * alpha * dR[x];
            }
            
            double overlap1 = 0.0;
            double overlap2 = 0.0;
            for(unsigned int x=0; x<mmin(p1.Size(), p2.Size()); x++)
            {
                overlap1 += -1.0 * ((((double) e1.Kappa()) - ((double) e2.Kappa()))/(((double) J) + 1.0)) * (sph_bessel_prime(J, k * R[x]) + (boost::math::sph_bessel(J, k * R[x])/(k *R[x]))) * ((p1.f[x] * p2.g[x]) + (p1.g[x] * p2.f[x])) * alpha * dR[x];
                overlap2 += ((p1.f[x] * p2.g[x]) - (p1.g[x] * p2.f[x])) * J * (boost::math::sph_bessel(J, k * R[x])/(k * R[x])) * alpha * dR[x];
            }
            
            matrix_element = coeff * overlap;
            // This is the matrix operator with dimensions in atomic units
            matrix_element = ((double) boost::math::double_factorial<double>((2.0 * ((double) J)) + 1.0)/pow((double) k, J)) * matrix_element;
            
            //*outstream << "Overlap = " << overlap << std::endl;
            //*outstream << "Overlap1 = " << overlap1 /pow((double) k, J) * (double) boost::math::double_factorial<double>((2.0 * ((double) J)) + 1.0) << std::endl;
            //*outstream << "Overlap2 = " << overlap2 /pow((double) k, J) * (double) boost::math::double_factorial<double>((2.0 * ((double) J)) + 1.0) << std::endl;
            //*outstream << "Overlap = " << overlap /pow((double) k, J) * (double) boost::math::double_factorial<double>((2.0 * ((double) J)) + 1.0) << std::endl;
            
            //*outstream << "Reduced matrix element = " << matrix_element/(pow(-1.0, e1.J() - e1.M()) * MathConstant::Instance()->Electron3j(e2.J(), e1.J(), J, e2.M(), -e1.M())) << std::endl;
        }
    }
    
    //*outstream << "matrix_element (dim) = " << matrix_element << std::endl;

    return matrix_element;
}

std::complex<double> RateCalculator::GetEJMatrixElementZenonas(unsigned int J, const ElectronInfo& e1, const ElectronInfo& e2, TransitionGaugeType::Enum aGaugeType)
{
    //*outstream << "Electron from " << e1.Name() << " to " << e2.Name() << std::endl;
    
    std::complex<double> coeff = 0.0;
    double matrix_element = 0.0;
    
    coeff = pow((std::complex<double>) -1.0, (std::complex<double>) (double(e2.L()) - double(J) - double(e1.L()) - 1.0)/(2.0 - double(e1.J())));
    //coeff *= pow(((2.0 * double(e1.J())) + 1.0) * ((2.0 * double(e2.J())) + 1.0) / ((2.0 * double(J)) + 1.0), 0.5);
    coeff *= pow(((2.0 * double(e1.J())) + 1.0) * ((2.0 * double(e2.J())) + 1.0), 0.5);
    coeff *= MathConstant::Instance()->Wigner3j(e2.J(), e1.J(), double(J), 0.5, -0.5, 0.0);
    //*outstream << "Wigner3j = " << MathConstant::Wigner3j(e2.J(), e1.J(), double(J), 0.5, -0.5, 0.0) << std::endl;
    //*outstream << "(" << e2.J() << " " << e1.J() <<  " " << double(J) << ") = " << MathConstant::Wigner3j(e2.J(), e1.J(), double(J), 0.5, -0.5, 0.0) << std::endl;
    //*outstream << "(" << 0.5 << " " << -0.5 <<  " " << 0 << ")" << std::endl;
    coeff *= pow(-1.0, e1.J() - e1.M()) * MathConstant::Instance()->Electron3j(e2.J(), e1.J(), J, e2.M(), -e1.M());
    coeff *= pow(-1.0, e1.L() + e2.J() + 0.5 + J);
    
    Integrator integrator(excited->GetLattice());
    const Orbital& p1 = *excited->GetState(e1);
    const Orbital& p2 = *excited->GetState(e2);
    double k = fabs(p1.Energy() - p2.Energy()) * PhysicalConstant::Instance()->GetAlpha();
    
    if(triangular_condition_even(e1.L(), e2.L(), J) == 1.0)
    {
        matrix_element += integrator.AIntegralUpper(p2, p1, J, k) - (((double(2.0 * J) + 1.0)/(double(J) + 1.0)) * integrator.BIntegralUpper(p2, p1, J, k));
        //*outstream << "Triangular 1 (" << e1.L() << ", " << e2.L() << ", " << J << ") = " << integrator.AIntegralUpper(p2, p1, J, k) - (((double(2.0 * J) + 1.0)/(double(J) + 1.0)) * integrator.BIntegralUpper(p2, p1, J, k)) << std::endl;
    }

    if(triangular_condition_even(e1.L_Prime(), e2.L_Prime(), J) == 1.0)
    {
        matrix_element += integrator.AIntegralLower(p2, p1, J, k) + (((double(2.0 * J) + 1.0)/(double(J) + 1.0)) * integrator.BIntegralLower(p2, p1, J, k));
        //*outstream << "Triangular 2 (" << e1.L_Prime() << ", " << e2.L_Prime() << ", " << J << ") = " << integrator.AIntegralLower(p2, p1, J, k) + (((double(2.0 * J) + 1.0)/(double(J) + 1.0)) * integrator.BIntegralLower(p2, p1, J, k)) << std::endl;
    }

    //*outstream << matrix_element * coeff << " = " << (std::real(matrix_element * coeff)/fabs(std::real(matrix_element * coeff))) * std::abs(matrix_element * coeff) << std::endl;
    if(std::real(matrix_element * coeff) == 0.0 && std::imag(matrix_element * coeff) == 0.0)
    {
        return 0.0;
    }
    
    return (std::real(matrix_element * coeff)/fabs(std::real(matrix_element * coeff))) * std::abs(matrix_element * coeff);
}

double RateCalculator::GetMJMatrixElement(Atom* A, unsigned int J, const ElectronInfo& e1, const ElectronInfo& e2)
{   
    //*outstream << "Electron from " << e1.Name() << " to " << e2.Name() << std::endl;
    double matrix_element = 0.0;
    double coeff = 0.0;
    
    // Determine the effective L of the transition
    // L initial is replaced by L initial - 1 if Kappa > 0
    // otherwise it is replaced by L initial + 1
    // Note in the end this function must be symmetric on exchange of e1 and e2
    /*
    if(e1.Kappa() > 0)
    {
        // This is the reduced matrix element of the spherical tensor CJ
        coeff = MathConstant::Wigner3j((e1.L() - 1.0), J, e2.L(), 0.0, 0.0, 0.0);
        coeff *= pow(-1.0, e1.L() - 1.0) * pow((2.0 * (e1.L() - 1.0)) + 1.0, 0.5) * pow((2.0 * e2.L()) + 1.0, 0.5);
        
        // This factor is (-1)^(L + S + J' + k) [J, J']^(1/2) {L  S  J }
        //                                                    {J' k  L'}
        // Because the operator depends only on L (kappa)/J
        coeff *= pow(-1.0, (e1.L() - 1.0) + e2.J() + J + 0.5) * pow((2.0 * e1.J()) + 1.0, 0.5) * pow((2.0 * e2.J()) + 1.0, 0.5) * MathConstant::Wigner6j((e1.L() - 1.0), 0.5, e1.J(), e2.J(), J, e2.L());
    }
    else if(e1.Kappa() < 0)
    {
        // This is the reduced matrix element of the spherical tensor CJ
        coeff = MathConstant::Wigner3j((e1.L() + 1.0), J, e2.L(), 0.0, 0.0, 0.0);
        coeff *= pow(-1.0, e1.L() + 1.0) * pow((2.0 * (e1.L() + 1.0)) + 1.0, 0.5) * pow((2.0 * e2.L()) + 1.0, 0.5);
        
        // This factor is (-1)^(L + S + J' + k) [J, J']^(1/2) {L  S  J }
        //                                                    {J' k  L'}
        // Because the operator depends only on L (kappa)/J
        coeff *= pow(-1.0, (e1.L() + 1.0) + e2.J() + J + 0.5) * pow((2.0 * e1.J()) + 1.0, 0.5) * pow((2.0 * e2.J()) + 1.0, 0.5) * MathConstant::Wigner6j((e1.L() + 1.0), 0.5, e1.J(), e2.J(), J, e2.L());
    }
    
    *outstream << "Coeff = " << coeff << " | Formula for same thing = " << SphericalTensorReducedMatrixElement(J, -e1.Kappa(), e2.Kappa()) << std::endl;*/
    
    coeff = SphericalTensorReducedMatrixElement(J, -e1.Kappa(), e2.Kappa());
    // Convert from reduced matrix element to full matrix element
    coeff *= pow(-1.0, e1.J() - e1.M()) * MathConstant::Instance()->Electron3j(e2.J(), e1.J(), J, e2.M(), -e1.M());
    
    // Don't bother computing the overlap if the angular part is zero
    if(coeff == 0)
    {
        return 0.0;
    }

    double overlap = 0.;
    
    const Orbital& p1 = *excited->GetState(e1);
    const Orbital& p2 = *excited->GetState(e2);
    const double* R = excited->GetLattice()->R();
    const double* dR = excited->GetLattice()->dR();
    const double alpha = PhysicalConstant::Instance()->GetAlpha();
    
    // k = omega/c, in atomic units
    // omega = Energy_1 - Energy_2 also in atomic units
    double k = fabs(p1.Energy() - p2.Energy()) * alpha;
    
    if(k == 0)
    {
        // For the same state to itself, a special case is needed to compute the
        // overlap, because dependence of the spherical bessel function on k
        // exactly cancels the k dependence added by reintroducing dimensionality
        for(unsigned int x=0; x<mmin(p1.Size(), p2.Size()); x++)
        {
            overlap += ((p1.f[x] * p2.g[x]) + (p1.g[x] * p2.f[x])) * alpha * sph_bessel_small_limit(J, R[x]) * ((((double) e1.Kappa()) + ((double) e2.Kappa()))/(((double) J) + 1.0)) * dR[x];
        }
        
        matrix_element = coeff * overlap;
        
        // This is the matrix operator with dimensions in atomic units
        matrix_element = boost::math::double_factorial<double>((2.0 * ((double) J)) + 1.0) * matrix_element * 2.0 * (1.0 / alpha);
        //*outstream << "k = 0, Matrix element (dim) = " << matrix_element << std::endl;
    }
    else
    {
        for(unsigned int x=0; x<mmin(p1.Size(), p2.Size()); x++)
        {
            overlap += ((p1.f[x] * p2.g[x]) + (p1.g[x] * p2.f[x])) * alpha * boost::math::sph_bessel(J, R[x] * k) * ((((double) e1.Kappa()) + ((double) e2.Kappa()))/(((double) J) + 1.0)) * dR[x];
        }
        
        matrix_element = coeff * overlap;
        
        // This is the matrix operator with dimensions in atomic units
        matrix_element = (boost::math::double_factorial<double>((2.0 * ((double) J)) + 1.0)/pow((double) k, J)) * matrix_element * 2.0 * (1.0 / alpha);
        //*outstream << "Matrix element (dim) = " << matrix_element << std::endl;
    }

    return matrix_element;
}
#include "RateCalculator.h"
#include "Include.h"
#include "Universal/Constant.h"
#include "HartreeFock/StateIntegrator.h"
#include "Universal/CoulombIntegrator.h"

void RateCalculator::CalculateAllDipoleStrengths(Atom* A, Symmetry sym1, unsigned int solution1)
{
    Eigenstates* eigenstates1 = A->GetEigenstates(sym1);
    if(!eigenstates1 || !eigenstates1->Restore())
    {   *outstream << sym1.GetString() << "eigenstates restore failed" << std::endl;
        return;
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

    for(unsigned int TwoJ = TwoJ_min; TwoJ <= TwoJ_max; TwoJ += 2)
    {
        Symmetry sym2(TwoJ, other_P);
        Eigenstates* eigenstates2 = A->GetEigenstates(sym2);
        if(!eigenstates2 || !eigenstates2->Restore())
        {   *outstream << sym2.GetString() << "eigenstates restore failed" << std::endl;
            return;
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
        
        *outstream << "(" << sym1.GetString() << ": " << solution1 << ") -> (" << sym2.GetString()
                   << ": i) oscillator strengths" << std::endl;
        double wigner_coeff = Constant::Electron3j(sym1.GetTwoJ(), TwoJ, 1, -sym1.GetTwoJ(), TwoJ);
        for(solution2 = 0; solution2 < num_solutions2; solution2++)
        {
            // Get reduced matrix element squared
            total[solution2] = total[solution2]/wigner_coeff;
            total[solution2] = total[solution2] * total[solution2];
            
            // Oscillator strength
            double deltaE = fabs(E1[solution1] - E2[solution2]);
            *outstream << "  " << solution2 << "  " << deltaE * total[solution2] << std::endl;
        }
        
        delete[] coeff;
        delete[] total;
    }
}

double RateCalculator::CalculateDipoleStrength(Atom* A, Symmetry sym1, unsigned int solution1, Symmetry sym2, unsigned int solution2)
{
    Eigenstates* eigenstates1 = A->GetEigenstates(sym1);
    if(!eigenstates1 || !eigenstates1->Restore())
    {   *errstream << "eigenstates1 restore failed" << std::endl;
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
    unsigned int num_electrons = configs1->front().NumParticles();
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
                            }
                        }
                        
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

    total = total/Constant::Electron3j(sym1.GetTwoJ(), sym2.GetTwoJ(), 1, -sym1.GetTwoJ(), sym2.GetTwoJ());
    total = total * total;
    *outstream << "E1 reduced matrix element squared = " << total << std::endl;
    
    double deltaE = fabs((eigenstates1->GetEigenvalues())[solution1] - (eigenstates2->GetEigenvalues())[solution2]);
    *outstream << "Oscillator strength f = " << deltaE * total << std::endl;

    double sig = deltaE * Constant::HartreeEnergy_cm;
    *outstream << "sigma = " << sig <<std::endl;
    *outstream << "Radiative rate (/sec) = " << total * sig * sig * sig * 2.0261e-6 << std::endl;

    return total;
}

double RateCalculator::GetE1MatrixElement(const ElectronInfo& e1, const ElectronInfo& e2) const
{
    double matrix_element = 0.0;

    // Check e1.L() + e2.L() + 1 is even
    if((e1.L() + e2.L())%2)
    {
        double coeff = Constant::Electron3j(e1.TwoJ(), e2.TwoJ(), 1, -e1.TwoM(), e2.TwoM());
        
        if(coeff)
        {   coeff = coeff * Constant::Electron3j(e1.TwoJ(), e2.TwoJ(), 1, 1, -1)
                          * sqrt(double(e1.MaxNumElectrons() * e2.MaxNumElectrons()));

            if((abs(e1.TwoM() + 1)/2)%2)
                coeff = -coeff;

            const DiscreteState& p1 = *excited->GetState(e1);
            const DiscreteState& p2 = *excited->GetState(e2);

            double overlap = 0.;
            const double* R = excited->GetLattice()->R();
            const double* dR = excited->GetLattice()->dR();
            for(unsigned int x=0; x<mmin(p1.Size(), p2.Size()); x++)
                overlap += (p1.f[x] * p2.f[x] + Constant::AlphaSquared * p1.g[x] * p2.g[x]) * R[x] * dR[x];

            matrix_element = coeff * overlap;
        }
    }
    
    return matrix_element;
}

double RateCalculator::CalculateAugerRate(Atom* ion, Symmetry gs, Atom* atom, Symmetry sym1, unsigned int solution1)
{
    // Assume gs is s_1/2 (for now)
    // Get energy of atom
    Eigenstates* ion_eigenstates = ion->GetEigenstates(gs);
    Eigenstates* atom_eigenstates = atom->GetEigenstates(sym1);
    double continuum_energy = atom_eigenstates->GetEigenvalues()[solution1]
                              - ion_eigenstates->GetEigenvalues()[0];

    const double* V1 = atom_eigenstates->GetEigenvectors();
    const double* E1 = atom_eigenstates->GetEigenvalues();
    unsigned int N1 = atom_eigenstates->GetEigenvectorLength();
    RelativisticConfigList* configs1 = atom_eigenstates->GetRelConfigs();

    Core* core = ion->GetCore();

    // Match parity
    Parity cs_parity;
    if(gs.GetParity() == sym1.GetParity())
        cs_parity = even;
    else
        cs_parity = odd;

    double total = 0.0;

    // change to step size 2!
    for(int cs_twoj = abs(int(sym1.GetTwoJ()) - 1); cs_twoj <= sym1.GetTwoJ() + 1; cs_twoj += 2)
    {
        // kappa = (j + 1/2) * (-1)^(j + 1/2 + l)
        int cs_kappa = (cs_twoj + 1)/2;
        if(cs_kappa%2 == 1)
            cs_kappa = -cs_kappa;
        if(cs_parity = odd)
            cs_kappa = -cs_kappa;

        // Calculate continuum state
        ContinuumState* cs = new ContinuumState(continuum_energy, cs_kappa);
        core->CalculateExcitedState(cs);
        //core->CalculateContinuumState(cs);

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
            p_ion.Add(ElectronInfo(2, -1, ion_TwoM));
            
            double ion_coeff = 1.;
            if(cs_twoj == sym1.GetTwoJ() + 1)
            {   if(ion_TwoM == -1)
                    ion_coeff = - sqrt(double(sym1.GetTwoJ() + 1)/double(sym1.GetTwoJ() + 2));
                else // ion_TwoM == 1
                    ion_coeff = 1./sqrt(double(sym1.GetTwoJ() + 2));
            }
            
            ElectronInfo cs_electron_info(99, cs_kappa, sym1.GetTwoJ() - ion_TwoM);

            // Loop over configurations
            unsigned int i=0, j;
            RelativisticConfigList::const_iterator list_it = configs1->begin();
            while(list_it != configs1->end())
            {
                const ProjectionSet& proj_i = list_it->GetProjections();
                unsigned int proj_i_size = proj_i.size();
                unsigned int num_states_i = list_it->NumJStates();
                const double* coefficients_i = list_it->GetJCoefficients();

                // Iterate over projections
                ProjectionSet::const_iterator pi_it = proj_i.begin();
                unsigned int pi = 0;
                while(pi_it != proj_i.end())
                {
                    double matrix_element = GetProjectionH(p_ion, *pi_it, cs, cs_electron_info);

                    if(fabs(matrix_element) > 1.e-16)
                    {
                        double coeff = 0.0;

                        // Summation over jstates
                        for(unsigned int jstate_i = 0; jstate_i < num_states_i; jstate_i++)
                        {
                            coeff += coefficients_i[jstate_i*proj_i_size + pi]
                                     * ion_coeff
                                     * V1[solution1*N1 + i + jstate_i];
                        }

                        total += coeff * matrix_element;
                    }

                    pi_it++;
                }

                list_it++; i+=num_states_i;
            }

            ion_TwoM += 2;
        }
    }

    return total;
}

/** Get the Hamiltonian matrix element <first + continuum | H | second> */
double RateCalculator::GetProjectionH(const Projection& first, const Projection& second, const ContinuumState* cs, const ElectronInfo& cs_electron) const
{
    unsigned int diff[4];   // Storage for projection differences.
    int numdiff = Projection::GetProjectionDifferences(first, second, diff);

    int sign;
    if(numdiff >= 0)
        sign = 1;
    else
        sign = -1;

    double value = 0.;

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
        if((cs_electron.M() == s1.M()) && (cs_electron.Kappa() == s1.Kappa()))
        {
            const DiscreteState* sj = excited->GetState(s1);
            value = SI.HamiltonianMatrixElement(*cs, *sj, *excited->GetCore()) * sign;
        }

        // Sum(e) <ae|g|be> - <ae|g|eb>
        for(unsigned int i=0; i<first.Size(); i++)
        {
            if(i != diff[0])
            {   const ElectronInfo& e(first[i]);
                value += sign * (CoulombMatrixElement(cs_electron, e, s1, e, cs) 
                                - CoulombMatrixElement(cs_electron, e, e, s1, cs));
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
        value = sign * (CoulombMatrixElement(cs_electron, f1, s2, s1, cs)
                        - CoulombMatrixElement(cs_electron, f1, s1, s2, cs));
    }

    return value;
}

/** Get the Coulomb matrix element < e1, e2 | 1/r | e3, e4 >.
    e1 is the continuum state.
 */
double RateCalculator::CoulombMatrixElement(const ElectronInfo& e1, const ElectronInfo& e2, const ElectronInfo& e3, const ElectronInfo& e4, const ContinuumState* cs) const
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

            //integrals.GetTwoElectronIntegral(k, e1, e2, e3, e4)
            double radial = 0;

            const State* s_1 = cs;
            const State* s_2 = excited->GetState(e2);
            const State* s_3 = excited->GetState(e3);
            const State* s_4 = excited->GetState(e4);

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
                density[p] = s_2->f[p] * s_4->f[p] + Constant::AlphaSquared * s_2->g[p] * s_4->g[p];
            }
            density.resize(excited->GetCore()->GetHFPotential().size());

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

            total += coeff * radial;
        }

        k = k+2;
    }

    return total;
}

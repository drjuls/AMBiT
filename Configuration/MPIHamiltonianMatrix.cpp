#ifdef _MPI
#include <mpi.h>
#include "Include.h"
#include "MPIHamiltonianMatrix.h"
#include "MPIMatrix.h"
#include "Universal/ScalapackMatrix.h"
#include "Universal/Eigensolver.h"
#include "Universal/Constant.h"

void MPIHamiltonianMatrix::GenerateMatrix()
{
    unsigned int i, j=0;

    // proc runs from (-NumProcessors) to (NumProcessors-1) cyclically.
    // Our part of the matrix corresponds to proc == rank1 or rank2
    int proc = - int(NumProcessors);
    int rank1 = int(ProcessorRank);
    int rank2 = -1 - int(ProcessorRank);

    if(M == NULL)
        M = new MPIMatrix(N, *configs);
    else
        M->Clear();

    M->WriteMode(true);

    // Loop through relativistic configurations
    RelativisticConfigList::const_iterator list_it = configs->begin();
    i = 0;
    while(list_it != configs->end())
    {
        if((proc == rank1) || (proc == rank2))
        {
            const ProjectionSet& proj_i = list_it->GetProjections();
            unsigned int proj_i_size = proj_i.size();
            unsigned int num_states_i = list_it->NumJStates();
            const double* coefficients_i = list_it->GetJCoefficients();

            RelativisticConfigList::const_iterator list_jt = list_it;
            j = i;
            
            while(list_jt != configs->end())
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
        }

        i += list_it->NumJStates();
        list_it++;

        proc++;
        if(proc == int(NumProcessors))
            proc = -proc;
    }
    *logstream << "Matrix Generated" << std::endl;
}

void MPIHamiltonianMatrix::WriteToFile(const std::string& filename)
{
    if(dynamic_cast<MPIMatrix*>(M))
        dynamic_cast<MPIMatrix*>(M)->WriteToFile(filename);
}

void MPIHamiltonianMatrix::PollMatrix()
{
    M->WriteMode(false);

    unsigned int i, j, count;
    double value;
    unsigned int proc_range[10];
    unsigned int range[10];
    for(i=0; i<10; i++)
        proc_range[i] = 0;

    // proc runs from (-NumProcessors) to (NumProcessors-1) cyclically.
    // Our part of the matrix corresponds to proc == rank1 or rank2
    int proc = - int(NumProcessors);
    int rank1 = int(ProcessorRank);
    int rank2 = -1 - int(ProcessorRank);

    // Loop through relativistic configurations
    RelativisticConfigList::const_iterator list_it = configs->begin();
    i = 0;
    while(list_it != configs->end())
    {
        if((proc == rank1) || (proc == rank2))
        {
            for(unsigned int jstate_i = 0; jstate_i < list_it->NumJStates(); jstate_i++)
            {
                for(j=i; j<N; j++)
                {
                    value = fabs(M->At(i, j));
                    if(value >= 1.)
                    {   proc_range[9]++;
                        if(i != j)
                            proc_range[9]++;
                    }
                    else
                    {   count = 0;
                        while(value > 1.e-16)
                        {   value = value/100.;
                            count++;
                        }
                        proc_range[count]++;
                        if(i != j)
                            proc_range[count]++;
                    }
                }
                i++;
            }
        }
        else
            i += list_it->NumJStates();

        list_it++;
        proc++;
        if(proc == int(NumProcessors))
            proc = -proc;
    }

    MPI::Intracomm& comm_world = MPI::COMM_WORLD;
    comm_world.Reduce(&proc_range, &range, 10, MPI::UNSIGNED, MPI::SUM, 0);

    if(ProcessorRank == 0)
    {   for(i=0; i<10; i++)
            *outstream << i << " " << range[i] << " " << double(range[i])/double(N*N)*100. << std::endl;
    }
}

void MPIHamiltonianMatrix::SolveMatrix(unsigned int num_solutions, Eigenstates& eigenstates, bool gFactors)
{
    if(N == 0)
    {   *outstream << "\nNo solutions" << std::endl;
        return;
    }

    M->WriteMode(false);
    ConfigFileGenerator* config_file_gen = dynamic_cast<ConfigFileGenerator*>(confgen);
    const std::set<Configuration>* leading_configs = confgen->GetLeadingConfigs();

    *outstream << "\nFinding solutions" << std::endl;

    unsigned int NumSolutions = mmin(num_solutions, N);
    
    double* V = new double[NumSolutions * N];
    double* E = new double[NumSolutions];

    Eigensolver solver;
    solver.MPISolveLargeSymmetric(M, E, V, N, NumSolutions);

    eigenstates.SetEigenvalues(E, NumSolutions);
    eigenstates.SetEigenvectors(V, NumSolutions);

    // Calculate g-Factors
    double* g_factors = NULL;
    if(gFactors && eigenstates.GetTwoJ())
    {   g_factors = new double[NumSolutions];
        GetgFactors(eigenstates, g_factors);
        eigenstates.SetgFactors(g_factors);
    }

    if(ProcessorRank == 0)
    {
        unsigned int i, j;

        *outstream << "Solutions for J = " << double(eigenstates.GetTwoJ())/2. << ", P = ";
        if(eigenstates.GetParity() == even)
            *outstream << "even:" << std::endl;
        else
            *outstream << "odd:" << std::endl;

        for(i=0; i<NumSolutions; i++)
        {
            unsigned int solution = i;

            *outstream << i << ": " << std::setprecision(8) << E[solution] << "    "
                << std::setprecision(12) << E[solution]*Constant::HartreeEnergy_cm << " /cm" << std::endl;

            // Get non-rel configuration percentages
            RelativisticConfigList::const_iterator list_it = configs->begin();
            std::map<Configuration, double> percentages;  // Map non-rel configurations to percentages

            j = 0;
            while(list_it != configs->end())
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

            // Find most important configuration, and print all leading configurations.
            std::map<Configuration, double>::const_iterator it_largest_percentage = percentages.begin();
            double largest_percentage = 0.0;

            std::map<Configuration, double>::const_iterator it = percentages.begin();
            while(it != percentages.end())
            {
                if(it->second > largest_percentage)
                {   it_largest_percentage = it;
                    largest_percentage = it->second;
                }

                if(it->second > 1.)
                    *outstream << std::setw(20) << it->first.Name() << "  "<< std::setprecision(2)
                        << it->second << "%" << std::endl;
                it++;
            }

            // If the most important configuration is a leading configuration, add this state to the config file.
            if(config_file_gen && (leading_configs->find(it_largest_percentage->first) != leading_configs->end()))
                config_file_gen->AddPercentages(percentages);

            if(g_factors)
                *outstream << "    g-factor = " << std::setprecision(5) << g_factors[solution] << std::endl;

            *outstream << std::endl;
        }
    }
}

void MPIHamiltonianMatrix::GetEigenvalues(const Eigenstates& eigenstates) const
{
    unsigned int NumSolutions = eigenstates.GetNumEigenvectors();
    const double* V = eigenstates.GetEigenvectors();

    double* my_total = new double[NumSolutions];
    double* coeff = new double[NumSolutions];
    unsigned int solution;
    for(solution = 0; solution < NumSolutions; solution++)
        my_total[solution] = 0.;

    // Iterate over different relativistic configurations
    unsigned int i = 0, j;

    int proc = - int(NumProcessors);
    int rank1 = int(ProcessorRank);
    int rank2 = -1 - int(ProcessorRank);

    RelativisticConfigList::const_iterator list_it = configs->begin();
    while(list_it != configs->end())
    {
        if((proc == rank1) || (proc == rank2))
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
                                my_total[solution] += coeff[solution] * matrix_element;
                        }

                        pj_it++; pj++;
                    }

                    pi_it++; pi++;
                }

                list_jt++; j+=num_states_j;
            }
        }

        i += list_it->NumJStates();
        list_it++;

        proc++;
        if(proc == int(NumProcessors))
            proc = -proc;
    }

    double* total;
    total = new double[NumSolutions];

    MPI::Intracomm& comm_world = MPI::COMM_WORLD;
    comm_world.Reduce(my_total, total, NumSolutions, MPI::DOUBLE, MPI::SUM, 0);

    if(ProcessorRank == 0)
    {   for(solution = 0; solution < NumSolutions; solution++)
            *outstream << solution << ": " 
                << std::setprecision(8) << total[solution] << "    "
                << total[solution]*Constant::HartreeEnergy_cm
                << " /cm" << std::endl;
    }

    delete[] total;
    delete[] my_total;
    delete[] coeff;
}

void MPIHamiltonianMatrix::GetgFactors(const Eigenstates& eigenstates, double* g_factors) const
{
    GetgFactors(eigenstates.GetNumEigenvectors(), eigenstates.GetEigenvectors(), eigenstates.GetTwoJ(), g_factors);
}

void MPIHamiltonianMatrix::GetgFactors(unsigned int NumSolutions, const double* V, unsigned int two_j, double* g_factors) const
{
    // Case where J=0
    if(two_j == 0)
        return;

    double* total = new double[NumSolutions];
    double* coeff = new double[NumSolutions];
    unsigned int solution;

    for(solution = 0; solution < NumSolutions; solution++)
        total[solution] = 0.;

    unsigned int diff[4];   // Storage for projection differences.
    unsigned int num_electrons = configs->front().NumParticles();

    // Iterate over different relativistic configurations
    unsigned int i = 0, j;

    int proc = - int(NumProcessors);
    int rank1 = int(ProcessorRank);
    int rank2 = -1 - int(ProcessorRank);

    RelativisticConfigList::const_iterator list_it = configs->begin();
    while(list_it != configs->end())
    {
        if((proc == rank1) || (proc == rank2))
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
        }

        i += list_it->NumJStates();
        list_it++;

        proc++;
        if(proc == int(NumProcessors))
            proc = -proc;
    }

    MPI::Intracomm& comm_world = MPI::COMM_WORLD;
    comm_world.Allreduce(total, g_factors, NumSolutions, MPI::DOUBLE, MPI::SUM);

    double J = double(two_j)/2.;
    for(solution = 0; solution < NumSolutions; solution++)
        g_factors[solution] = g_factors[solution]/J + 1.;

    delete[] coeff;
    delete[] total;
}

#ifdef _SCALAPACK
void MPIHamiltonianMatrix::SolveScalapack(const std::string& filename, double eigenvalue_limit, Eigenstates& eigenstates, bool gFactors)
{
    if(N == 0)
    {   *outstream << "\nNo solutions" << std::endl;
        return;
    }

    if(M)
        delete M;

    M = new ScalapackMatrix(N);
    ScalapackMatrix* SM = dynamic_cast<ScalapackMatrix*>(M);
    
    SM->ReadTriangle(filename);
    
    ConfigFileGenerator* config_file_gen = dynamic_cast<ConfigFileGenerator*>(confgen);
    const std::set<Configuration>* leading_configs = confgen->GetLeadingConfigs();

    *outstream << "\nFinding solutions" << std::endl;

    unsigned int NumSolutions = mmin(50, N);

    double* E = new double[N];
    double* V = new double[NumSolutions * N];

    bool multiple_chunks = false;

    SM->Diagonalise(E);

    eigenstates.SetEigenvalues(E, N);

    if(ProcessorRank == 0)
    {
        *outstream << "Solutions for J = " << double(eigenstates.GetTwoJ())/2. << ", P = ";
        if(eigenstates.GetParity() == even)
            *outstream << "even:" << std::endl;
        else
            *outstream << "odd:" << std::endl;
    }

    // Calculate g-Factors
    double* g_factors = NULL;
    if(gFactors && eigenstates.GetTwoJ())
        g_factors = new double[NumSolutions];

    unsigned int i, j, count;

    i = 0;
    while((i < N) && (E[i] < eigenvalue_limit))
    {
        // Second pass: keep the first chunk, and continue with the others
        if((i > 0) && !multiple_chunks)
        {
            multiple_chunks = true;
            V = new double[NumSolutions * N];
        }

        // Do a chunk of eigenvectors, maximum size NumSolutions
        for(count = 0; (count < NumSolutions) && (i+count < N) &&
                                    (E[i+count] < eigenvalue_limit); count++)
        {
            SM->GetColumn(i+count, &V[count*N]);
        }

        // Change NumSolutions for GetgFactors(), in case we aren't doing a maximum-sized chunk
        NumSolutions = count;

        // Save first chunk
        if(i == 0)
            eigenstates.SetEigenvectors(V, NumSolutions);

        if(g_factors)
            GetgFactors(NumSolutions, V, eigenstates.GetTwoJ(), g_factors);

        for(count = 0; count < NumSolutions; count++)
        {
            if(ProcessorRank == 0)
            {
                *outstream << i << ": " << std::setprecision(8) << E[i] << "    "
                    << std::setprecision(12) << E[i]*Constant::HartreeEnergy_cm << " /cm" << std::endl;

                // Get non-rel configuration percentages
                RelativisticConfigList::const_iterator list_it = configs->begin();
                std::map<Configuration, double> percentages;  // Map non-rel configurations to percentages

                j = 0;
                while(list_it != configs->end())
                {
                    Configuration nrconfig(list_it->GetNonRelConfiguration());
                    if(percentages.find(nrconfig) == percentages.end())
                        percentages[nrconfig] = 0.;

                    for(unsigned int Jstate = 0; Jstate < list_it->NumJStates(); Jstate++)
                    {
                        double coeff = V[count*N + j];
                        coeff = coeff * coeff * 100;

                        percentages[nrconfig] += coeff;
                        j++;
                    }

                    list_it++;
                }

                // Find most important configuration, and print all leading configurations.
                std::map<Configuration, double>::const_iterator it_largest_percentage = percentages.begin();
                double largest_percentage = 0.0;

                std::map<Configuration, double>::const_iterator it = percentages.begin();
                while(it != percentages.end())
                {
                    if(it->second > largest_percentage)
                    {   it_largest_percentage = it;
                        largest_percentage = it->second;
                    }

                    if(it->second > 5.)
                        *outstream << std::setw(20) << it->first.Name() << "  "<< std::setprecision(2)
                            << it->second << "%" << std::endl;
                    it++;
                }

                // If the most important configuration is a leading configuration, add this state to the config file.
                if(config_file_gen && (leading_configs->find(it_largest_percentage->first) != leading_configs->end()))
                    config_file_gen->AddPercentages(percentages);

                if(g_factors)
                    *outstream << "    g-factor = " << std::setprecision(5) << g_factors[count] << std::endl;

                *outstream << std::endl;
            }

            i++;
        }
    }

    if(g_factors)
        delete[] g_factors;

    if(multiple_chunks)
        delete[] V;
}
#endif

#endif


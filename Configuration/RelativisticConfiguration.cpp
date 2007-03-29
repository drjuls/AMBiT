#include "Include.h"
#include "RelativisticConfiguration.h"
#include "ConfigGenerator.h"
#include "Universal/Eigensolver.h"

RelativisticConfiguration::RelativisticConfiguration(const RelativisticConfiguration& other):
    Configuration(other), projections(other.projections),
    num_states(other.num_states)
{
    if(num_states)
    {   j_coefficients = new double[num_states*projections.size()];
        for(unsigned int i=0; i<num_states*projections.size(); i++)
            j_coefficients[i] = other.j_coefficients[i];
    }
}

StateInfo RelativisticConfiguration::GetInfo() const
{
    return StateInfo(Configuration::GetInfo());
}

bool RelativisticConfiguration::AddSingleParticle(const StateInfo& info)
{
    std::map<StateInfo, unsigned int>::iterator a_it = Config.find(info);
    if(a_it != Config.end())
    {
        StateInfo rel_info(a_it->first);
        if(a_it->second >= (unsigned int)(2*abs(rel_info.Kappa()))) // maximum number of electrons
        {   a_it->second = 2*abs(rel_info.Kappa());
            return false;
        }
        else
            a_it->second = a_it->second + 1;
    }
    else
    {
        Config[info] = 1;
    }
    return true;
}

bool RelativisticConfiguration::GenerateProjections(int two_m)
{
    if(num_states)
    {   delete[] j_coefficients;
        j_coefficients = NULL;
        num_states = 0;
    }
    projections.clear();

    // Get electron vector
    std::vector<ElectronInfo> electrons;
    First();
    while(!AtEnd())
    {
        StateInfo rinfo = GetInfo();
        for(unsigned int i=0; i<GetOccupancy(); i++)
        {   electrons.push_back(ElectronInfo(rinfo.PQN(), rinfo.Kappa(), 0));
        }

        Next();
    }

    // Get projections by recursion
    DoElectron(electrons, 0);
    projections.sort();
    projections.unique();

    // Eliminate projections with the wrong value of M
    ProjectionSet::iterator p = projections.begin();

    while(p != projections.end())
    {   if(p->GetTwoM() != two_m)
            p = projections.erase(p);
        else
            p++;
    }

    return (projections.size() != 0);
}

void RelativisticConfiguration::DoElectron(std::vector<ElectronInfo>& electrons, unsigned int index)
{
    if(index >= electrons.size())
    {
        Projection proj;
        for(unsigned int i=0; i<electrons.size(); i++)
        {   proj.Add(electrons[i]);
        }
        proj.Sort();
        projections.push_back(proj);
    }
    else
    {
        const ElectronInfo& e_info = electrons[index];

        int two_m = 0;
        if(index != 0)
        {
            const ElectronInfo& prev_info = electrons[index-1];
            if(StateInfo(e_info) == StateInfo(prev_info))
            {   
                two_m = prev_info.TwoM() + 2;
                if(two_m > int(e_info.TwoJ()))
                    return;
            }
            else
            {   two_m = - int(e_info.TwoJ());
            }
        }
        else
            two_m = - int(e_info.TwoJ());

        while(two_m <= int(e_info.TwoJ()))
        {
            electrons[index].SetTwoM(two_m);
            DoElectron(electrons, index+1);
            two_m += 2;
        }
    }
}

int RelativisticConfiguration::GetTwiceMaxProjection() const
{
    std::map<StateInfo, unsigned int>::const_iterator m_it = Config.begin();
    int proj = 0;
    while(m_it != Config.end())
    {
        int j = m_it->first.TwoJ();
        for(unsigned int i=0; i<m_it->second; i++)
        {   proj += j;
            j = j - 2;
        }

        m_it++;
    }

    return proj;
}

std::string RelativisticConfiguration::Name() const
{
    std::map<StateInfo, unsigned int>::const_iterator m_it = Config.begin();
    std::string name;
    while(m_it != Config.end())
    {
        name.append(" " + StateInfo(m_it->first).Name());

        char buffer[20];
        sprintf(buffer, "%d", m_it->second);
        std::string occupancy(buffer);
        name.append(buffer);

        m_it++;
    }
    return name;
}

void RelativisticConfiguration::SetJCoefficients(unsigned int num_Jstates, double* coefficients)
{
    if(num_states)
    {   delete[] j_coefficients;
    }

    if(num_Jstates)
        j_coefficients = coefficients;
    else
        j_coefficients = NULL;

    num_states = num_Jstates;
}

bool RelativisticConfiguration::GenerateJCoefficients(double J)
{
    unsigned int N = projections.size();

    if(N == 0)
        return false;

    if(num_states)
    {   delete[] j_coefficients;
        j_coefficients = NULL;
        num_states = 0;
    }

    // Generate the matrix
    unsigned int i, j;
    double* M = new double[N*N];

    ProjectionSet::const_iterator i_it = projections.begin();
    ProjectionSet::const_iterator j_it;
    i = 0;
    while(i_it != projections.end())
    {
        j_it = i_it;
        j = i;
        while(j_it != projections.end())
        {
            double matrix_element = GetJSquared(*i_it, *j_it);
            M[i*N + j] = M[j*N + i] = matrix_element;
            j_it++; j++;
        }

        i_it++; i++;
    }

    // Solve the matrix
    double* V = new double[N];  // eigenvalues
    for(i=0; i<N; i++)
        V[i] = 0.;
    
    Eigensolver E;
    E.SolveSmallSymmetric(M, V, N);

    // Count number of good eigenvalues
    double JSquared = J * (J+1.);
    num_states = 0;
    for(i=0; i<N; i++)
        if(fabs(V[i] - JSquared) < 1.e-6)
            num_states++;

    // Transfer eigenvalues
    if(num_states)
    {
        j_coefficients = new double[num_states * N];
        unsigned int count = 0;
        for(i=0; i<N; i++)
        {
            if(fabs(V[i] - JSquared) < 1.e-6)
            {
                for(j = 0; j < N; j++)
                    j_coefficients[count*N + j] = M[i*N + j];

                count++;
            }
        }
    }

    delete[] M;
    delete[] V;

    return (num_states != 0);
}

double RelativisticConfiguration::GetJSquared(const Projection& first, const Projection& second) const
{
    double ret = 0.;
    unsigned int i, j;
    unsigned int diff[4];

    int numdiff = Projection::GetProjectionDifferences(first, second, diff);
    if(numdiff == 0)
    {   // <J^2> += Sum_i (j_i)(j_i + 1)  + Sum_(i<j) 2(m_i)(m_j)
        i=0;
        while(i < first.Size())
        {
            ret += first[i].J() * (first[i].J() + 1);

            j = i+1;
            while(j < second.Size())
            {   ret += 2. * first[i].M() * second[j].M();
                j++;
            }
            
            i++;
        }
        // <J^2> += Sum_(i<j) (j+_i)(j-_j) + (j-_i)(j+_j)
        for(i=0; i<second.Size(); i++)
        {   for(j=i+1; j<second.Size(); j++)
            {
                if((first[i].Kappa() == first[j].Kappa()) && (first[i].PQN() == first[j].PQN())
                   && (abs(first[i].TwoM() - first[j].TwoM()) == 2))
                {
                    double J = first[i].J();
                    double M;
                    if(first[i].TwoM() > first[j].TwoM())
                        M = first[j].M();
                    else
                        M = first[i].M();

                    double value = (J * (J + 1) - M * (M + 1));
                    if((j-i)%2)
                        ret -= value;
                    else
                        ret += value;
                }
                else
                    break;  // WARNING: Assumes ordering of electrons in projection
            }               //          is grouped by PQN and Kappa.
        }
    }
    else if(abs(numdiff) == 2)
    {
        const ElectronInfo& f1 = first[diff[0]];
        const ElectronInfo& s1 = second[diff[1]];
        const ElectronInfo& f2 = first[diff[2]];
        const ElectronInfo& s2 = second[diff[3]];

        if((f1.Kappa() == s1.Kappa()) && (f1.PQN() == s1.PQN())
           && (f2.Kappa() == s2.Kappa()) && (f2.PQN() == s2.PQN())
           && (abs(f1.TwoM() - s1.TwoM()) == 2) && (abs(f2.TwoM() - s2.TwoM()) == 2)
           && (f1.TwoM() + f2.TwoM() == s1.TwoM() + s2.TwoM()))
        {
            ret = sqrt(f1.J() * (f1.J() + 1.) - f1.M() * s1.M())*
                  sqrt(f2.J() * (f2.J() + 1.) - f2.M() * s2.M());
        }
        else if((f1.Kappa() == s2.Kappa()) && (f1.PQN() == s2.PQN())
           && (f2.Kappa() == s1.Kappa()) && (f2.PQN() == s1.PQN())
           && (abs(f1.TwoM() - s2.TwoM()) == 2) && (abs(f2.TwoM() - s1.TwoM()) == 2)
           && (f1.TwoM() + f2.TwoM() == s1.TwoM() + s2.TwoM()))
        {
            ret = sqrt(f1.J() * (f1.J() + 1.) - f1.M() * s2.M())*
                  sqrt(f2.J() * (f2.J() + 1.) - f2.M() * s1.M());
        }

        if(numdiff < 0)
            ret = -ret;
    }

    return ret;
}

Configuration RelativisticConfiguration::GetNonRelConfiguration() const
{
    Configuration ret;

    std::map<StateInfo, unsigned int>::const_iterator c_it = Config.begin();
    while(c_it != Config.end())
    {
        ret.SetIterator(NonRelInfo(c_it->first));
        if(ret.AtEnd())
        {   ret.AddSingleParticle(NonRelInfo(c_it->first));
            ret.SetOccupancy(NonRelInfo(c_it->first), c_it->second);
        }
        else
        {   unsigned int occ = ret.GetOccupancy();
            ret.SetOccupancy(occ + c_it->second);
        }

        c_it++;
    }

    return ret;
}

void RelativisticConfiguration::Write(FILE* fp) const
{
    // Write config
    unsigned int size = Config.size();
    fwrite(&size, sizeof(unsigned int), 1, fp);    

    std::map<StateInfo, unsigned int>::const_iterator cit;
    cit = Config.begin();

    while(cit != Config.end())
    {
        unsigned int pqn = cit->first.PQN();
        int kappa = cit->first.Kappa();
        unsigned int occupancy = cit->second;
        
        // Write PQN, Kappa, occupancy
        fwrite(&pqn, sizeof(unsigned int), 1, fp);    
        fwrite(&kappa, sizeof(int), 1, fp);    
        fwrite(&occupancy, sizeof(unsigned int), 1, fp);    

        cit++;
    }

    // Write number of projections. Projections will be generated again upon read,
    // so this serves as a sanity check.
    size = projections.size();
    fwrite(&size, sizeof(unsigned int), 1, fp);

    if(size)
    {   // Write twoM (all projections should have the same value).
        int twoM = projections.front().GetTwoM();
        fwrite(&twoM, sizeof(int), 1, fp);

        // Write JCoefficients. (These can be costly to generate in some cases.)
        fwrite(&num_states, sizeof(unsigned int), 1, fp);
        if(num_states)
            fwrite(j_coefficients, sizeof(double), num_states * projections.size(), fp);
    }
}

void RelativisticConfiguration::Read(FILE* fp)
{
    // Clear the current configuration
    Config.clear();
    if(num_states)
        delete[] j_coefficients;
    
    // Read config
    unsigned int size;
    fread(&size, sizeof(unsigned int), 1, fp);
    
    for(unsigned int i = 0; i < size; i++)
    {
        unsigned int pqn;
        int kappa;
        unsigned int occupancy;

        // Read PQN, Kappa, occupancy
        fread(&pqn, sizeof(unsigned int), 1, fp);    
        fread(&kappa, sizeof(int), 1, fp);    
        fread(&occupancy, sizeof(unsigned int), 1, fp);    

        SetOccupancy(StateInfo(pqn, kappa), occupancy);
    }

    // Generate projections and check that the same number were stored.
    fread(&size, sizeof(unsigned int), 1, fp);
    
    if(size)
    {   // Read twoM.
        int twoM;
        fread(&twoM, sizeof(int), 1, fp);

        GenerateProjections(twoM);
        if(projections.size() != size)
        {   *errstream << "RelativisticConfiguration::Read(): " << Name()
                       << " has incorrect number of projections stored." << std::endl;
            exit(1);
        }

        // Read JCoefficients.
        fread(&num_states, sizeof(unsigned int), 1, fp);
        if(num_states)
        {   j_coefficients = new double[num_states * projections.size()];
            fread(j_coefficients, sizeof(double), num_states * projections.size(), fp);
        }
    }
}

#include "Include.h"
#include "RelativisticConfiguration.h"
#include "ConfigGenerator.h"
#include "Universal/Eigensolver.h"

RelativisticInfo RelativisticConfiguration::GetInfo() const
{
    return RelativisticInfo(Configuration::GetInfo());
}

bool RelativisticConfiguration::AddSingleParticle(const SingleParticleInfo& info)
{
    std::map<SingleParticleInfo, unsigned int>::iterator a_it = Config.find(info);
    if(a_it != Config.end())
    {
        RelativisticInfo rel_info(a_it->first);
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
    j_coefficients.clear();
    projections.clear();

    // Get electron vector
    std::vector<ElectronInfo> electrons;
    First();
    while(!AtEnd())
    {
        RelativisticInfo rinfo = GetInfo();
        for(unsigned int i=0; i<GetOccupancy(); i++)
        {   electrons.push_back(ElectronInfo(rinfo.PQN(), rinfo.Kappa(), 0));
        }

        Next();
    }

    // Get projections by recursion
    DoElectron(electrons, 0);

    // Eliminate projections with the wrong value of M
    ProjectionSet::iterator p = projections.begin();
    ProjectionSet::iterator nextp = p;

    while(p != projections.end())
    {
        nextp++;
        if(p->GetTwoM() != two_m)
            projections.erase(p);

        p = nextp;
    }

    if(GetProjectionCoefficients(double(two_m)/2.))
    {
        //std::cout << " ---- " << std::endl;
        //p = projections.begin();
        //unsigned int i, j = 0;
        //while(p != projections.end())
        //{   std::cout << p->Name();
        //    for(i = 0; i < j_coefficients.size(); i++)
        //        std::cout << "\t" << j_coefficients[i][j];
        //    std::cout << std::endl;
        //    p++; j++;
        //}
        return true;
    }
    else
        return false;
}

void RelativisticConfiguration::DoElectron(std::vector<ElectronInfo> electrons, unsigned int index)
{
    if(index >= electrons.size())
    {
        Projection proj;
        for(unsigned int i=0; i<electrons.size(); i++)
        {   proj.Add(electrons[i]);
        }
        proj.Sort();
        projections.insert(proj);
    }
    else
    {
        ElectronInfo e_info = electrons[index];

        int two_m = 0;
        if(index != 0)
        {
            ElectronInfo prev_info = electrons[index-1];
            if(RelativisticInfo(e_info) == RelativisticInfo(prev_info))
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
    std::map<SingleParticleInfo, unsigned int>::const_iterator m_it = Config.begin();
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
    std::map<SingleParticleInfo, unsigned int>::const_iterator m_it = Config.begin();
    std::string name;
    while(m_it != Config.end())
    {
        name.append(" " + RelativisticInfo(m_it->first).Name());

        char buffer[20];
        sprintf(buffer, "%d", m_it->second);
        std::string occupancy(buffer);
        name.append(buffer);

        m_it++;
    }
    return name;
}

bool RelativisticConfiguration::GetProjectionCoefficients(double J)
{
    unsigned int N = projections.size();

    if(N == 0)
        return false;

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

    // Transfer eigenvalues
    double JSquared = J * (J+1.);
    for(i=0; i<N; i++)
    {
        if(fabs(V[i] - JSquared) < 1.e-6)
        {
            std::vector<double> coeff(N);
            for(j = 0; j < N; j++)
                coeff[j] = M[i*N + j];

            j_coefficients.push_back(coeff);
        }
    }

    delete[] M;
    delete[] V;

    return (j_coefficients.size() != 0);
}

double RelativisticConfiguration::GetJSquared(const Projection& first, const Projection& second) const
{
    double ret = 0.;
    unsigned int i, j;

    if(first == second)
    {
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
    }
    
    for(i=0; i<second.Size(); i++)
        for(j=0; j<second.Size(); j++)
        {
            if(i != j)
            {
                Projection proj(second);
                ElectronInfo e1(proj[i]);
                ElectronInfo e2(proj[j]);

                proj[i].SetTwoM(proj[i].TwoM()+2);
                proj[j].SetTwoM(proj[j].TwoM()-2);

                if((proj[i].TwoM() <= int(proj[i].TwoJ())) && (proj[j].TwoM() >= -int(proj[j].TwoJ())))
                {
                    bool sort = proj.Sort();
                    if(first == proj)
                    {
                        double coeff = sqrt(e1.J() * (e1.J() + 1.) - e1.M() * (e1.M() + 1.))*
                                       sqrt(e2.J() * (e2.J() + 1.) - e2.M() * (e2.M() - 1.));
                        if(sort)
                            ret = ret - coeff;
                        else
                            ret = ret + coeff;
                    }
                }
            }
        }

    return ret;
}

Configuration RelativisticConfiguration::GetNonRelConfiguration() const
{
    Configuration ret;

    std::map<SingleParticleInfo, unsigned int>::const_iterator c_it = Config.begin();
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

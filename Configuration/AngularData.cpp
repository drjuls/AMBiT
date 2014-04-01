#include "AngularData.h"
#include "Include.h"
#include "RelativisticConfiguration.h"
#include "Universal/Eigensolver.h"
#include <numeric>

AngularData::AngularData(const RelativisticConfiguration& config, int two_m):
    two_m(two_m), two_j(-1), num_CSFs(0), CSFs(nullptr)
{
    GenerateProjections(config, two_m);
}

AngularData::AngularData(const RelativisticConfiguration& config, int two_m, int two_j):
    two_m(two_m), two_j(-1), num_CSFs(0), CSFs(nullptr)
{
    GenerateProjections(config, two_m);
    GenerateCSFs(config, two_j);
}

AngularData::~AngularData()
{
    if(CSFs)
        delete[] CSFs;
}

bool AngularData::ProjectionCompare(const std::vector<int>& first, const std::vector<int>& second)
{
    auto it_first = first.begin();
    auto it_second = second.begin();
    while(it_first != first.end() && it_second != second.end())
    {
        if(*it_first < *it_second)
            return true;
        else if(*it_first > *it_second)
            return false;

        it_first++;
        it_second++;
    }

    if(first.size() < second.size())
        return true;
    else
        return false;
}

int AngularData::GenerateProjections(const RelativisticConfiguration& config, int two_m)
{
    // Start with maximum M.
    // At each step populate all possible projections with M-1.
    // Stop when desired M is reached.

    // Create set of boxes, one for each orbital, each with a list of projections to be filled.
    int projection_size = config.ExcitationNumber();
    std::vector<int> box_number(projection_size);
    std::vector<int> min_projection(projection_size);
    std::vector<int> max_projection(projection_size);

    int box_count = config.size();
    std::vector< std::pair<int, int> > box_particle_bounds(box_count);

    int count = 0;
    int box = 0;
    for(auto config_it = config.begin(); config_it != config.end(); config_it++)
    {
        int current_two_m = config_it->first.TwoJ();
        int num_particles = abs(config_it->second);
        box_particle_bounds[box].first = count;
        box_particle_bounds[box].second = count + num_particles;

        for(int i = 0; i < num_particles; i++)
        {
            max_projection[count + i] = current_two_m;
            min_projection[count + num_particles-1 - i] = -current_two_m;
            box_number[count + i] = box;
            current_two_m -= 2;
        }

        count += num_particles;
        box++;
    }

    max_projection.shrink_to_fit();
    int current_two_m = std::accumulate(max_projection.begin(), max_projection.end(), 0);

    std::vector< std::list< std::vector<int> > > boxes(box_count);
    if(current_two_m >= two_m)
        boxes[0].push_back(max_projection);

    // Iterative method to find all projections with correct m.
    while(current_two_m > two_m)
    {
        std::vector< std::list< std::vector<int> > > next_boxes(box_count);

        // For each projection in boxes[i], apply M_j-- for i <= j and
        // store the result in next_boxes[j].
        for(int box_i = 0; box_i < box_count; box_i++)
        {
            for(std::vector<int>& existing_proj: boxes[box_i])
            {
                int i = box_particle_bounds[box_i].first;

                for(int j = i; j < projection_size; j++)
                {
                    std::vector<int> new_proj(existing_proj);
                    if(j == projection_size-1 ||            // Last particle or
                       box_number[j] != box_number[j+1] ||  // next particle is in a different orbital or
                       new_proj[j] - new_proj[j+1] > 2)     // there is room to apply J^-
                    {
                        new_proj[j] -= 2;
                        if(new_proj[j] >= min_projection[j])
                        {
                            // new_proj is legit, store it
                            next_boxes[box_number[j]].push_back(new_proj);
                        }
                    }
                }
            }

            next_boxes[box_i].sort(ProjectionCompare);
            next_boxes[box_i].unique();
        }

        boxes.swap(next_boxes);
        current_two_m -= 2;
    }

    // Sort and merge projections lists
    projections.clear();
    for(auto& list: boxes)
    {   list.sort(ProjectionCompare);
        projections.merge(list, ProjectionCompare);
    }
    projections.unique();

    return projections.size();
}

int AngularData::GenerateCSFs(const RelativisticConfiguration& config, int two_j)
{
    unsigned int N = projections.size();
    this->two_j = two_j;

    // Clear existing
    if(num_CSFs)
    {   delete[] CSFs;
        CSFs = nullptr;
        num_CSFs = 0;
    }

    if(N == 0 || two_j < abs(two_m))
        return 0;

    // Generate the matrix
    unsigned int i, j;
    double* M = new double[N*N];

    // Make list of Projections.
    std::list< Projection > real_Projection_list;
    for(auto& p: projections)
        real_Projection_list.push_back(Projection(config, p));

    auto i_it = real_Projection_list.begin();
    auto j_it = i_it;

    i = 0;
    while(i_it != real_Projection_list.end())
    {
        j_it = i_it;
        j = i;
        while(j_it != real_Projection_list.end())
        {
            double matrix_element = GetJSquared(*i_it, *j_it);
            M[i*N + j] = M[j*N + i] = matrix_element;
            j_it++; j++;
        }

        i_it++; i++;
    }

    // Solve the matrix
    double* V = new double[N];  // eigenvalues
    memset(V, 0, N * sizeof(double));

    Eigensolver E;
    E.SolveSmallSymmetric(M, V, N);

    // Count number of good eigenvalues
    double JSquared = double(two_j * (two_j + 2.)) / 4.;
    num_CSFs = 0;
    for(i=0; i<N; i++)
    {   // Check that all eigenvalues are good
        // j^2 + j - V = 0
        double TwoJ = (std::sqrt(1. + 4 * V[i]) - 1.);
        if(fabs(floor(TwoJ + 0.5) - TwoJ) > 1.e-6)
        {   *errstream << "AngularData::GenerateCSFs(): generated noninteger TwoJ:\n"
                       << "    config: " << config.Name()
                       << "    eigenvalue = " << V[i] << std::endl;
        }

        if(fabs(V[i] - JSquared) < 1.e-6)
            num_CSFs++;
    }
    *outstream << std::endl;

    // Transfer eigenvalues
    if(num_CSFs)
    {
        CSFs = new double[N * num_CSFs];
        unsigned int count = 0;
        for(i=0; i<N; i++)
        {
            if(fabs(V[i] - JSquared) < 1.e-6)
            {
                for(j = 0; j < N; j++)
                    CSFs[j * num_CSFs + count] = M[i*N + j];

                count++;
            }
        }
    }

    delete[] M;
    delete[] V;

    return num_CSFs;
}

double AngularData::GetJSquared(const Projection& first, const Projection& second) const
{
    double ret = 0.;
    unsigned int i, j;
    unsigned int diff[4];

    int numdiff = Projection::GetProjectionDifferences(first, second, diff);
    if(numdiff == 0)
    {   // <J^2> += Sum_i (j_i)(j_i + 1)  + Sum_(i<j) 2(m_i)(m_j)
        i=0;
        while(i < first.size())
        {
            ret += first[i].J() * (first[i].J() + 1);

            j = i+1;
            while(j < second.size())
            {   ret += 2. * first[i].M() * second[j].M();
                j++;
            }

            i++;
        }
        // <J^2> += Sum_(i<j) (j+_i)(j-_j) + (j-_i)(j+_j)
        for(i=0; i<second.size(); i++)
        {
            for(j=i+1; j<second.size(); j++)
            {
                // Assume ordering from large M to small M.
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
                    ret -= value;
//                    if((j-i)%2)
//                        ret -= value;
//                    else
//                        ret += value;
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
           && (f1.TwoM() + f2.TwoM() == s1.TwoM() + s2.TwoM())) // This should be automatic?
        {
            ret = sqrt(f1.J() * (f1.J() + 1.) - f1.M() * s1.M())*
                  sqrt(f2.J() * (f2.J() + 1.) - f2.M() * s2.M());
        }
        else if((f1.Kappa() == s2.Kappa()) && (f1.PQN() == s2.PQN())
                && (f2.Kappa() == s1.Kappa()) && (f2.PQN() == s1.PQN())
                && (abs(f1.TwoM() - s2.TwoM()) == 2) && (abs(f2.TwoM() - s1.TwoM()) == 2)
                && (f1.TwoM() + f2.TwoM() == s1.TwoM() + s2.TwoM()))
        {
            ret = - sqrt(f1.J() * (f1.J() + 1.) - f1.M() * s2.M())*
                    sqrt(f2.J() * (f2.J() + 1.) - f2.M() * s1.M());
        }

        if(numdiff < 0)
            ret = -ret;
    }
    
    return ret;
}

int AngularData::GenerateCSFs(const AngularData::ConfigKeyType& key, int two_j)
{
    this->two_j = two_j;

    // Create equivalent RelativisticConfiguration and use it.
    RelativisticConfiguration rconfig;
    int prev_kappa = 0;
    int pqn = 1;

    for(auto pair: key)
    {   // Set PQN to some integer
        if(pair.first == prev_kappa)
            pqn++;
        else
            pqn = 1;

        rconfig[OrbitalInfo(pqn, pair.first)] = pair.second;
    }

    return GenerateCSFs(rconfig, two_j);
}

AngularDataLibrary::AngularDataLibrary(int particle_number, int two_m, int two_j):
    two_m(two_m), two_j(two_j)
{
    filename = itoa(particle_number) + "." + itoa(two_m) + "." + itoa(two_j) + ".angular";
}

pAngularData AngularDataLibrary::operator[](const RelativisticConfiguration& config)
{
    KeyType mykey(GenerateKey(config));

    auto ret = library.find(mykey);
    if(ret != library.end())
        return ret->second;

    // Check whether this config can have correct M
    if(config.GetTwiceMaxProjection() < abs(two_m))
        return nullptr;

    pAngularData ang(new AngularData(config, two_m));
    library[mykey] = ang;
    return ang;
}

AngularDataLibrary::KeyType AngularDataLibrary::GenerateKey(const RelativisticConfiguration& config) const
{
    KeyType key;
    for(auto& config_it: config)
    {   key.push_back(std::make_pair(config_it.first.Kappa(), config_it.second));
    }
    return key;
}

void AngularDataLibrary::GenerateCSFs()
{
    for(auto& pair: library)
    {   auto& pAng = pair.second;
        pAng->GenerateCSFs(pair.first, two_j);
    }
}


#include "AngularData.h"
#include "Include.h"
#include "RelativisticConfiguration.h"
#include "Universal/Eigensolver.h"
#include "ManyBodyOperator.h"
#include <numeric>
#include <dirent.h>

AngularData::AngularData(int two_m):
    two_m(two_m), two_j(-1), num_CSFs(0), CSFs(nullptr), have_CSFs(false)
{}

AngularData::AngularData(const RelativisticConfiguration& config, int two_m):
    two_m(two_m), two_j(-1), num_CSFs(0), CSFs(nullptr), have_CSFs(false)
{
    GenerateProjections(config, two_m);
}

AngularData::AngularData(const RelativisticConfiguration& config, int two_m, int two_j):
    two_m(two_m), two_j(-1), num_CSFs(0), CSFs(nullptr), have_CSFs(false)
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

    int projection_size = config.ParticleNumber();

    // Skip if vacuum
    if(projection_size == 0)
    {
        projections.clear();
        projections.push_back(std::vector<int>());
        return 1;
    }

    // Create set of boxes, one for each orbital, each with a list of projections to be filled.
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
    if(CSFs)
    {   delete[] CSFs;
        CSFs = nullptr;
        have_CSFs = false;
    }

    if(N == 0 || two_j < abs(two_m))
        return 0;

    // Generate the matrix
    unsigned int i, j;
    double* M = new double[N*N];

    // Make list of Projections.
    std::list< Projection > real_Projection_list;
    for(const auto& p: projections)
        real_Projection_list.push_back(Projection(config, p));

    auto i_it = real_Projection_list.begin();
    auto j_it = i_it;

    ManyBodyOperator<const JSquaredOperator*, const JSquaredOperator*> J_squared(&J_squared_operator, &J_squared_operator);

    i = 0;
    while(i_it != real_Projection_list.end())
    {
        j_it = i_it;
        j = i;
        while(j_it != real_Projection_list.end())
        {
            double matrix_element = J_squared.GetMatrixElement(*i_it, *j_it);
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
        double TwoJ = (std::sqrt(1. + 4. * V[i]) - 1.);
        if(fabs(floor(TwoJ + 0.5) - TwoJ) > 1.e-6)
        {   *errstream << "AngularData::GenerateCSFs(): generated noninteger TwoJ:\n"
                       << "    config: " << config.Name()
                       << "    eigenvalue = " << V[i] << std::endl;
        }

        if(fabs(V[i] - JSquared) < 1.e-6)
            num_CSFs++;
    }

    have_CSFs = true;

    // Transfer eigenvalues
    if(num_CSFs)
    {
        CSFs = new double[N * num_CSFs];
        unsigned int count = 0;
        for(i=0; i < N && count < num_CSFs; i++)
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
        {   pqn = 1;
            prev_kappa = pair.first;
        }

        rconfig[OrbitalInfo(pqn, pair.first)] = pair.second;
    }

    return GenerateCSFs(rconfig, two_j);
}

AngularDataLibrary::AngularDataLibrary(int electron_number, const Symmetry& sym, int two_m, std::string lib_directory):
    electron_number(electron_number), two_m(two_m), two_j(sym.GetTwoJ())
{
    if(lib_directory == "")
        filename.clear();

    else
    {   // Check library directory exists
        if(!opendir(lib_directory.c_str()))
        {   *errstream << "AngularDataLibarary: unknown directory path " << lib_directory << std::endl;
            exit(1);
        }

        filename = lib_directory;
        if(filename[filename.length()-1] != '/')
            filename.append("/");

        std::string holelike = (electron_number<0)?"h":"";
        filename += itoa(abs(electron_number)) + holelike + "." + sym.GetString() + "." + itoa(two_m) + ".angular";
    }
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
        if(pAng->CSFs_calculated() == false)
            pAng->GenerateCSFs(pair.first, two_j);
    }
}

/** Structure of *.angular files:
    (Note: number of particles and symmetry is stored in filename)
    - (int) number of stored AngularData objects
    then for each AngularData object
        - (int) key size = number of pairs (kappa, num particles)
        - key
        - (int) number of projections = N
        - projections
        - (int) numCSFs
        - CSFs: double* (numCSFs * N)
 */
void AngularDataLibrary::Read()
{
    if(filename.empty())
        return;

    FILE* fp = fopen(filename.c_str(), "rb");
    if(!fp)
        return;

    int num_angular_data_objects = 0;

    fread(&num_angular_data_objects, sizeof(int), 1, fp);

    int count = 0;

    while(count < num_angular_data_objects)
    {
        // Key
        int key_size = 0;
        fread(&key_size, sizeof(int), 1, fp);

        KeyType key;    // vector of pair(kappa, num particles)
        key.reserve(key_size);
        for(int i = 0; i < key_size; i++)
        {
            int kappa, num_particles;
            fread(&kappa, sizeof(int), 1, fp);
            fread(&num_particles, sizeof(int), 1, fp);

            key.push_back(std::make_pair(kappa, num_particles));
        }

        pAngularData ang(new AngularData(two_m));
        ang->two_j = two_j;

        // Projections
        int num_projections = 0;
        fread(&num_projections, sizeof(int), 1, fp);

        int particle_number = 0;
        fread(&particle_number, sizeof(int), 1, fp);
        int projection_array[particle_number];

        for(int i = 0; i < num_projections; i++)
        {
            fread(projection_array, sizeof(int), particle_number, fp);
            ang->projections.push_back(std::vector<int>(projection_array, projection_array + particle_number));
        }

        // CSFs
        fread(&ang->num_CSFs, sizeof(int), 1, fp);
        if(ang->num_CSFs)
        {   ang->CSFs = new double[num_projections * ang->num_CSFs];
            fread(ang->CSFs, sizeof(double), ang->num_CSFs * num_projections, fp);
        }
        ang->have_CSFs = true;

        library[key] = ang;
        count++;
    }

    fclose(fp);
}

void AngularDataLibrary::Write() const
{
    if(filename.empty())
        return;

    FILE* fp = fopen(filename.c_str(), "wb");

    if(!fp)
    {   *errstream << "AngularData::Couldn't open file " << filename << " for writing." << std::endl;
    }

    int num_angular_data_objects = library.size();
    fwrite(&num_angular_data_objects, sizeof(int), 1, fp);

    for(const auto& pair: library)
    {
        // Key
        int key_size = pair.first.size();
        fwrite(&key_size, sizeof(int), 1, fp);

        for(auto& key_pair: pair.first)
        {
            fwrite(&key_pair.first, sizeof(int), 1, fp);
            fwrite(&key_pair.second, sizeof(int), 1, fp);
        }

        // Projections
        int num_projections = pair.second->projections.size();
        fwrite(&num_projections, sizeof(int), 1, fp);

        int particle_number = 0;
        if(num_projections)
            particle_number = pair.second->projections.front().size();
        fwrite(&particle_number, sizeof(int), 1, fp);

        for(const auto& projection: pair.second->projections)
        {
            fwrite(&projection[0], sizeof(int), particle_number, fp);
        }

        // CSFs
        fwrite(&pair.second->num_CSFs, sizeof(int), 1, fp);
        if(pair.second->num_CSFs)
            fwrite(pair.second->CSFs, sizeof(double), pair.second->num_CSFs * num_projections, fp);
    }

    fclose(fp);
}

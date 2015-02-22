#include "AngularData.h"
#include "Include.h"
#include "RelativisticConfiguration.h"
#include "Eigen/Eigen"
#include "ManyBodyOperator.h"
#include <numeric>
#include <dirent.h>
#ifdef _MPI
#include <mpi.h>
#endif

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

AngularData::AngularData(const AngularData& other):
    two_m(other.two_m), two_j(other.two_j), num_CSFs(other.num_CSFs), have_CSFs(other.have_CSFs), projections(other.projections)
{
    if(have_CSFs)
    {
        CSFs = new double[projection_size() * num_CSFs];
        memcpy(CSFs, other.CSFs, projection_size() * num_CSFs * sizeof(double));
    }
}

AngularData::AngularData(AngularData&& other):
    two_m(other.two_m), two_j(other.two_j), num_CSFs(other.num_CSFs), have_CSFs(other.have_CSFs), projections(other.projections)
{
    if(have_CSFs)
    {
        CSFs = other.CSFs;
        other.have_CSFs = false;
        other.CSFs = nullptr;
    }
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

    int particle_number = config.ParticleNumber();

    // Skip if vacuum
    if(particle_number == 0)
    {
        projections.clear();
        projections.push_back(std::vector<int>());
        return 1;
    }

    // Create set of boxes, one for each orbital, each with a list of projections to be filled.
    std::vector<int> box_number(particle_number);
    std::vector<int> min_projection(particle_number);
    std::vector<int> max_projection(particle_number);

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

                for(int j = i; j < particle_number; j++)
                {
                    if(j == particle_number-1 ||            // Last particle or
                       box_number[j] != box_number[j+1] ||  // next particle is in a different orbital or
                       existing_proj[j] - existing_proj[j+1] > 2)   // there is room to apply J^-
                    {
                        std::vector<int> new_proj(existing_proj);
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
    Eigen::MatrixXd M = Eigen::MatrixXd::Zero(N, N);

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
            M(i, j) = M(j, i) = matrix_element;
            j_it++; j++;
        }

        i_it++; i++;
    }

    // Solve the matrix
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(M);
    const Eigen::MatrixXd& eigenvectors = es.eigenvectors();
    const Eigen::VectorXd& V = es.eigenvalues();

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
                    CSFs[j * num_CSFs + count] = eigenvectors(j, i);

                count++;
            }
        }
    }

    return num_CSFs;
}

void AngularData::LadderLowering(const RelativisticConfiguration& config, const AngularData& parent)
{
    unsigned int N = projections.size();
    this->two_j = parent.two_j;

    // Clear existing
    if(CSFs)
    {   delete[] CSFs;
        CSFs = nullptr;
        have_CSFs = false;
    }

    if(N == 0 || two_j < abs(two_m))
        return;

    // Create new CSFs and initialise to zero
    num_CSFs = parent.num_CSFs;
    CSFs = new double[N * num_CSFs];
    memset(CSFs, 0, sizeof(double) * N * num_CSFs);

    // Create mapping from projections to CSF index
    std::map<std::vector<int>, int> projection_positions;
    int index = 0;
    for(auto proj_it = projections.begin(); proj_it != projections.end(); proj_it++)
    {   projection_positions[*proj_it] = index;
        index++;
    }

    int particle_number = config.ParticleNumber();
    
    // Create set of boxes, one for each orbital, each with a list of projections to be filled.
    std::vector<int> box_number(particle_number);
    std::vector<int> min_projection(particle_number);
    std::vector<int> max_projection(particle_number);
    std::vector<int> box_twoj(particle_number);
    std::vector<bool> box_ishole(particle_number);

    int count = 0;
    int box = 0;
    for(auto config_it = config.begin(); config_it != config.end(); config_it++)
    {
        int current_two_m = config_it->first.TwoJ();
        int num_particles = abs(config_it->second);
        
        for(int i = 0; i < num_particles; i++)
        {
            max_projection[count + i] = current_two_m;
            min_projection[count + num_particles-1 - i] = -current_two_m;
            box_number[count + i] = box;
            box_twoj[count + i] = config_it->first.TwoJ();
            box_ishole[count + i] = (config_it->second < 0);
            current_two_m -= 2;
        }
        
        count += num_particles;
        box++;
    }

    // For each projection in parent, apply J- for all orbitals, find corresponding projection,
    // and update CSFs
    int parent_index = 0;
    for(auto& parent_proj: parent.projections)
    {
        for(int i = 0; i < particle_number; i++)
        {
            int parent_twom = parent_proj[i];

            if(parent_twom > min_projection[i] &&
               (i == particle_number-1 ||               // Last particle or
                box_number[i] != box_number[i+1] ||     // next particle is in a different orbital or
                parent_proj[i] - parent_proj[i+1] > 2)) // there is room to apply J^-
            {
                std::vector<int> new_proj(parent_proj);
                new_proj[i] -= 2;
                auto it = projection_positions.find(new_proj);
                if(it != projection_positions.end())
                {
                    int child_index = it->second;
                    double prefactor = sqrt((box_twoj[i] + parent_twom) * (box_twoj[i] - parent_twom + 2))/2.;
                    if(box_ishole[i])
                        prefactor = -prefactor;

                    for(int j = 0; j < num_CSFs; j++)
                    {
                        CSFs[child_index * num_CSFs + j] += prefactor * parent.CSFs[parent_index * num_CSFs + j];
                    }
                }
                else
                    *errstream << "AngularData::LadderLowering(): failed to find child projection." << std::endl;
            }
        }

        parent_index++;
    }

    // Renormalise CSFs
    double prefactor = 2./sqrt((two_j + two_m + 2) * (two_j - two_m));
    for(int i = 0; i < num_CSFs * N; i++)
        CSFs[i] *= prefactor;

    have_CSFs = true;
}

AngularDataLibrary::AngularDataLibrary(int electron_number, const Symmetry& sym, int two_m, const std::string& lib_directory):
    electron_number(electron_number), sym(sym), two_m(two_m), lib_directory(lib_directory), parent(nullptr)
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

AngularDataLibrary::KeyType AngularDataLibrary::GenerateKey(const RelativisticConfiguration& config)
{
    KeyType key;
    for(auto& config_it: config)
    {   key.push_back(std::make_pair(config_it.first.Kappa(), config_it.second));
    }
    return key;
}

RelativisticConfiguration AngularDataLibrary::GenerateRelConfig(const KeyType& key)
{
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

    return rconfig;
}

void AngularDataLibrary::GenerateCSFs()
{
    // Generate from ladder operators if M < J
    if(two_m < sym.GetTwoJ())
    {
        int parent_twom = two_m + 2;
        if(parent == nullptr)
        {   parent = pAngularDataLibrary(new AngularDataLibrary(electron_number, sym, parent_twom, lib_directory));
            parent->Read();
        }

        // Make sure parent has all necessary keys
        for(auto& pair: library)
        {
            if(pair.second->CSFs_calculated() == false)
            {
                auto it = parent->library.find(pair.first);
                if(it == parent->library.end())
                {
                    pAngularData ang(new AngularData(GenerateRelConfig(pair.first), parent_twom));
                    parent->library[pair.first] = ang;
                }
            }
        }

        // Generate parent CSFs (recursive)
        parent->GenerateCSFs();

        // Use ladder operators
        for(auto& pair: library)
        {
            if(pair.second->CSFs_calculated() == false)
            {
                pAngularDataConst parent_angular_data = parent->library[pair.first];
                pair.second->LadderLowering(GenerateRelConfig(pair.first), *parent_angular_data);
            }
        }
    }
    else
    {
    #ifndef _MPI
        for(auto& pair: library)
        {   auto& pAng = pair.second;
            if(pAng->CSFs_calculated() == false)
            {
                pAng->GenerateCSFs(GenerateRelConfig(pair.first), sym.GetTwoJ());
            }
        }
    #else
        // Distribute AngularData objects with lots of projections (large matrix) using MPI
        const unsigned int SHARING_SIZE_LIM = 200;

        std::set<std::pair<KeyType, pAngularData>, ProjectionSizeFirstComparator> big_library;

        for(auto& pair: library)
        {
            auto& pAng = pair.second;
            if(pAng->CSFs_calculated() == false)
            {
                if(pAng->projection_size() < SHARING_SIZE_LIM)
                    pAng->GenerateCSFs(GenerateRelConfig(pair.first), two_j);
                else
                    big_library.insert(pair);
            }
        }

        // Find big CSFs
        int root = 0;
        for(auto& pair: big_library)
        {
            auto& pAng = pair.second;
            if(root == ProcessorRank)
                pAng->GenerateCSFs(GenerateRelConfig(pair.first), two_j);

            root++;
            if(root >= NumProcessors)
                root = 0;
        }

        // Share CSFs with all processors
        root = 0;
        int ierr = 0;
        for(auto& pair: big_library)
        {
            auto& pAng = pair.second;
            ierr = MPI_Barrier(MPI_COMM_WORLD);
            ierr = MPI_Bcast(&(pAng->num_CSFs), 1, MPI_INT, root, MPI_COMM_WORLD);

            int buffer_size = pAng->num_CSFs * pAng->projection_size();

            // Allocate space
            if(root != ProcessorRank)
            {
                if(pAng->CSFs)
                    delete[] pAng->CSFs;

                pAng->CSFs = new double[buffer_size];
            }

            ierr = MPI_Bcast(pAng->CSFs, buffer_size, MPI_DOUBLE, root, MPI_COMM_WORLD);
            pAng->have_CSFs = true;

            root++;
            if(root >= NumProcessors)
                root = 0;
        }
    #endif
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
        ang->two_j = sym.GetTwoJ();

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
    if(filename.empty() || ProcessorRank != 0)
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

    // Write parent (recursive)
    if(parent)
        parent->Write();
}

void AngularDataLibrary::PrintKeys() const
{
    for(const auto& pair: library)
    {
        *outstream << GenerateRelConfig(pair.first) << ": projection size = " << pair.second->projection_size();
        if(pair.second->have_CSFs)
            *outstream << "; num CSFs = " << pair.second->num_CSFs << std::endl;
    }
}

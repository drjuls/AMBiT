#include "AngularData.h"
#include "Include.h"
#include "RelativisticConfiguration.h"
#include <Eigen/Eigen>
#include "ManyBodyOperator.h"
#include <boost/interprocess/sync/file_lock.hpp>
#include <boost/interprocess/sync/scoped_lock.hpp>
#include <boost/interprocess/sync/sharable_lock.hpp>
#include <numeric>
#ifdef AMBIT_USE_MPI
    #include <mpi.h>
#endif

namespace Ambit
{
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
    int reserve_size = 0;
    for(auto& list: boxes)
        reserve_size += list.size();
    projections.reserve(reserve_size);

    for(auto& list: boxes)
    {   projections.insert(projections.end(), list.begin(), list.end());
    }
    std::sort(projections.begin(), projections.end(), ProjectionCompare);
    std::unique(projections.begin(), projections.end());
    projections.shrink_to_fit();

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

    // Make vector of Projections.
    std::vector<Projection> real_Projection_list;
    real_Projection_list.reserve(projections.size());
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


AngularDataLibrary::AngularDataLibrary(const std::string& lib_directory):
    directory(lib_directory)
{
    if(!directory.empty())
    {
        // Check library directory exists
        if(!boost::filesystem::exists(directory) || !boost::filesystem::is_directory(directory))
        {
            // Attempt to create, else fail
            if(!boost::filesystem::create_directory(directory))
            {   *errstream << "AngularDataLibarary::AngularDataLibarary() cannot create directory " << directory << std::endl;
                exit(1);
            }
        }

        directory = boost::filesystem::canonical(directory);
    }
}

pAngularData AngularDataLibrary::GetData(const RelativisticConfiguration& config, const Symmetry& sym, int two_m)
{
    KeyType mykey(GenerateKey(config, sym, two_m));
    return GetData(mykey);
}

pAngularData AngularDataLibrary::GetData(const KeyType& key)
{
    pAngularData& ret = library[key];
    if(ret != nullptr)
        return ret;

    auto file_info_key = std::make_tuple(GetElectronNumber(key), key[0].first, key[0].second);
    if(file_info[file_info_key].second.empty())
    {
        Read(file_info_key);
    }

    if(ret != nullptr)
        return ret;

    ret = std::make_shared<AngularData>(GenerateRelConfig(key), key[0].second);
    return ret;
}

AngularDataLibrary::KeyType AngularDataLibrary::GenerateKey(const RelativisticConfiguration& config, const Symmetry& sym, int two_m)
{
    KeyType key;
    key.reserve(config.size() + 1);
    key.push_back(std::make_pair(sym.GetJpi(), two_m));
    for(auto& config_it: config)
    {   key.push_back(std::make_pair(config_it.first.Kappa(), config_it.second));
    }
    return key;
}

RelativisticConfiguration AngularDataLibrary::GenerateRelConfig(const KeyType& key)
{
    // Create equivalent RelativisticConfiguration and use it.
    RelativisticConfiguration rconfig;
    int pqn = 1;

    auto it = key.begin();
    it++;   // Skip pair<symmetry.Jpi, two_m>

    while(it != key.end())
    {   // Set PQN to some integer
        rconfig.insert(std::make_pair(OrbitalInfo(pqn, it->first), it->second));
        pqn++;
        it++;
    }

    return rconfig;
}

int AngularDataLibrary::GetElectronNumber(const KeyType& key)
{
    int electron_number = 0;
    for(int i = 1; i < key.size(); i++)
        electron_number += key[i].second;
    return electron_number;
}

void AngularDataLibrary::GenerateCSFs()
{
    // First, for all keys with no CSFs, generate entries in the library with higher M,
    // recursively until CSFs are found or M = J.
    for(auto& pair: library)
    {
        pAngularData pAng = pair.second;
        Symmetry sym(pair.first[0].first);
        int two_m = pair.first[0].second;

        // If no CSFs and M < J
        while(pAng->CSFs_calculated() == false && two_m < sym.GetTwoJ())
        {
            KeyType new_key = pair.first;
            two_m += 2;
            new_key[0].second = two_m;
            pAng = GetData(new_key);
        }
    }

    // Get CSFs for M = J
#ifndef AMBIT_USE_MPI
    for(auto& pair: library)
    {
        Symmetry sym(pair.first[0].first);
        int two_m = pair.first[0].second;
        auto& pAng = pair.second;

        if(pAng->CSFs_calculated() == false && sym.GetTwoJ() == two_m)
        {
            RelativisticConfiguration rconfig(GenerateRelConfig(pair.first));
            pAng->GenerateCSFs(rconfig, two_m);

            // Set write_needed to true;
            file_info[std::make_tuple(GetElectronNumber(pair.first), sym.GetJpi(), two_m)].first = true;
        }
    }
#else
    // Distribute AngularData objects with lots of projections (large matrix) using MPI
    const unsigned int SHARING_SIZE_LIM = 200;

    std::set<std::pair<KeyType, pAngularData>, ProjectionSizeFirstComparator> big_library;

    for(auto& pair: library)
    {
        Symmetry sym(pair.first[0].first);
        int two_m = pair.first[0].second;
        auto& pAng = pair.second;

        if(pAng->CSFs_calculated() == false && sym.GetTwoJ() == two_m)
        {
            RelativisticConfiguration rconfig(GenerateRelConfig(pair.first));

            if(pAng->projection_size() < SHARING_SIZE_LIM)
                pAng->GenerateCSFs(rconfig, sym.GetTwoJ());
            else
                big_library.insert(pair);

            // Set write_needed to true;
            file_info[std::make_tuple(GetElectronNumber(pair.first), sym.GetJpi(), two_m)].first = true;
        }
    }

    // Find big CSFs with M = J
    int jobcount = 0;   // jobcount will run from 0 -> 2 * NumProcessors - 1
    int otherRank = 2 * NumProcessors - 1 - ProcessorRank;

    for(auto& pair: big_library)
    {
        Symmetry sym(pair.first[0].first);
        auto& pAng = pair.second;

        if(jobcount == ProcessorRank || jobcount == otherRank)
        {
            RelativisticConfiguration rconfig(GenerateRelConfig(pair.first));
            *logstream << "Calculating CSF: N = " << pAng->projection_size() << " " << rconfig << std::endl;

            pAng->GenerateCSFs(rconfig, sym.GetTwoJ());
        }

        jobcount++;
        if(jobcount >= 2 * NumProcessors)
            jobcount = 0;
    }

    // Share CSFs with all processors
    jobcount = 0;
    int ierr = 0;
    for(auto& pair: big_library)
    {
        int root = jobcount;
        if(root >= NumProcessors)
            root = 2 * NumProcessors - 1 - root;

        auto& pAng = pair.second;
        ierr = MPI_Barrier(MPI_COMM_WORLD);
        ierr = MPI_Bcast(&(pAng->num_CSFs), 1, MPI_INT, root, MPI_COMM_WORLD);

        int buffer_size = pAng->num_CSFs * pAng->projection_size();

        // Allocate space
        if(jobcount != ProcessorRank && jobcount != otherRank)
        {
            if(pAng->CSFs)
                delete[] pAng->CSFs;

            pAng->CSFs = new double[buffer_size];
        }

        ierr = MPI_Bcast(pAng->CSFs, buffer_size, MPI_DOUBLE, root, MPI_COMM_WORLD);
        pAng->have_CSFs = true;

        jobcount++;
        if(jobcount >= 2 * NumProcessors)
            jobcount = 0;
    }
#endif

    // Get the rest recursively using ladder operators
    for(auto& pair: library)
    {
        if(pair.second->CSFs_calculated() == false)
        {
            KeyType new_key = pair.first;
            RelativisticConfiguration rconfig(GenerateRelConfig(pair.first));

            pAngularData pAng = pair.second;
            Symmetry sym(pair.first[0].first);
            int two_m = pair.first[0].second;

            // If no CSFs and M < J, get first parent with CSFs
            while(pAng->CSFs_calculated() == false && two_m < sym.GetTwoJ())
            {
                two_m += 2;
                new_key[0].second = two_m;
                pAng = GetData(new_key);
            }

            // Work backwards to current pair
            while(two_m > pair.first[0].second)
            {
                pAngularData parent = pAng;
                two_m -= 2;
                new_key[0].second = two_m;
                pAng = GetData(new_key);
                pAng->LadderLowering(rconfig, *parent);

                // Set write_needed to true;
                file_info[std::make_tuple(GetElectronNumber(pair.first), sym.GetJpi(), two_m)].first = true;

                pAng = parent;
            }
        }
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
void AngularDataLibrary::Read(int electron_number, const Symmetry& sym, int two_m)
{
    if(directory.empty())
        return;

    auto& filedata = file_info[std::make_tuple(electron_number, sym.GetJpi(), two_m)];
    boost::filesystem::path& filepath = filedata.second;

    if(filepath.empty())
    {
        std::string holelike = (electron_number<0)?"h":"";
        std::string filename = itoa(abs(electron_number)) + holelike + "." + sym.GetString() + "." + itoa(two_m) + ".angular";

        filepath = directory / filename;
    }

    // Only try to get a shared file lock if the file actually exists, otherwise there's no point reading
    if(!boost::filesystem::exists(filedata.second.c_str()))
        return;

    boost::interprocess::file_lock f_lock(filedata.second.c_str());
    boost::interprocess::sharable_lock<boost::interprocess::file_lock> shlock(f_lock);

    FILE* fp = file_err_handler->fopen(filepath.string().c_str(), "rb");
    if(!fp)
        return;

    int num_angular_data_objects = 0;
    file_err_handler->fread(&num_angular_data_objects, sizeof(int), 1, fp);

    int count = 0;

    while(count < num_angular_data_objects)
    {
        // Key
        int key_size = 0;
        file_err_handler->fread(&key_size, sizeof(int), 1, fp);

        KeyType key;    // vector of pair(kappa, num particles)
        key.reserve(key_size);
        for(int i = 0; i < key_size; i++)
        {
            int kappa, num_particles;
            file_err_handler->fread(&kappa, sizeof(int), 1, fp);
            file_err_handler->fread(&num_particles, sizeof(int), 1, fp);

            key.push_back(std::make_pair(kappa, num_particles));
        }

        pAngularData ang(new AngularData(two_m));
        ang->two_j = sym.GetTwoJ();

        // Projections
        int num_projections = 0;
        file_err_handler->fread(&num_projections, sizeof(int), 1, fp);

        int particle_number = 0;
        file_err_handler->fread(&particle_number, sizeof(int), 1, fp);
        int projection_array[particle_number];

        for(int i = 0; i < num_projections; i++)
        {
            file_err_handler->fread(projection_array, sizeof(int), particle_number, fp);
            ang->projections.push_back(std::vector<int>(projection_array, projection_array + particle_number));
        }

        // CSFs
        file_err_handler->fread(&ang->num_CSFs, sizeof(int), 1, fp);
        if(ang->num_CSFs)
        {   ang->CSFs = new double[num_projections * ang->num_CSFs];
            file_err_handler->fread(ang->CSFs, sizeof(double), ang->num_CSFs * num_projections, fp);
        }
        ang->have_CSFs = true;

        library[key] = ang;
        count++;
    }

    file_err_handler->fclose(fp);
    filedata.first = false;
}

void AngularDataLibrary::Read(const std::tuple<int, int, int>& file_info_key)
{
    return Read(std::get<0>(file_info_key), Symmetry(std::get<1>(file_info_key)), std::get<2>(file_info_key));
}

void AngularDataLibrary::Write(int electron_number, const Symmetry& sym, int two_m)
{
    if(directory.empty())
        return;

    auto file_info_it = file_info.find(std::make_tuple(electron_number, sym.GetJpi(), two_m));

    // If symmetry not found or write not needed, stop.
    if(file_info_it == file_info.end() || file_info_it->second.first == false)
        return;

#ifdef AMBIT_USE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif

    if(ProcessorRank == 0)
    {
        // Read any extra existing AngularData objects
        auto& filedata = file_info_it->second;
        if(filedata.second.empty())
            Read(electron_number, sym, two_m);

        FILE* fp = nullptr;
        // Create the AngularData file if it doesn't already exist
        if(!boost::filesystem::exists(filedata.second.c_str()))
        {
            fp = file_err_handler->fopen(filedata.second.c_str(), "wb");
            if(!fp)
            {   *errstream << "AngularDataLibrary::Couldn't open file " << filedata.second << " for writing." << std::endl;
                return;
            }
        }

        // Wait for exclusive file lock
        boost::interprocess::file_lock f_lock(filedata.second.c_str());
        boost::interprocess::scoped_lock<boost::interprocess::file_lock> shylock(f_lock);

        // Open file
        if(!fp)
        {
            fp = file_err_handler->fopen(filedata.second.c_str(), "wb");
            if(!fp)
            {   *errstream << "AngularDataLibrary::Couldn't open file " << filedata.second << " for writing." << std::endl;
                return;
            }
        }

        // Count number of elements with same symmetry
        int num_angular_data_objects = 0;
        auto symmetry_pair = std::make_pair(sym.GetJpi(), two_m);
        for(auto& pair: library)
        {
            if(pair.first[0] == symmetry_pair)
                num_angular_data_objects++;
        }
        file_err_handler->fwrite(&num_angular_data_objects, sizeof(int), 1, fp);

        for(const auto& pair: library)
        {
            if(pair.first[0] == symmetry_pair)
            {
                // Key
                int key_size = pair.first.size();
                file_err_handler->fwrite(&key_size, sizeof(int), 1, fp);

                for(auto& key_pair: pair.first)
                {
                    file_err_handler->fwrite(&key_pair.first, sizeof(int), 1, fp);
                    file_err_handler->fwrite(&key_pair.second, sizeof(int), 1, fp);
                }

                // Projections
                int num_projections = pair.second->projections.size();
                file_err_handler->fwrite(&num_projections, sizeof(int), 1, fp);

                int particle_number = 0;
                if(num_projections)
                    particle_number = pair.second->projections.front().size();
                file_err_handler->fwrite(&particle_number, sizeof(int), 1, fp);

                for(const auto& projection: pair.second->projections)
                {
                    file_err_handler->fwrite(&projection[0], sizeof(int), particle_number, fp);
                }

                // CSFs
                file_err_handler->fwrite(&pair.second->num_CSFs, sizeof(int), 1, fp);
                if(pair.second->num_CSFs)
                    file_err_handler->fwrite(pair.second->CSFs, sizeof(double), pair.second->num_CSFs * num_projections, fp);
            }
        }
        file_err_handler->fclose(fp);
    }

#ifdef AMBIT_USE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif

    // Set write_needed to false
    file_info_it->second.first = false;
}

void AngularDataLibrary::Write()
{
    if(directory.empty())
        return;

    for(auto& filedata: file_info)
    {
        if(filedata.second.first)
            Write(std::get<0>(filedata.first), Symmetry(std::get<1>(filedata.first)), std::get<2>(filedata.first));
    }
}

void AngularDataLibrary::RemoveUnused()
{
    auto it = library.begin();
    while(it != library.end())
    {
        // Remove if library holds the only pointer to the AngularData object.
        if(it->second.use_count() == 1)
        {
            // Remove filename (force Read() next time)
            auto file_info_key = std::make_tuple(GetElectronNumber(it->first), it->first[0].first, it->first[0].second);
            file_info[file_info_key].second.clear();
            it = library.erase(it);
        }
        else
            it++;
    }
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
}

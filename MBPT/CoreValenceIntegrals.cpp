#include "Include.h"
#include <boost/filesystem.hpp>
#include <chrono>
#ifdef AMBIT_USE_MPI
#include <mpi.h>
#endif

// Below purposely not included: this file is for a template class and should be included in the header.
// #include "CoreValenceIntegrals.h"

namespace Ambit
{
template <class MapType>
CoreValenceIntegrals<MapType>::CoreValenceIntegrals(pOrbitalManagerConst orbitals, pHFIntegrals one_body, pHartreeY hartreeY_op, const std::string& write_file):
    CoreValenceIntegrals(one_body, pSlaterIntegrals(new SlaterIntegrals<MapType>(orbitals, hartreeY_op)), write_file)
{}

template <class MapType>
CoreValenceIntegrals<MapType>::CoreValenceIntegrals(pOrbitalManagerConst orbitals, pHFIntegrals one_body, pSlaterIntegrals bare_integrals, const std::string& write_file):
    SlaterIntegrals<MapType>(orbitals, false), write_file(write_file), core_PT(nullptr),
    include_core(true), include_core_subtraction(true), include_core_extra_box(true),
    include_valence(false), include_valence_subtraction(false), include_valence_extra_box(false)
{
    core_PT.reset(new CoreMBPTCalculator(this->orbitals, one_body, bare_integrals));
    valence_PT.reset(new ValenceMBPTCalculator(this->orbitals, one_body, bare_integrals));
}

template <class MapType>
CoreValenceIntegrals<MapType>::CoreValenceIntegrals(pOrbitalManagerConst orbitals, pCoreMBPTCalculator core_mbpt_calculator, pValenceMBPTCalculator valence_mbpt_calculator, const std::string& write_file):
    SlaterIntegrals<MapType>(orbitals, false), write_file(write_file),
    core_PT(core_mbpt_calculator), valence_PT(valence_mbpt_calculator),
    include_core(false), include_core_subtraction(false), include_core_extra_box(false),
    include_valence(false), include_valence_subtraction(false), include_valence_extra_box(false)
{
    if(core_PT)
    {   include_core = true;
        include_core_subtraction = true;
        include_core_extra_box = true;
    }
    if(valence_PT)
    {   include_valence = true;
        include_valence_subtraction = true;
        include_valence_extra_box = true;
    }
}

template <class MapType>
CoreValenceIntegrals<MapType>::~CoreValenceIntegrals()
{}

template <class MapType>
unsigned int CoreValenceIntegrals<MapType>::CalculateTwoElectronIntegrals(pOrbitalMapConst orbital_map_1, pOrbitalMapConst orbital_map_2, pOrbitalMapConst orbital_map_3, pOrbitalMapConst orbital_map_4, bool check_size_only)
{
    unsigned int i1, i2, i3, i4;
    int k, kmax;

    std::set<KeyType> previous_keys;   // Keep track of integrals we've already read from a file
    for(auto integral : this->TwoElectronIntegrals)
        previous_keys.insert(integral.first);
#ifdef AMBIT_USE_MPI
    int count = 0;  // Processor index

    my_calculations_done = false;
    root_complete = false;
#endif

    // NOTE: these have to be local variables because OpenMP doesn't support class-members in
    // data-sharing clauses. This means we need to pass these in as arguments to Write()
    std::vector<std::tuple<int, unsigned, unsigned, unsigned, unsigned>> expanded_keys;
    std::vector<KeyType> keys;
    std::vector<double> values;

    if(!check_size_only)
    {
        if(include_core || include_core_subtraction || include_core_extra_box)
            core_PT->UpdateIntegrals();
        if(include_valence || include_valence_subtraction || include_valence_extra_box)
            valence_PT->UpdateIntegrals();
    }
    
    // Populate a table of key-value pairs (basically a vector of pairs). This doesn't actually
    // calculate the integrals
    auto it_1 = orbital_map_1->begin();
    while(it_1 != orbital_map_1->end())
    {
        i1 = this->orbitals->state_index.at(it_1->first);
        const auto& s1 = it_1->first;

        auto it_3 = orbital_map_3->begin();
        while(it_3 != orbital_map_3->end())
        {
            i3 = this->orbitals->state_index.at(it_3->first);
            const auto& s3 = it_3->first;

            auto it_2 = orbital_map_2->begin();
            while(it_2 != orbital_map_2->end())
            {
                i2 = this->orbitals->state_index.at(it_2->first);
                const auto& s2 = it_2->first;

                auto it_4 = orbital_map_4->begin();
                while(it_4 != orbital_map_4->end())
                {
                    i4 = this->orbitals->state_index.at(it_4->first);
                    const auto& s4 = it_4->first;

                    // Check parity conservation
                    if((s1.L() + s2.L() + s3.L() + s4.L())%2 == 0)
                    {
                        // Limits on k
                        k = mmax(abs(s1.L() - s3.L()), abs(s2.L() - s4.L()));
                        if((abs(s1.TwoJ() - s3.TwoJ()) > 2 * k) || abs(s2.TwoJ() - s4.TwoJ()) > 2 * k)
                        {
                            if(include_core_extra_box || include_valence_extra_box)
                                k++;
                            else
                                k+=2;
                        }

                        kmax = mmin(s1.L() + s3.L(), s2.L() + s4.L());
                        if((s1.TwoJ() + s3.TwoJ() <  2 * kmax) || (s2.TwoJ() + s4.TwoJ() <  2 * kmax))
                        {
                            if(include_core_extra_box || include_valence_extra_box)
                                kmax--;
                            else
                                kmax-=2;
                        }

                        while(k <= kmax)
                        {
                            // Usual multipolarity rules
                            bool usual_parity = ((s2.L() + s4.L() + k)%2 == 0);
                            if(usual_parity || include_core_extra_box || include_valence_extra_box)
                            {
                                KeyType key = this->GetKey(k, i1, i2, i3, i4);

                                    if(!previous_keys.count(key) && (include_core || include_core_subtraction || include_valence || include_valence_subtraction || !usual_parity))
                                    {
                                        // Keep track of the key even if it's not ours. This is necessary
                                        // for check-sizes to work correctly across multiple processes.
                                        previous_keys.insert(key);
                                    #ifdef AMBIT_USE_MPI
                                        // Check if this is our job
                                        if(count == ProcessorRank || check_size_only)
                                        {
                                    #endif
                                            keys.push_back(key);
                                            // We need to store the actual orbital indices as well as the keys, since translating 
                                            // between key and orbitals is not its own inverse (i.e. it doesn't preserve the order of
                                            // the orbitals in the integral) and this can slightly affect the result. This will 
                                            // increase the memory consumption within this subroutine, but these arrays will be 
                                            // freed once we move on, si the overall effect on memory footprint should be small.
                                            // TODO: Check with Julian as to whether the ordering of the orbitals *should* be 
                                            // significant or is this a bug?
                                            expanded_keys.push_back(std::make_tuple<int, unsigned, unsigned, unsigned, unsigned>\
                                                    (std::move(k), std::move(i1), std::move(i2), std::move(i3), std::move(i4)));

                                    #ifdef AMBIT_USE_MPI
                                        }
                                            count++;
                                            if(count == NumProcessors)
                                                count = 0;
                                    #endif
                                    }
                            }

                            if(include_core_extra_box || include_valence_extra_box)
                                k++;
                            else
                                k+=2;
                        }
                    }
                    it_4++;
                }
                it_2++;
            }
            it_3++;
        }
        it_1++;
    }
    
    if(check_size_only)
        return(previous_keys.size());

    // Make sure to allocate enough storage for the calculated integrals
    values.resize(expanded_keys.size());

    // Now run through all of this process's keys and calculate the corresponding integrals
#ifdef AMBIT_USE_OPENMP
    #pragma omp parallel for default(none) shared(expanded_keys, keys, values, stderr) schedule(dynamic, 8)
#endif
    for(int n = 0; n < expanded_keys.size(); n++)
    {   // Expand the key into a set of orbital indices and multipolarity k
        auto expanded_key = expanded_keys[n];
        int k = std::get<0>(expanded_key);
        unsigned int i1 = std::get<1>(expanded_key);
        // reverse_state_index is a typedef for std::vector<OrbitalInfo>
        const auto& s1 = this->orbitals->reverse_state_index.at(i1);

        unsigned int i2 = std::get<2>(expanded_key);
        const auto& s2 = this->orbitals->reverse_state_index.at(i2);

        unsigned int i3 = std::get<3>(expanded_key);
        const auto& s3 = this->orbitals->reverse_state_index.at(i3);

        unsigned int i4 = std::get<4>(expanded_key);
        const auto& s4 = this->orbitals->reverse_state_index.at(i4);
     
        // Now actually calculate the various diagrams
        double radial = 0;
        bool usual_parity = ((s2.L() + s4.L() + k)%2 == 0);

        if(usual_parity)
        {   if(include_core)
                radial += core_PT->GetTwoElectronDiagrams(k, s1, s2, s3, s4);
            if(include_core_subtraction)
                radial += core_PT->GetTwoElectronSubtraction(k, s1, s2, s3, s4);
            if(include_valence)
                radial += valence_PT->GetTwoElectronValence(k, s1, s2, s3, s4);
            if(include_valence_subtraction)
                radial += valence_PT->GetTwoElectronSubtraction(k, s1, s2, s3, s4);
        }
        else
        {   if(include_core_extra_box)
                radial += core_PT->GetTwoElectronBoxDiagrams(k, s1, s2, s3, s4);
            if(include_valence_extra_box)
                radial += valence_PT->GetTwoElectronBoxValence(k, s1, s2, s3, s4);
        }

        values[n] = radial;
    }

    #ifdef AMBIT_USE_MPI
    if(!check_size_only)
    {   // Gather to root node, write to file, and read back in
        my_calculations_done = true;
        if(ProcessorRank == 0)
            root_complete = true;

        // Do loop protects processes that finished before root and therefore
        // may have to write multiple times. No need to explicitly collect the key-value pairs, since
        // this happens when we read the file back in.
        do{
            this->Write(write_file, keys, values);
        } while(!root_complete);

        this->clear();
        keys.clear();
        values.clear();
        this->Read(write_file);

        return this->TwoElectronIntegrals.size();
    }
    #else
    {   
        // Collect all of our key-value pairs and put them in the class' map for writing
        auto key_it = keys.begin();
        auto value_it = values.begin();
        while(key_it != keys.end() && value_it != values.end())
        {   
            this->TwoElectronIntegrals.insert(std::pair<KeyType, double>(*key_it, *value_it));
            key_it++;
            value_it++;
        } 

        this->Write(write_file);
        return this->TwoElectronIntegrals.size();
    }
    #endif
    // We should never reach this point unless there's been a compilation error, so bail with an error
    exit(EXIT_FAILURE);
}

#ifdef AMBIT_USE_MPI
template <class MapType>
void CoreValenceIntegrals<MapType>::Write(const std::string& filename, const std::vector<KeyType>& keys, const std::vector<double>& values) const
{
    unsigned int KeyType_size = sizeof(KeyType);
    MPI_Datatype mpikeytype;
    switch(KeyType_size)
    {
        case 4:
            mpikeytype = MPI_UNSIGNED;
            break;
        case 8:
        default:
            mpikeytype = MPI_UNSIGNED_LONG_LONG;
            break;
    }
    
    if(ProcessorRank == 0)
    {
        // Send completion status
        int complete = (my_calculations_done? 1: 0);
        MPI_Bcast(&complete, 1, MPI_INT, 0, MPI_COMM_WORLD);

        // Get number of integrals from all processes. 
        std::vector<unsigned int> num_integrals(NumProcessors);
        num_integrals[0] = values.size();

        MPI_Status status;
        for(int proc = 1; proc < NumProcessors; proc++)
        {
            MPI_Recv(&num_integrals[proc], 1, MPI_UNSIGNED, proc, 1, MPI_COMM_WORLD, &status);
        }

        // Open after it is clear that all processes are on board
        FILE* fp = file_err_handler->fopen(filename.c_str(), "wb");

        // Write state index
        WriteOrbitalIndexes(this->orbitals->state_index, fp);

        file_err_handler->fwrite(&KeyType_size, sizeof(unsigned int), 1, fp);

        unsigned int total_integrals = std::accumulate(num_integrals.begin(), num_integrals.end(), this->size());
        file_err_handler->fwrite(&total_integrals, sizeof(unsigned int), 1, fp);

        // Write root integrals
        for(auto& pair: this->TwoElectronIntegrals)
        {
            const double value = pair.second;   // Convert to double
            if(!std::isnan(value))
            {
                file_err_handler->fwrite(&pair.first, sizeof(KeyType), 1, fp);
                file_err_handler->fwrite(&value, sizeof(double), 1, fp);
            }
        }

        for(int i = 0; i < num_integrals[0]; i++)
        {
            const double value = values[i];   // Convert to double
            if(!std::isnan(value))
            {
                file_err_handler->fwrite(&keys[i], sizeof(KeyType), 1, fp);
                file_err_handler->fwrite(&value, sizeof(double), 1, fp);
            }
        }

        // Receive and write data from other processes
        std::vector<KeyType> mpi_keys;
        std::vector<double> mpi_values;
        for(int proc = 1; proc < NumProcessors; proc++)
        {
            mpi_keys.resize(num_integrals[proc]);
            mpi_values.resize(num_integrals[proc]);

            MPI_Recv(mpi_keys.data(), num_integrals[proc], mpikeytype, proc, 2, MPI_COMM_WORLD, &status);
            MPI_Recv(mpi_values.data(), num_integrals[proc], MPI_DOUBLE, proc, 3, MPI_COMM_WORLD, &status);

            for(int i = 0; i < num_integrals[proc]; i++)
            {
                if(!std::isnan(mpi_values[i]))
                {
                    file_err_handler->fwrite(&mpi_keys[i], sizeof(KeyType), 1, fp);
                    file_err_handler->fwrite(&mpi_values[i], sizeof(double), 1, fp);
                }
            }
        }

        file_err_handler->fclose(fp);
        *logstream << "Written " << total_integrals << " two-body MBPT integrals." << std::endl;

        // The *.two.int.save file is created to protect the work while
        // the *.two.int file is being re-written.
        // It is removed after the last write.
        boost::filesystem::path current_path(filename);
        boost::filesystem::path save_path(filename + ".save");

        if(boost::filesystem::exists(save_path))
            boost::filesystem::remove(save_path);
        if(!my_calculations_done)
            boost::filesystem::copy(current_path, save_path);
    }
    else
    {   // Receive root completion status
        if(!root_complete)
        {   int buffer;
            MPI_Bcast(&buffer, 1, MPI_INT, 0, MPI_COMM_WORLD);
            root_complete = (buffer != 0);
        }

        // If root is complete, but we are not, continue until our work is complete
        if(root_complete && !my_calculations_done)
            return;

        // Send data to root. First number of keys.
        unsigned int num_integrals = values.size();
        MPI_Send(&num_integrals, 1, MPI_UNSIGNED, 0, 1, MPI_COMM_WORLD);

        // Send keys and data to root
        MPI_Send(keys.data(), num_integrals, mpikeytype, 0, 2, MPI_COMM_WORLD);
        MPI_Send(values.data(), num_integrals, MPI_DOUBLE, 0, 3, MPI_COMM_WORLD);
    }

    MPI_Barrier(MPI_COMM_WORLD);
}
#endif

template <class MapType>
void CoreValenceIntegrals<MapType>::IncludeCore(bool include_mbpt, bool include_subtraction, bool include_wrong_parity_box_diagrams)
{
    include_core = include_mbpt && core_PT;
    include_core_subtraction = include_subtraction && core_PT;
    include_core_extra_box = include_wrong_parity_box_diagrams && core_PT;
}

template <class MapType>
void CoreValenceIntegrals<MapType>::IncludeValence(bool include_mbpt, bool include_subtraction, bool include_wrong_parity_box_diagrams)
{
    include_valence = include_mbpt && valence_PT;
    include_valence_subtraction = include_subtraction && valence_PT;
    include_valence_extra_box = include_wrong_parity_box_diagrams && valence_PT;
}
}

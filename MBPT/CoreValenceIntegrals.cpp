#include "Include.h"
#include <boost/filesystem.hpp>
#include <chrono>
#ifdef AMBIT_USE_MPI
#include <mpi.h>
#endif

// Below purposely not included: this file is for a template class and should be included in the header.
// #include "CoreValenceIntegrals.h"

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

    std::set<KeyType> found_keys;   // For check_size_only or MPI checking

#ifdef AMBIT_USE_MPI
    int count = 0;  // Processor index

    if(!check_size_only)
    {
        int new_keys_per_processor = CalculateTwoElectronIntegrals(orbital_map_1, orbital_map_2, orbital_map_3, orbital_map_4, true);
        new_keys_per_processor = (new_keys_per_processor + NumProcessors - 1)/NumProcessors;

        new_keys.clear();
        new_keys.reserve(new_keys_per_processor);
        new_values.clear();
        new_values.reserve(new_keys_per_processor);

        for(const auto& pair: this->TwoElectronIntegrals)
            found_keys.insert(pair.first);

        my_calculations_done = false;
        root_complete = false;
    }
#endif

    if(!check_size_only)
    {
        if(include_core || include_core_subtraction || include_core_extra_box)
            core_PT->UpdateIntegrals();
        if(include_valence || include_valence_subtraction || include_valence_extra_box)
            valence_PT->UpdateIntegrals();
    }

    // Save state every hour or so unless nearly_done is true
    std::chrono::steady_clock::time_point mark_time = std::chrono::steady_clock::now();
    std::chrono::minutes gap(47);   // 47 minutes
    bool nearly_done = false;

    auto it_1 = orbital_map_1->begin();
    while(it_1 != orbital_map_1->end())
    {
        if(std::next(it_1) == orbital_map_1->end())
            nearly_done = true;

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

                                if(check_size_only)
                                {
                                    if((this->TwoElectronIntegrals.find(key) == this->TwoElectronIntegrals.end()) &&
                                       (include_core || include_core_subtraction || include_valence || include_valence_subtraction || !usual_parity))
                                        found_keys.insert(key);
                                }
                                else
                                {   // Check that this integral doesn't already exist
                                    #ifdef AMBIT_USE_MPI
                                    if(!found_keys.count(key))
                                    {
                                        // Check if this is our job
                                        if(count == ProcessorRank)
                                        {
                                    #else
                                    if(this->TwoElectronIntegrals.find(key) == this->TwoElectronIntegrals.end())
                                    {
                                    #endif

                                        double radial = 0;
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

                                    #ifdef AMBIT_USE_MPI
                                        new_keys.push_back(key);
                                        new_values.push_back(radial);
                                        }

                                        found_keys.insert(key);
                                        count++;
                                        if(count == NumProcessors)
                                            count = 0;
                                    #else
                                        this->TwoElectronIntegrals.insert(std::pair<KeyType, double>(key, radial));
                                    #endif
                                    }
                                }
                            }

                            if(include_core_extra_box || include_valence_extra_box)
                                k++;
                            else
                                k+=2;
                        }
                    }

                    // Save state if lots of time has passed
                    std::chrono::steady_clock::duration t = std::chrono::steady_clock::now() - mark_time;
                    if(!nearly_done && t > gap)
                    {   this->Write(write_file);    // Write has an MPI_Barrier, so all processes are at the same time
                        mark_time = std::chrono::steady_clock::now();
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
        return found_keys.size();
    else
    #ifdef AMBIT_USE_MPI
    {   // Gather to root node, write to file, and read back in
        my_calculations_done = true;
        if(ProcessorRank == 0)
            root_complete = true;

        // Do loop protects processes that finished before root and therefore
        // may have to write multiple times.
        do{
            this->Write(write_file);
        } while(!root_complete);

        this->clear();
        new_keys.clear();
        new_values.clear();
        this->Read(write_file);

        return this->TwoElectronIntegrals.size();
    }
    #else
    {   this->Write(write_file);
        return this->TwoElectronIntegrals.size();
    }
    #endif
}

#ifdef AMBIT_USE_MPI
template <class MapType>
void CoreValenceIntegrals<MapType>::Write(const std::string& filename) const
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

        // Get number of keys from all processes
        std::vector<unsigned int> num_integrals(NumProcessors);
        num_integrals[0] = new_keys.size();

        MPI_Status status;
        for(int proc = 1; proc < NumProcessors; proc++)
        {
            MPI_Recv(&num_integrals[proc], 1, MPI_UNSIGNED, proc, 1, MPI_COMM_WORLD, &status);
        }

        // Open after it is clear that all processes are on board
        FILE* fp = fopen(filename.c_str(), "wb");

        // Write state index
        WriteOrbitalIndexes(this->orbitals->state_index, fp);

        fwrite(&KeyType_size, sizeof(unsigned int), 1, fp);

        unsigned int total_integrals = std::accumulate(num_integrals.begin(), num_integrals.end(), this->size());
        fwrite(&total_integrals, sizeof(unsigned int), 1, fp);

        // Write root integrals
        for(auto& pair: this->TwoElectronIntegrals)
        {
            const double value = pair.second;   // Convert to double
            fwrite(&pair.first, sizeof(KeyType), 1, fp);
            fwrite(&value, sizeof(double), 1, fp);
        }

        for(int i = 0; i < num_integrals[0]; i++)
        {
            const double value = new_values[i];   // Convert to double
            fwrite(&new_keys[i], sizeof(KeyType), 1, fp);
            fwrite(&value, sizeof(double), 1, fp);
        }

        // Receive and write data from other processes
        std::vector<KeyType> keys;
        std::vector<double> values;
        for(int proc = 1; proc < NumProcessors; proc++)
        {
            keys.resize(num_integrals[proc]);
            values.resize(num_integrals[proc]);

            MPI_Recv(keys.data(), num_integrals[proc], mpikeytype, proc, 2, MPI_COMM_WORLD, &status);
            MPI_Recv(values.data(), num_integrals[proc], MPI_DOUBLE, proc, 3, MPI_COMM_WORLD, &status);

            for(int i = 0; i < num_integrals[proc]; i++)
            {
                fwrite(&keys[i], sizeof(KeyType), 1, fp);
                fwrite(&values[i], sizeof(double), 1, fp);
            }
        }

        fclose(fp);
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
        unsigned int num_keys = new_keys.size();
        MPI_Send(&num_keys, 1, MPI_UNSIGNED, 0, 1, MPI_COMM_WORLD);

        // Send keys and data to root
        MPI_Send(new_keys.data(), num_keys, mpikeytype, 0, 2, MPI_COMM_WORLD);
        MPI_Send(new_values.data(), num_keys, MPI_DOUBLE, 0, 3, MPI_COMM_WORLD);
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

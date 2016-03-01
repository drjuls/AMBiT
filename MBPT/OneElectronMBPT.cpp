#include "OneElectronMBPT.h"
#ifdef AMBIT_USE_MPI
#include <mpi.h>
#endif

OneElectronMBPT::OneElectronMBPT(pOrbitalManagerConst orbitals, pHFIntegrals bare_one_body, pSlaterIntegrals bare_two_body, const std::string& write_file):
    OneElectronMBPT(orbitals, bare_one_body->GetOperator(),
                    std::make_shared<CoreMBPTCalculator>(orbitals, bare_one_body, bare_two_body),
                    std::make_shared<ValenceMBPTCalculator>(orbitals, bare_one_body, bare_two_body), write_file)
{}

OneElectronMBPT::OneElectronMBPT(pOrbitalManagerConst orbitals, pSpinorMatrixElementConst pOperator, pCoreMBPTCalculator core_mbpt_calculator, pValenceMBPTCalculator valence_mbpt_calculator, const std::string& write_file):
    OneElectronIntegrals(orbitals, pOperator), core_PT(core_mbpt_calculator), valence_PT(valence_mbpt_calculator), write_file(write_file),
    include_core(false), include_core_subtraction(false), include_valence_subtraction(false)
{
    if(core_PT)
    {   include_core = true;
        include_core_subtraction = true;
    }
    if(valence_PT)
        include_valence_subtraction = true;
}

unsigned int OneElectronMBPT::CalculateOneElectronIntegrals(pOrbitalMapConst orbital_map_1, pOrbitalMapConst orbital_map_2, bool check_size_only)
{
    unsigned int i1, i2;

    std::set<unsigned int> found_keys;   // For check_size_only

#ifdef AMBIT_USE_MPI
    int count = 0;  // Processor index

    if(!check_size_only)
    {
        int new_keys_per_processor = CalculateOneElectronIntegrals(orbital_map_1, orbital_map_2, true);
        new_keys_per_processor = (new_keys_per_processor + NumProcessors - 1)/NumProcessors;

        new_keys.clear();
        new_keys.reserve(new_keys_per_processor);
        new_values.clear();
        new_values.reserve(new_keys_per_processor);

        for(const auto& pair: integrals)
            found_keys.insert(pair.first);
    }
#endif

    if(!check_size_only)
    {
        if(include_core || include_core_subtraction)
            core_PT->UpdateIntegrals();
        if(include_valence_subtraction)
            valence_PT->UpdateIntegrals();
    }

    auto it_1 = orbital_map_1->begin();
    while(it_1 != orbital_map_1->end())
    {
        i1 = orbitals->state_index.at(it_1->first);
        const OrbitalInfo& s1 = it_1->first;

        auto it_2 = orbital_map_2->begin();
        while(it_2 != orbital_map_2->end())
        {
            if(it_1->first.Kappa() == it_2->first.Kappa())
            {
                i2 = orbitals->state_index.at(it_2->first);
                const OrbitalInfo& s2 = it_2->first;

                unsigned int key = GetKey(i1, i2);

                if(check_size_only)
                {
                    if(include_core || include_core_subtraction || include_valence_subtraction)
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
                    if(integrals.find(key) == integrals.end())
                    {
                    #endif
                        double value = 0.;
                        if(include_core)
                            value += core_PT->GetOneElectronDiagrams(s1, s2);
                        if(include_core_subtraction)
                            value += core_PT->GetOneElectronSubtraction(s1, s2);
                        if(include_valence_subtraction)
                            value += valence_PT->GetOneElectronSubtraction(s1, s2);

                    #ifdef AMBIT_USE_MPI
                        new_keys.push_back(key);
                        new_values.push_back(value);
                        }

                        found_keys.insert(key);
                        count++;
                        if(count == NumProcessors)
                            count = 0;
                    #else
                        integrals[key] = value;
                    #endif
                    }
                }
            }
            it_2++;
        }
        it_1++;
    }

    if(check_size_only)
        return found_keys.size();
    else
    #ifdef AMBIT_USE_MPI
    {   // Gather to root node, write to file, and read back in
        Write(write_file);
        clear();
        new_keys.clear();
        new_values.clear();
        Read(write_file);

        return integrals.size();
    }
    #else
    {   Write(write_file);
        return integrals.size();
    }
    #endif
}

#ifdef AMBIT_USE_MPI
void OneElectronMBPT::Write(const std::string& filename) const
{
    if(ProcessorRank == 0)
    {
        FILE* fp = fopen(filename.c_str(), "wb");

        // Write state index
        WriteOrbitalIndexes(orbitals->state_index, fp);

        // Get number of keys from all processes
        std::vector<unsigned int> num_integrals(NumProcessors);
        num_integrals[0] = new_keys.size();

        MPI_Status status;
        for(int proc = 1; proc < NumProcessors; proc++)
        {
            MPI_Recv(&num_integrals[proc], 1, MPI_UNSIGNED, proc, 1, MPI_COMM_WORLD, &status);
        }

        unsigned int total_integrals = std::accumulate(num_integrals.begin(), num_integrals.end(), size());
        fwrite(&total_integrals, sizeof(unsigned int), 1, fp);

        // Write root integrals
        for(auto& pair: integrals)
        {
            fwrite(&pair.first, sizeof(unsigned int), 1, fp);
            fwrite(&pair.second, sizeof(double), 1, fp);
        }

        for(int i = 0; i < num_integrals[0]; i++)
        {
            fwrite(&new_keys[i], sizeof(unsigned int), 1, fp);
            fwrite(&new_values[i], sizeof(double), 1, fp);
        }

        // Receive and write data from other processes
        std::vector<unsigned int> keys;
        std::vector<double> values;
        for(int proc = 1; proc < NumProcessors; proc++)
        {
            keys.resize(num_integrals[proc]);
            values.resize(num_integrals[proc]);

            MPI_Recv(keys.data(), num_integrals[proc], MPI_UNSIGNED, proc, 2, MPI_COMM_WORLD, &status);
            MPI_Recv(values.data(), num_integrals[proc], MPI_DOUBLE, proc, 3, MPI_COMM_WORLD, &status);

            for(int i = 0; i < num_integrals[proc]; i++)
            {
                fwrite(&keys[i], sizeof(unsigned int), 1, fp);
                fwrite(&values[i], sizeof(double), 1, fp);
            }
        }

        fclose(fp);
    }
    else
    {   // Send data to root. First number of keys.
        unsigned int num_keys = new_keys.size();
        MPI_Send(&num_keys, 1, MPI_UNSIGNED, 0, 1, MPI_COMM_WORLD);

        // Send keys and data to root
        MPI_Send(new_keys.data(), num_keys, MPI_UNSIGNED, 0, 2, MPI_COMM_WORLD);
        MPI_Send(new_values.data(), num_keys, MPI_DOUBLE, 0, 3, MPI_COMM_WORLD);
    }

    MPI_Barrier(MPI_COMM_WORLD);
}
#endif

void OneElectronMBPT::IncludeCore(bool include_mbpt, bool include_subtraction)
{
    include_core = include_mbpt && core_PT;
    include_core_subtraction = include_subtraction && core_PT;
}

void OneElectronMBPT::IncludeValence(bool include_subtraction)
{
    include_valence_subtraction = include_subtraction && valence_PT;
}

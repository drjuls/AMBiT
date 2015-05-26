#include "Atom.h"
#include "Include.h"
#include "HartreeFock/HartreeFocker.h"
#include "Configuration/ManyBodyOperator.h"
#include "HartreeFock/ConfigurationParser.h"
#include "Configuration/ConfigGenerator.h"
#ifdef AMBIT_USE_MPI
#include <mpi.h>
#endif

void Atom::Autoionization(pLevelConst target)
{
    double ionization_energy = user_input("DR/IonizationEnergy", 0.0);
    double energy_limit = user_input("DR/EnergyLimit", 0.0);

    auto math = MathConstant::Instance();

    ConfigGenerator gen(orbitals, user_input);
    pOrbitalMap valence(new OrbitalMap(*orbitals->valence));

    // Select orbitals in target wavefunction
    std::set<OrbitalInfo> target_info_set;
    pOrbitalMap target_map(new OrbitalMap(lattice));
    for(const auto& rconfig: *target->GetRelativisticConfigList())
    {
        for(const auto& orb_pair: rconfig)
            target_info_set.insert(orb_pair.first);
    }
    for(const auto& info: target_info_set)
    {
        target_map->AddState(valence->GetState(info));
    }

    // Target information
    Symmetry target_symmetry(target->GetTwoJ(), target->GetParity());
    const std::vector<double>& target_eigenvector = target->GetEigenvector();

    // Create operator for construction of continuum field
    pHFOperator hf_continuum;
    std::string config = user_input("DR/ContinuumFieldConfiguration", "");
    if(config == "")
        hf_continuum = hf_open;
    else
    {   hf_continuum.reset(hf_open->Clone());
        pCore continuum_core(new Core(lattice, config));
        for(auto& orbital: *continuum_core)
        {
            pOrbital basis_state = orbitals->all->GetState(orbital.first);
            continuum_core->AddState(basis_state);
        }

        hf_continuum->SetCore(continuum_core);
    }

    *outstream << "\nAutoionization rates:"
               << "\n E(eV)   A(ns)   J   P   LevelID" << std::endl;

    double energy_unit_conversion = math->HartreeEnergyIneV();
    double rate_unit_conversion = 0.413e8;

    // Loop over all HamiltonianIDs
    for(auto& key: *levels)
    {
        Symmetry sym = key->GetSymmetry();

        LevelVector levelvec = levels->GetLevels(key);
        pRelativisticConfigListConst compound_configs = key->GetRelativisticConfigList();
        if(!compound_configs && levelvec.size())
            compound_configs = levelvec[0]->GetRelativisticConfigList();

        // Loop over all levels
        int level_index = 0;
        for(auto level_it = levelvec.begin(); level_it != levelvec.end(); level_it++)
        {
            const auto& compound_eigenvector = (*level_it)->GetEigenvector();

            // Build continuum
            double eps_energy = (*level_it)->GetEnergy() - ionization_energy;
            if(eps_energy <= 0.)
            {
                level_index++;
                continue;
            }
            else if(energy_limit && (eps_energy > energy_limit))
                break;

            double rate = 0.0;

            // Get continuum wave (eps)
            Parity eps_parity = sym.GetParity() * target_symmetry.GetParity();
            for(int eps_twoJ = abs(sym.GetTwoJ() - target_symmetry.GetTwoJ()); eps_twoJ <= sym.GetTwoJ() + target_symmetry.GetTwoJ(); eps_twoJ += 2)
            {
                double partial = 0.0;

                // Kappa = (-1)^(j+1/2 + l) (j + 1/2) = (-1)^(j+1/2).P.(j+1/2)
                int eps_kappa = (eps_twoJ + 1)/2;
                eps_kappa = math->minus_one_to_the_power(eps_kappa) * Sign(eps_parity) * eps_kappa;

                pODESolver ode_solver(new AdamsSolver(hf->GetOPIntegrator()));
                HartreeFocker HF_Solver(ode_solver);
                pContinuumWave eps(new ContinuumWave(eps_kappa, 100, eps_energy));
                HF_Solver.CalculateContinuumWave(eps, hf_continuum);

                // Create two-body integrals
                pOrbitalMap continuum_map(new OrbitalMap(lattice));
                continuum_map->AddState(eps);

                pOrbitalManager all_orbitals(new OrbitalManager(*orbitals));
                all_orbitals->all->AddState(eps);
                all_orbitals->MakeStateIndexes();
                pOneElectronIntegrals one_body_integrals(new OneElectronIntegrals(all_orbitals, hf));
                pSlaterIntegrals two_body_integrals(new SlaterIntegralsDenseHash(all_orbitals, hartreeY));
                one_body_integrals->clear();
                one_body_integrals->CalculateOneElectronIntegrals(valence, continuum_map);
                two_body_integrals->clear();
                two_body_integrals->CalculateTwoElectronIntegrals(continuum_map, target_map, valence, valence);

                // Create operator for autoionization
                pTwoElectronCoulombOperator two_body_operator(new TwoElectronCoulombOperator<pSlaterIntegrals>(two_body_integrals));
                ManyBodyOperator<pOneElectronIntegrals, pTwoElectronCoulombOperator> H(one_body_integrals, two_body_operator);

                pRelativisticConfigList target_configs(new RelativisticConfigList(*target->GetRelativisticConfigList()));

                // Couple target and epsilon; start with maximum target M
                for(int target_twoM = target_symmetry.GetTwoJ(); target_twoM >= -target_symmetry.GetTwoJ(); target_twoM -= 2)
                {
                    int eps_twoM = sym.GetTwoJ() - target_twoM;
                    if(eps_twoM < -eps_twoJ)
                        continue;
                    else if (eps_twoM > eps_twoJ)
                        break;
                    ElectronInfo eps_info(100, eps_kappa, eps_twoM);

                    gen.GenerateProjections(target_configs, target_symmetry, target_twoM, angular_library);

                    // Coupling constant for |(eps_jlm; target_jlm) JJ >
                    double coupling = math->Wigner3j(eps_twoJ, target_symmetry.GetTwoJ(), sym.GetTwoJ(), eps_twoM, target_twoM);
                    coupling *= sqrt(double(sym.GetTwoJ() + 1)) * math->minus_one_to_the_power((eps_twoJ - target_symmetry.GetTwoJ() + eps_twoM + target_twoM)/2);

                    // Iterate over (generally larger) compound states first, for better multiporocessor load balance
                    double local_partial_j = 0.;
                    auto proj_jt = compound_configs->projection_begin();
                    int proj_index = 0;
                    while(proj_jt != compound_configs->projection_end())
                    {
                        if(proj_index%NumProcessors == ProcessorRank)
                        {
                            auto proj_it = target_configs->projection_begin();
                            while(proj_it != target_configs->projection_end())
                            {
                                double matrix_element = H.GetMatrixElement(*proj_it, *proj_jt, &eps_info);

                                if(matrix_element)
                                {
                                    // Summation over CSFs
                                    double coeff = 0.;

                                    for(auto coeff_i = proj_it.CSF_begin(); coeff_i != proj_it.CSF_end(); coeff_i++)
                                    {
                                        for(auto coeff_j = proj_jt.CSF_begin(); coeff_j != proj_jt.CSF_end(); coeff_j++)
                                        {
                                            coeff += (*coeff_i) * (*coeff_j)
                                                    * target_eigenvector[coeff_i.index()]
                                                    * compound_eigenvector[coeff_j.index()];
                                        }
                                    }

                                    local_partial_j += coupling * coeff * matrix_element;
                                }
                                proj_it++;
                            }
                        }
                        proj_jt++;
                        proj_index++;
                    }
                    #ifdef AMBIT_USE_MPI
                        double reduced_partial_j;
                        MPI_Allreduce(&local_partial_j, &reduced_partial_j, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                        partial += reduced_partial_j;
                    #else
                        partial += local_partial_j;
                    #endif
                }

                rate += 2. * math->Pi() * partial * partial;
            }

            char sep = ' ';
            *outstream << eps_energy * energy_unit_conversion << sep
                       << rate * rate_unit_conversion << sep
                       << sym.GetJ() <<  sep << Sign(sym.GetParity()) << sep
                       << key->Name() << ':' << level_index << std::endl;

            level_index++;
        }
    }
}

void Atom::AutoionizationEnergyGrid(pLevelConst target)
{
    double ionization_energy = user_input("DR/IonizationEnergy", 0.0);
    double grid_max = user_input("DR/EnergyLimit", 3.675);      // Default 100eV
    double grid_step = user_input("DR/EnergyGridStep", 0.03675); // Default 1eV
    double grid_min = user_input("DR/EnergyGridMin", 0.003675); // Default 0.1eV
    std::vector<double> energy_grid;
    int max_continuum_l = user_input("DR/ContinuumLMax", 6);    // Default maximum L = 6

    int pqn_offset = 100;
    energy_grid.reserve((grid_max - grid_min)/grid_step + 1);

    auto math = MathConstant::Instance();

    *outstream << "\nAutoionization rates:"
               << "\n E(eV)   A(ns)   J   P   LevelID" << std::endl;

    double energy_unit_conversion = math->HartreeEnergyIneV();
    double rate_unit_conversion = 0.413e8;

    ConfigGenerator gen(orbitals, user_input);
    pOrbitalMap valence(new OrbitalMap(*orbitals->valence));

    // Select orbitals in target wavefunction
    std::set<OrbitalInfo> target_info_set;
    pOrbitalMap target_map(new OrbitalMap(lattice));
    for(const auto& rconfig: *target->GetRelativisticConfigList())
    {
        for(const auto& orb_pair: rconfig)
            target_info_set.insert(orb_pair.first);
    }
    for(const auto& info: target_info_set)
    {
        target_map->AddState(valence->GetState(info));
    }

    // Target information
    Symmetry target_symmetry(target->GetTwoJ(), target->GetParity());
    const std::vector<double>& target_eigenvector = target->GetEigenvector();

    // Create operator for construction of continuum field
    pHFOperator hf_continuum;
    std::string config = user_input("DR/ContinuumFieldConfiguration", "");
    if(config == "")
        hf_continuum = hf_open;
    else
    {   hf_continuum.reset(hf_open->Clone());
        pCore continuum_core(new Core(lattice, config));
        for(auto& orbital: *continuum_core)
        {
            pOrbital basis_state = orbitals->all->GetState(orbital.first);
            continuum_core->AddState(basis_state);
        }

        hf_continuum->SetCore(continuum_core);
    }

    // Create continuum waves
    pOrbitalManager all_orbitals(new OrbitalManager(*orbitals));
    pOrbitalMap continuum_map(new OrbitalMap(lattice));
    pODESolver ode_solver(new AdamsSolver(hf->GetOPIntegrator()));
    HartreeFocker HF_Solver(ode_solver);

    int pqn = pqn_offset;
    for(double eps_energy = grid_min; eps_energy <= grid_max; eps_energy += grid_step)
    {
        for(int eps_l = 0; eps_l <= max_continuum_l; eps_l++)
        {
            for(int eps_kappa = -eps_l -1; eps_kappa <= eps_l; eps_kappa += 2 * mmax(eps_l, 1) + 1)
            {
                pContinuumWave eps(new ContinuumWave(eps_kappa, pqn, eps_energy));
                HF_Solver.CalculateContinuumWave(eps, hf_continuum);

                continuum_map->AddState(eps);
                all_orbitals->all->AddState(eps);
            }
        }

        energy_grid.push_back(eps_energy);
        pqn++;
    }

    // Create two-body integrals
    all_orbitals->MakeStateIndexes();
    pOneElectronIntegrals one_body_integrals(new OneElectronIntegrals(all_orbitals, hf));
    pSlaterIntegrals two_body_integrals(new SlaterIntegralsDenseHash(all_orbitals, hartreeY));
    one_body_integrals->clear();
    one_body_integrals->CalculateOneElectronIntegrals(valence, continuum_map);
    two_body_integrals->clear();
    two_body_integrals->CalculateTwoElectronIntegrals(continuum_map, target_map, valence, valence);

    // Create operator for autoionization
    pTwoElectronCoulombOperator two_body_operator(new TwoElectronCoulombOperator<pSlaterIntegrals>(two_body_integrals));
    ManyBodyOperator<pOneElectronIntegrals, pTwoElectronCoulombOperator> H(one_body_integrals, two_body_operator);

    // Create copies of target_configs with different two_M
    std::vector<pRelativisticConfigList> target_configs;
    target_configs.reserve(2 * target_symmetry.GetTwoJ() + 1);
    for(int two_m = target_symmetry.GetTwoJ(); two_m >= -target_symmetry.GetTwoJ(); two_m -= 2)
    {
        pRelativisticConfigList ptarget_config = std::make_shared<RelativisticConfigList>(*target->GetRelativisticConfigList());
        target_configs.push_back(ptarget_config);
        gen.GenerateProjections(ptarget_config, target_symmetry, two_m, angular_library);
    }

    // Loop over all HamiltonianIDs
    for(auto& key: *levels)
    {
        Symmetry sym = key->GetSymmetry();

        LevelVector levelvec = levels->GetLevels(key);
        pRelativisticConfigListConst compound_configs = key->GetRelativisticConfigList();
        if(!compound_configs && levelvec.size())
            compound_configs = levelvec[0]->GetRelativisticConfigList();

        // Loop over all levels
        int level_index = 0;
        for(auto level_it = levelvec.begin(); level_it != levelvec.end(); level_it++)
        {
            // Choose continuum energy
            double eps_energy = (*level_it)->GetEnergy() - ionization_energy;
            if(eps_energy > grid_max)
                break;

            auto energy_grid_it = std::lower_bound(energy_grid.begin(), energy_grid.end(), eps_energy);

            if(energy_grid_it == energy_grid.begin())
                pqn = pqn_offset;
            else if(energy_grid_it == energy_grid.end())
                pqn = pqn_offset + energy_grid.size() - 1;
            else
            {   pqn = energy_grid_it - energy_grid.begin() + pqn_offset;
                // Check it is closest
                if(*energy_grid_it - eps_energy > eps_energy - *(energy_grid_it-1))
                    pqn--;
            }

            const auto& compound_eigenvector = (*level_it)->GetEigenvector();

            double rate = 0.0;

            // Get continuum wave (eps)
            Parity eps_parity = sym.GetParity() * target_symmetry.GetParity();
            int max_eps_twoJ = mmin(sym.GetTwoJ() + target_symmetry.GetTwoJ(), 2 * max_continuum_l + 1);
            for(int eps_twoJ = abs(sym.GetTwoJ() - target_symmetry.GetTwoJ()); eps_twoJ <= max_eps_twoJ; eps_twoJ += 2)
            {
                double partial = 0.0;

                // Kappa = (-1)^(j+1/2 + l) (j + 1/2) = (-1)^(j+1/2).P.(j+1/2)
                int eps_kappa = (eps_twoJ + 1)/2;
                eps_kappa = math->minus_one_to_the_power(eps_kappa) * Sign(eps_parity) * eps_kappa;

                // Couple target and epsilon; start with maximum target M
                auto target_config_it = target_configs.begin();
                for(int target_twoM = target_symmetry.GetTwoJ(); target_twoM >= -target_symmetry.GetTwoJ(); target_twoM -= 2)
                {
                    pRelativisticConfigList& ptarget_config = *target_config_it;
                    int eps_twoM = sym.GetTwoJ() - target_twoM;
                    if(eps_twoM < -eps_twoJ)
                        continue;
                    else if (eps_twoM > eps_twoJ)
                        break;
                    ElectronInfo eps_info(pqn, eps_kappa, eps_twoM);

                    // Coupling constant for |(eps_jlm; target_jlm) JJ >
                    double coupling = math->Wigner3j(eps_twoJ, target_symmetry.GetTwoJ(), sym.GetTwoJ(), eps_twoM, target_twoM);
                    coupling *= sqrt(double(sym.GetTwoJ() + 1)) * math->minus_one_to_the_power((eps_twoJ - target_symmetry.GetTwoJ() + eps_twoM + target_twoM)/2);

                    // Iterate over (generally larger) compound states first, for better multiprocessor load balance
                    double local_partial_j = 0.;
                    auto proj_jt = compound_configs->projection_begin();
                    int proj_index = 0;
                    while(proj_jt != compound_configs->projection_end())
                    {
                        if(proj_index%NumProcessors == ProcessorRank)
                        {
                            auto proj_it = ptarget_config->projection_begin();
                            while(proj_it != ptarget_config->projection_end())
                            {
                                double matrix_element = H.GetMatrixElement(*proj_it, *proj_jt, &eps_info);

                                if(matrix_element)
                                {
                                    // Summation over CSFs
                                    double coeff = 0.;

                                    for(auto coeff_i = proj_it.CSF_begin(); coeff_i != proj_it.CSF_end(); coeff_i++)
                                    {
                                        for(auto coeff_j = proj_jt.CSF_begin(); coeff_j != proj_jt.CSF_end(); coeff_j++)
                                        {
                                            coeff += (*coeff_i) * (*coeff_j)
                                            * target_eigenvector[coeff_i.index()]
                                            * compound_eigenvector[coeff_j.index()];
                                        }
                                    }

                                    local_partial_j += coupling * coeff * matrix_element;
                                }
                                proj_it++;
                            }
                        }
                        proj_jt++;
                        proj_index++;
                    }
                    #ifdef AMBIT_USE_MPI
                        double reduced_partial_j;
                        MPI_Allreduce(&local_partial_j, &reduced_partial_j, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                        partial += reduced_partial_j;
                    #else
                        partial += local_partial_j;
                    #endif

                    target_config_it++;
                }

                rate += 2. * math->Pi() * partial * partial;
            }

            char sep = ' ';
            *outstream << eps_energy * energy_unit_conversion << sep
                       << rate * rate_unit_conversion << sep
                       << sym.GetJ() <<  sep << Sign(sym.GetParity()) << sep
                       << key->Name() << ':' << level_index << std::endl;

            level_index++;
        }
    }
}

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
    double energy_limit = user_input("DR/EnergyLimit", -1.0);   // Default no limit
    int min_continuum_l = user_input("DR/ContinuumLMin", 0);
    int max_continuum_l = user_input("DR/ContinuumLMax", 6);    // Default maximum L = 6

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
    pLevel target_copy = std::make_shared<Level>(*target);

    // Create copies of target_configs with different two_M
    std::vector<pRelativisticConfigList> target_configs;
    target_configs.reserve(2 * target_symmetry.GetTwoJ() + 1);
    for(int two_m = target_symmetry.GetTwoJ(); two_m >= -target_symmetry.GetTwoJ(); two_m -= 2)
    {
        pRelativisticConfigList ptarget_config = std::make_shared<RelativisticConfigList>(*target->GetRelativisticConfigList());
        target_configs.push_back(ptarget_config);
        gen.GenerateProjections(ptarget_config, target_symmetry, two_m, angular_library);
    }

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

        // Loop over all levels
        int level_index = 0;
        for(auto level_it = levelvec.begin(); level_it != levelvec.end(); level_it++)
        {
            // Build continuum
            double eps_energy = (*level_it)->GetEnergy() - ionization_energy;
            if(eps_energy <= 0.)
            {
                level_index++;
                continue;
            }
            else if(energy_limit > 0.0 && (eps_energy > energy_limit))
                break;

            double rate = 0.0;

            // Get continuum wave (eps)
            Parity eps_parity = sym.GetParity() * target_symmetry.GetParity();

            // Minimum L and Parity -> Minimum J
            int min_eps_twoJ;
            if(eps_parity == (min_continuum_l%2? Parity::odd: Parity::even))
                min_eps_twoJ = mmax(2 * min_continuum_l - 1, 1);
            else
                min_eps_twoJ = 2 * min_continuum_l + 1;

            // continuum J should not be less than (sym.J - target.J)
            min_eps_twoJ = mmax(min_eps_twoJ, abs(sym.GetTwoJ() - target_symmetry.GetTwoJ()));

            // Maximum L and Parity -> Maximum J
            int max_eps_twoJ;
            if(eps_parity == (max_continuum_l%2? Parity::odd: Parity::even))
                max_eps_twoJ = 2 * max_continuum_l + 1;
            else
                max_eps_twoJ = 2 * max_continuum_l - 1;

            // continuum J should not exceed sym.J + target.J
            max_eps_twoJ = mmin(max_eps_twoJ, sym.GetTwoJ() + target_symmetry.GetTwoJ());

            for(int eps_twoJ = min_eps_twoJ; eps_twoJ <= max_eps_twoJ; eps_twoJ += 2)
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

                // Couple target and epsilon; start with maximum target M
                auto target_config_it = target_configs.begin();
                for(int target_twoM = target_symmetry.GetTwoJ(); target_twoM >= -target_symmetry.GetTwoJ(); target_twoM -= 2)
                {
                    int eps_twoM = sym.GetTwoJ() - target_twoM;
                    if(eps_twoM < -eps_twoJ)
                        continue;
                    else if (eps_twoM > eps_twoJ)
                        break;
                    ElectronInfo eps_info(100, eps_kappa, eps_twoM);

                    target_copy->SetRelativisticConfigList(*target_config_it);

                    // Coupling constant for |(eps_jlm; target_jlm) JJ >
                    double coupling = math->Wigner3j(eps_twoJ, target_symmetry.GetTwoJ(), sym.GetTwoJ(), eps_twoM, target_twoM);
                    coupling *= sqrt(double(sym.GetTwoJ() + 1)) * math->minus_one_to_the_power((eps_twoJ - target_symmetry.GetTwoJ() + eps_twoM + target_twoM)/2);

                    partial += coupling * H.GetMatrixElement(**level_it, *target_copy, &eps_info);
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

void Atom::AutoionizationEnergyGrid(pLevelConst target)
{
    double ionization_energy = user_input("DR/IonizationEnergy", 0.0);
    double energy_limit = user_input("DR/EnergyLimit", -1.0);     // Default no limit
    double grid_min = user_input("DR/EnergyGrid/Min", 0.003675);  // Default 0.1eV
    double grid_max = user_input("DR/EnergyGrid/Max", 3.675);     // Default 100eV
    double grid_step = user_input("DR/EnergyGrid/Step", 0.03675); // Default 1eV
    std::vector<double> energy_grid;
    int min_continuum_l = user_input("DR/ContinuumLMin", 0);
    int max_continuum_l = user_input("DR/ContinuumLMax", 6);    // Default maximum L = 6

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
    LevelVector target_levelvector;
    target_levelvector.push_back(std::make_shared<Level>(*target));

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

    int pqn_offset = 100;
    if(energy_limit > 0.0 && energy_limit < grid_max)
        grid_max = energy_limit;
    energy_grid.reserve((grid_max - grid_min)/grid_step + 1);

    int pqn = pqn_offset;
    for(double eps_energy = grid_min; eps_energy <= grid_max; eps_energy += grid_step)
    {
        for(int eps_l = min_continuum_l; eps_l <= max_continuum_l; eps_l++)
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

        // Group levels that use the same continuum energy
        auto level_it = levelvec.begin();
        unsigned int level_index = 0;
        for(int energy_grid_point = 0; energy_grid_point < energy_grid.size(); energy_grid_point++)
        {
            if(level_it == levelvec.end() ||
               (energy_limit > 0.0 && (*level_it)->GetEnergy() > energy_limit + ionization_energy))
                break;
            pqn = pqn_offset + energy_grid_point;

            // Get levels in the energy range
            auto level_end_of_energy_range = level_it;
            if((energy_grid_point == energy_grid.size() - 1) && (energy_limit <= 0.0)) // Last point and no energy limit
                level_end_of_energy_range = levelvec.end();
            else
            {   double max_energy;
                if(energy_grid_point == energy_grid.size() - 1)
                    max_energy = energy_limit + ionization_energy;
                else
                    max_energy = 0.5 * (energy_grid[energy_grid_point] + energy_grid[energy_grid_point+1]) + ionization_energy;

                while(level_end_of_energy_range != levelvec.end()
                      && (*level_end_of_energy_range)->GetEnergy() < max_energy)
                {
                    level_end_of_energy_range++;
                }
            }

            unsigned int num_levels = level_end_of_energy_range - level_it;
            if(num_levels != 0)
            {
                LevelVector level_subset(level_it, level_end_of_energy_range);
                unsigned int solution;
                std::vector<double> rate(num_levels, 0.0);

                // Get continuum wave (eps)
                Parity eps_parity = sym.GetParity() * target_symmetry.GetParity();

                // Minimum L and Parity -> Minimum J
                int min_eps_twoJ;
                if(eps_parity == (min_continuum_l%2? Parity::odd: Parity::even))
                    min_eps_twoJ = mmax(2 * min_continuum_l - 1, 1);
                else
                    min_eps_twoJ = 2 * min_continuum_l + 1;

                // continuum J should not be less than (sym.J - target.J)
                min_eps_twoJ = mmax(min_eps_twoJ, abs(sym.GetTwoJ() - target_symmetry.GetTwoJ()));
                
                // Maximum L and Parity -> Maximum J
                int max_eps_twoJ;
                if(eps_parity == (max_continuum_l%2? Parity::odd: Parity::even))
                    max_eps_twoJ = 2 * max_continuum_l + 1;
                else
                    max_eps_twoJ = 2 * max_continuum_l - 1;

                // continuum J should not exceed sym.J + target.J
                max_eps_twoJ = mmin(max_eps_twoJ, sym.GetTwoJ() + target_symmetry.GetTwoJ());

                for(int eps_twoJ = min_eps_twoJ; eps_twoJ <= max_eps_twoJ; eps_twoJ += 2)
                {
                    std::vector<double> partial(num_levels, 0.0);
                    std::vector<double> coeff(num_levels);

                    // Kappa = (-1)^(j+1/2 + l) (j + 1/2) = (-1)^(j+1/2).P.(j+1/2)
                    int eps_kappa = (eps_twoJ + 1)/2;
                    eps_kappa = math->minus_one_to_the_power(eps_kappa) * Sign(eps_parity) * eps_kappa;

                    // Couple target and epsilon; start with maximum target M
                    auto target_config_it = target_configs.begin();
                    for(int target_twoM = target_symmetry.GetTwoJ(); target_twoM >= -target_symmetry.GetTwoJ(); target_twoM -= 2)
                    {
                        int eps_twoM = sym.GetTwoJ() - target_twoM;
                        if(eps_twoM < -eps_twoJ)
                            continue;
                        else if (eps_twoM > eps_twoJ)
                            break;
                        ElectronInfo eps_info(pqn, eps_kappa, eps_twoM);

                        target_levelvector[0]->SetRelativisticConfigList(*target_config_it);

                        // Coupling constant for |(eps_jlm; target_jlm) JJ >
                        double coupling = math->Wigner3j(eps_twoJ, target_symmetry.GetTwoJ(), sym.GetTwoJ(), eps_twoM, target_twoM);
                        coupling *= sqrt(double(sym.GetTwoJ() + 1)) * math->minus_one_to_the_power((eps_twoJ - target_symmetry.GetTwoJ() + eps_twoM + target_twoM)/2);

                        // Iterate over (generally larger) compound states first, for better multiprocessor load balance
                        std::vector<double> partial_contribution = H.GetMatrixElement(level_subset, target_levelvector, &eps_info);
                        for(solution = 0; solution < num_levels; solution++)
                            partial[solution] += coupling * partial_contribution[solution];

                        target_config_it++;
                    }

                    for(solution = 0; solution < num_levels; solution++)
                        rate[solution] += 2. * math->Pi() * partial[solution] * partial[solution];
                }

                if(ProcessorRank == 0)
                {
                    char sep = ' ';
                    for(solution = 0; solution < num_levels; solution++)
                    {
                        double eps_energy = level_subset[solution]->GetEnergy() - ionization_energy;
                        *outstream << std::setprecision(5);
                        *outstream << eps_energy * energy_unit_conversion << sep
                                   << rate[solution] * rate_unit_conversion << sep
                                   << sym.GetJ() <<  sep << Sign(sym.GetParity()) << sep
                                   << key->Name() << ':' << level_index + solution << std::endl;
                    }
                }
            }

            level_it = level_end_of_energy_range;
            level_index += num_levels;
        }
    }
}

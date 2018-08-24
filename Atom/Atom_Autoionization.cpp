#include "Atom.h"
#include "Include.h"
#include "HartreeFock/HartreeFocker.h"
#include "Configuration/ManyBodyOperator.h"
#include "HartreeFock/ConfigurationParser.h"
#include "Configuration/ConfigGenerator.h"
#ifdef AMBIT_USE_MPI
#include <mpi.h>
#endif
#include "ExternalField/Hyperfine.h"

namespace Ambit
{
void Atom::InternalConversion(const LevelVector& source)
{
    auto math = MathConstant::Instance();
    double nuclear_energy = user_input("IC/NuclearEnergyEV", -1.0)/math->HartreeEnergyIneV();
    int min_continuum_l = user_input("IC/ContinuumLMin", 0);
    int max_continuum_l = user_input("IC/ContinuumLMax", 6);    // Default maximum L = 6

    if(nuclear_energy <= 0.0)
    {   *outstream << "Required to set IC/NuclearEnergyEV to a positive value." << std::endl;
        exit(1);
    }

    // Add nuclear_energy to source level energy
    nuclear_energy += source.levels[0]->GetEnergy();

    user_input.set_prefix(std::string("IC"));
    GeneralisedHyperfineCalculator calculon(user_input, *this);
    user_input.set_prefix("");
    auto op = calculon.GetOperator();

    pODESolver ode_solver = std::make_shared<AdamsSolver>(hf_open->GetIntegrator());
    HartreeFocker HF_Solver(ode_solver);

    char sep = ' ';
    *outstream << std::setprecision(5);

    ConfigGenerator gen(orbitals, user_input);

    // Source information
    Symmetry source_symmetry(source.hID->GetSymmetry());

    // Create operator for construction of continuum field
    pHFOperator hf_continuum;
    std::string config = user_input("IC/Residue", "");
    if(config == "")
        hf_continuum = hf_open;
    else
    {   hf_continuum = hf_open->Clone();
        pCore continuum_core(new Core(lattice, config));
        for(auto& orbital: *continuum_core)
        {
            pOrbital basis_state = orbitals->all->GetState(orbital.first);
            continuum_core->AddState(basis_state);
        }

        hf_continuum->SetCore(continuum_core);
    }

    *outstream << "\nInternal conversion rates/B(M->G):"
               << "\n E(eV)   IC (a.u.)   J   P   LevelID" << std::endl;

    double energy_unit_conversion = math->HartreeEnergyIneV();
    double total_rate = 0.0;

    // Loop over all HamiltonianIDs
    for(auto& key: levels->keys)
    {
        Symmetry sym = key->GetSymmetry();
        auto levelvec = levels->GetLevels(key);

        // One final level at a time. Create a copy of the final level for summation over M.
        LevelVector final_copy(levelvec.configs, levelvec.levels[0]);
        final_copy.hID = key;

        // Create copies of final_configs with different two_M
        std::vector<pRelativisticConfigList> final_configs;
        final_configs.reserve(2 * sym.GetTwoJ() + 1);
        for(int two_m = sym.GetTwoJ(); two_m >= -sym.GetTwoJ(); two_m -= 2)
        {
            pRelativisticConfigList pfinal_config = std::make_shared<RelativisticConfigList>(*levelvec.configs);
            final_configs.push_back(pfinal_config);
            gen.GenerateProjections(pfinal_config, sym, two_m, angular_library);
        }

        // Loop over all levels
        int level_index = 0;
        for(auto level_it = levelvec.levels.begin(); level_it != levelvec.levels.end(); level_it++)
        {
            // Create final
            final_copy.levels[0] = *level_it;

            // Build continuum
            double eps_energy = nuclear_energy - (*level_it)->GetEnergy();
            if(eps_energy <= 1.e-4)
                break;

            double rate = 0.0;

            // Get continuum wave (eps)
            Parity eps_parity = sym.GetParity() * source_symmetry.GetParity() * op->GetParity();

            // Minimum L and Parity -> Minimum J
            int min_eps_twoJ;
            if(eps_parity == (min_continuum_l%2? Parity::odd: Parity::even))
                min_eps_twoJ = mmax(2 * min_continuum_l - 1, 1);
            else
                min_eps_twoJ = 2 * min_continuum_l + 1;

            // continuum J should not be less than |J_i - K| - j_nu
            min_eps_twoJ = mmax(min_eps_twoJ, abs(source_symmetry.GetTwoJ() - 2*op->GetK()) - sym.GetTwoJ());

            // Maximum L and Parity -> Maximum J
            int max_eps_twoJ;
            if(eps_parity == (max_continuum_l%2? Parity::odd: Parity::even))
                max_eps_twoJ = 2 * max_continuum_l + 1;
            else
                max_eps_twoJ = 2 * max_continuum_l - 1;

            // continuum J should not exceed J_i + K + J_nu
            max_eps_twoJ = mmin(max_eps_twoJ, source_symmetry.GetTwoJ() + 2*op->GetK() + sym.GetTwoJ());

            for(int eps_twoJ = min_eps_twoJ; eps_twoJ <= max_eps_twoJ; eps_twoJ += 2)
            {
                double partial = 0.0;

                // Kappa = (-1)^(j+1/2 + l) (j + 1/2) = (-1)^(j+1/2).P.(j+1/2)
                int eps_kappa = (eps_twoJ + 1)/2;
                eps_kappa = math->minus_one_to_the_power(eps_kappa) * Sign(eps_parity) * eps_kappa;

                pContinuumWave eps = std::make_shared<ContinuumWave>(eps_kappa, 100, eps_energy);
                HF_Solver.CalculateContinuumWave(eps, hf_continuum);

                // Create operator integrals
                pOrbitalMap continuum_map = std::make_shared<OrbitalMap>(lattice);
                continuum_map->AddState(eps);
                pOrbitalManager all_orbitals(new OrbitalManager(*orbitals));
                all_orbitals->all->AddState(eps);
                all_orbitals->MakeStateIndexes();
                pTransitionIntegrals one_body_integrals = std::make_shared<TransitionIntegrals>(all_orbitals, op);
                one_body_integrals->clear();
                one_body_integrals->CalculateOneElectronIntegrals(orbitals->valence, continuum_map);

                // Create many body operator
                ManyBodyOperator<pTransitionIntegrals> T(one_body_integrals);

                // Sum over M_nu (final)
                auto final_config_it = final_configs.begin();
                for(int final_twoM = sym.GetTwoJ(); final_twoM >= -sym.GetTwoJ(); final_twoM -= 2)
                {
                    final_copy.configs = *final_config_it++;

                    // Sum over M_epsilon
                    int eps_twoM_min = mmax(source_symmetry.GetTwoJ() - final_twoM - 2 * op->GetK(), -eps_twoJ);
                    int eps_twoM_max = mmin(source_symmetry.GetTwoJ() - final_twoM + 2 * op->GetK(), eps_twoJ);

                    for(int eps_twoM = eps_twoM_min; eps_twoM <= eps_twoM_max; eps_twoM += 2)
                    {
                        ElectronInfo eps_info(100, eps_kappa, eps_twoM);
                        double matrixelement = T.GetMatrixElement(source, final_copy, &eps_info)[0];
                        partial += gsl_pow_2(matrixelement);
                    }
                }

                *outstream << "    " << eps->Name() << sep << partial << "\n";
                rate += partial;
            }

            *outstream << eps_energy * energy_unit_conversion << sep
                       << rate << sep
                       << sym.GetJ() <<  sep << Sign(sym.GetParity()) << sep
                       << key->Name() << ':' << level_index << std::endl;

            total_rate += rate;
            level_index++;
        }
    }

    total_rate *= 8. * gsl_pow_2(math->Pi()/(2*op->GetK()+1));

    *outstream << "\nInternal conversion rate/B(M->G) (a.u.) = " << total_rate << std::endl;
}

void Atom::InternalConversionConfigurationAveraged(const LevelVector& source)
{
    auto math = MathConstant::Instance();
    double nuclear_energy = user_input("IC/NuclearEnergyEV", -1.0)/math->HartreeEnergyIneV();
    int min_continuum_l = user_input("IC/ContinuumLMin", 0);
    int max_continuum_l = user_input("IC/ContinuumLMax", 6);    // Default maximum L = 6

    if(nuclear_energy <= 0.0)
    {   *outstream << "Required to set IC/NuclearEnergyEV to a positive value." << std::endl;
        exit(1);
    }

    user_input.set_prefix(std::string("IC"));
    GeneralisedHyperfineCalculator calculon(user_input, *this);
    user_input.set_prefix("");
    auto op = calculon.GetOperator();

    pODESolver ode_solver = std::make_shared<AdamsSolver>(hf_open->GetIntegrator());
    HartreeFocker HF_Solver(ode_solver);

    char sep = ' ';
    *outstream << std::setprecision(5);

    // Get source level occupancies
    Configuration<OrbitalInfo, double> source_occupancies;

    auto& level = source.levels.front();
    const std::vector<double>& eigenvector = level->GetEigenvector();

    const double* eigenvector_csf = eigenvector.data();
    for(auto& rconfig: *source.configs)
    {
        double contribution = 0.0;
        for(unsigned int j = 0; j < rconfig.NumCSFs(); j++)
        {   contribution += gsl_pow_2(*eigenvector_csf);
            eigenvector_csf++;
        }

        // Convert rconfig occupancy to double
        Configuration<OrbitalInfo, double> double_config;
        double_config += rconfig;

        source_occupancies += double_config * contribution;
    }

    // Print source configuration
    *outstream << "IC Source: " << source_occupancies << std::endl;

    double total_rate = 0.0;
    for(auto& initialpair: source_occupancies)
    {
        // Create operator for construction of continuum field
        pHFOperator hf_continuum;
        std::string config = user_input("IC/Residue", "");
        if(config == "")
            hf_continuum = hf_open;
        else
        {   hf_continuum = hf_open->Clone();
            pCore continuum_core(new Core(lattice, config));
            for(auto& orbital: *continuum_core)
            {
                pOrbital basis_state = orbitals->all->GetState(orbital.first);
                continuum_core->AddState(basis_state);
            }

            hf_continuum->SetCore(continuum_core);
        }

        auto initial = orbitals->valence->GetState(initialpair.first);

        // Get continuum wave (eps)
        Parity eps_parity = op->GetParity() * initial->GetParity();

        // Minimum L and Parity -> Minimum J
        int min_eps_twoJ;
        if(eps_parity == (min_continuum_l%2? Parity::odd: Parity::even))
            min_eps_twoJ = mmax(2 * min_continuum_l - 1, 1);
        else
            min_eps_twoJ = 2 * min_continuum_l + 1;

        // continuum J should not be less than |J_i - K|
        min_eps_twoJ = mmax(min_eps_twoJ, abs(2 * op->GetK() - initial->TwoJ()));

        // Maximum L and Parity -> Maximum J
        int max_eps_twoJ;
        if(eps_parity == (max_continuum_l%2? Parity::odd: Parity::even))
            max_eps_twoJ = 2 * max_continuum_l + 1;
        else
            max_eps_twoJ = 2 * max_continuum_l - 1;

        // continuum J should not exceed J_i + K
        max_eps_twoJ = mmin(max_eps_twoJ, 2 * op->GetK() + initial->TwoJ());

        double eps_energy = nuclear_energy + initial->Energy();
        if(eps_energy <= 1.e-4)
            continue;

        double orbital_rate = 0.0;
        for(int eps_twoJ = min_eps_twoJ; eps_twoJ <= max_eps_twoJ; eps_twoJ += 2)
        {
            double partial = 0.0;

            // Kappa = (-1)^(j+1/2 + l) (j + 1/2) = (-1)^(j+1/2).P.(j+1/2)
            int eps_kappa = (eps_twoJ + 1)/2;
            eps_kappa = math->minus_one_to_the_power(eps_kappa) * Sign(eps_parity) * eps_kappa;

            pContinuumWave eps = std::make_shared<ContinuumWave>(eps_kappa, 100, eps_energy);
            HF_Solver.CalculateContinuumWave(eps, hf_continuum);

            partial = op->GetReducedMatrixElement(*initial, *eps);
            *outstream << "    " << initial->Name() << sep << eps->Name() << " => " << gsl_pow_2(partial) << '\n';

            orbital_rate += gsl_pow_2(partial) * initialpair.second/abs(2 * initial->Kappa());
        }

        *outstream << "  Total " << initial->Name() << " = " << orbital_rate << "\n";
        total_rate += orbital_rate;
    }

    total_rate *= 8. * gsl_pow_2(math->Pi()/(2*op->GetK()+1));

    *outstream << "\nInternal conversion rate/B(M->G) (a.u.) = " << total_rate << std::endl;
}

void Atom::Autoionization(const LevelVector& target)
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
    for(const auto& rconfig: *target.configs)
    {
        for(const auto& orb_pair: rconfig)
            target_info_set.insert(orb_pair.first);
    }
    for(const auto& info: target_info_set)
    {
        target_map->AddState(valence->GetState(info));
    }
    // Add holes from compound states
    target_map->AddStates(*orbitals->hole);

    // Target information
    Symmetry target_symmetry(target.hID->GetSymmetry());
    LevelVector target_copy(target.configs, target.levels[0]);

    // Create copies of target_configs with different two_M
    std::vector<pRelativisticConfigList> target_configs;
    target_configs.reserve(2 * target_symmetry.GetTwoJ() + 1);
    for(int two_m = target_symmetry.GetTwoJ(); two_m >= -target_symmetry.GetTwoJ(); two_m -= 2)
    {
        pRelativisticConfigList ptarget_config = std::make_shared<RelativisticConfigList>(*target.configs);
        target_configs.push_back(ptarget_config);
        gen.GenerateProjections(ptarget_config, target_symmetry, two_m, angular_library);
    }

    // Create operator for construction of continuum field
    pHFOperator hf_continuum;
    std::string config = user_input("DR/ContinuumResidue", "");
    if(config == "")
        hf_continuum = hf_open;
    else
    {   hf_continuum = hf_open->Clone();
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
    double rate_unit_conversion = math->AtomicFrequencySI() * 1.e-9;    // (ns);

    // Loop over all HamiltonianIDs
    for(auto& key: levels->keys)
    {
        Symmetry sym = key->GetSymmetry();
        auto levelvec = levels->GetLevels(key);

        // Compound target
        LevelVector levelvec_single(levelvec.configs, levelvec.levels[0]);
        levelvec_single.hID = key;

        // Loop over all levels
        int level_index = 0;
        for(auto level_it = levelvec.levels.begin(); level_it != levelvec.levels.end(); level_it++)
        {
            // Create target
            levelvec_single.levels[0] = *level_it;

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

                pODESolver ode_solver(new AdamsSolver(hf->GetIntegrator()));
                HartreeFocker HF_Solver(ode_solver);
                pContinuumWave eps(new ContinuumWave(eps_kappa, 100, eps_energy));
                HF_Solver.CalculateContinuumWave(eps, hf_continuum);

                // Create two-body integrals
                pOrbitalMap continuum_map(new OrbitalMap(lattice));
                continuum_map->AddState(eps);

                pOrbitalManager all_orbitals(new OrbitalManager(*orbitals));
                all_orbitals->all->AddState(eps);
                all_orbitals->MakeStateIndexes();
                pHFIntegrals one_body_integrals(new HFIntegrals(all_orbitals, hf));
                pSlaterIntegrals two_body_integrals(new SlaterIntegralsDenseHash(all_orbitals, hartreeY));
                one_body_integrals->clear();
                one_body_integrals->CalculateOneElectronIntegrals(continuum_map, valence);
                two_body_integrals->clear();
                two_body_integrals->CalculateTwoElectronIntegrals(continuum_map, target_map, valence, valence);

                // Create operator for autoionization
                pTwoElectronCoulombOperator two_body_operator = std::make_shared<TwoElectronCoulombOperator>(two_body_integrals);
                ManyBodyOperator<pHFIntegrals, pTwoElectronCoulombOperator> H(one_body_integrals, two_body_operator);

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

                    target_copy.configs = *target_config_it;

                    // Coupling constant for |(eps_jlm; target_jlm) JJ >
                    double coupling = math->Wigner3j(eps_twoJ, target_symmetry.GetTwoJ(), sym.GetTwoJ(), eps_twoM, target_twoM);
                    coupling *= sqrt(double(sym.GetTwoJ() + 1)) * math->minus_one_to_the_power((eps_twoJ - target_symmetry.GetTwoJ() + eps_twoM + target_twoM)/2);

                    partial += coupling * H.GetMatrixElement(levelvec_single, target_copy, &eps_info)[0];
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

void Atom::AutoionizationEnergyGrid(const LevelVector& target)
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
    double rate_unit_conversion = math->AtomicFrequencySI() * 1.e-9;    // (ns);

    ConfigGenerator gen(orbitals, user_input);
    pOrbitalMap valence(new OrbitalMap(*orbitals->valence));

    // Select orbitals in target wavefunction
    std::set<OrbitalInfo> target_info_set;
    pOrbitalMap target_map(new OrbitalMap(lattice));
    for(const auto& rconfig: *target.configs)
    {
        for(const auto& orb_pair: rconfig)
            target_info_set.insert(orb_pair.first);
    }
    for(const auto& info: target_info_set)
    {
        target_map->AddState(valence->GetState(info));
    }
    // Add holes from compound states
    target_map->AddStates(*orbitals->hole);

    // Target information
    Symmetry target_symmetry(target.hID->GetSymmetry());
    LevelVector target_levelvector(target);

    // Create operator for construction of continuum field
    pHFOperator hf_continuum;
    std::string config = user_input("DR/ContinuumResidue", "");
    if(config == "")
        hf_continuum = hf_open;
    else
    {   hf_continuum = hf_open->Clone();
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
    pODESolver ode_solver(new AdamsSolver(hf->GetIntegrator()));
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
    pHFIntegrals one_body_integrals(new HFIntegrals(all_orbitals, hf));
    pSlaterIntegrals two_body_integrals(new SlaterIntegralsDenseHash(all_orbitals, hartreeY));
    one_body_integrals->clear();
    one_body_integrals->CalculateOneElectronIntegrals(continuum_map, valence);
    two_body_integrals->clear();
    two_body_integrals->CalculateTwoElectronIntegrals(continuum_map, target_map, valence, valence);

    // Create operator for autoionization
    pTwoElectronCoulombOperator two_body_operator = std::make_shared<TwoElectronCoulombOperator>(two_body_integrals);
    ManyBodyOperator<pHFIntegrals, pTwoElectronCoulombOperator> H(one_body_integrals, two_body_operator);

    // Create copies of target_configs with different two_M
    std::vector<pRelativisticConfigList> target_configs;
    target_configs.reserve(2 * target_symmetry.GetTwoJ() + 1);
    for(int two_m = target_symmetry.GetTwoJ(); two_m >= -target_symmetry.GetTwoJ(); two_m -= 2)
    {
        pRelativisticConfigList ptarget_config = std::make_shared<RelativisticConfigList>(*target.configs);
        target_configs.push_back(ptarget_config);
        gen.GenerateProjections(ptarget_config, target_symmetry, two_m, angular_library);
    }

    // Loop over all HamiltonianIDs
    for(auto& key: levels->keys)
    {
        Symmetry sym = key->GetSymmetry();
        LevelVector levelvec = levels->GetLevels(key);

        // Group levels that use the same continuum energy
        auto level_it = levelvec.levels.begin();
        unsigned int level_index = 0;
        for(int energy_grid_point = 0; energy_grid_point < energy_grid.size(); energy_grid_point++)
        {
            if(level_it == levelvec.levels.end() ||
               (energy_limit > 0.0 && (*level_it)->GetEnergy() > energy_limit + ionization_energy))
                break;
            pqn = pqn_offset + energy_grid_point;

            // Get levels in the energy range
            auto level_end_of_energy_range = level_it;
            if((energy_grid_point == energy_grid.size() - 1) && (energy_limit <= 0.0)) // Last point and no energy limit
                level_end_of_energy_range = levelvec.levels.end();
            else
            {   double max_energy;
                if(energy_grid_point == energy_grid.size() - 1)
                    max_energy = energy_limit + ionization_energy;
                else
                    max_energy = 0.5 * (energy_grid[energy_grid_point] + energy_grid[energy_grid_point+1]) + ionization_energy;

                while(level_end_of_energy_range != levelvec.levels.end()
                      && (*level_end_of_energy_range)->GetEnergy() < max_energy)
                {
                    level_end_of_energy_range++;
                }
            }

            unsigned int num_levels = level_end_of_energy_range - level_it;
            if(num_levels != 0)
            {
                LevelVector level_subset;
                level_subset.hID = levelvec.hID;
                level_subset.configs = levelvec.configs;
                level_subset.levels.insert(level_subset.levels.begin(), level_it, level_end_of_energy_range);

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

                        target_levelvector.configs = *target_config_it;

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
                        double eps_energy = level_subset.levels[solution]->GetEnergy() - ionization_energy;
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

void Atom::AutoionizationConfigurationAveraged(const LevelVector& target)
{
    // Get fractional occupations of target orbitals
    OccupationMap target_occ;

    const std::vector<double>& eigenvector = target.levels[0]->GetEigenvector();
    const double* eigenvector_csf = eigenvector.data();
    for(auto& rconfig: *target.configs)
    {
        double contribution = 0.0;
        for(unsigned int j = 0; j < rconfig.NumCSFs(); j++)
        {   contribution += (*eigenvector_csf) * (*eigenvector_csf);
            eigenvector_csf++;
        }

        target_occ += rconfig * contribution;
    }

    AutoionizationConfigurationAveraged(target_occ);
}

void Atom::AutoionizationConfigurationAveraged(const OccupationMap& target)
{
    double ionization_energy = user_input("DR/IonizationEnergy", 0.0);
    double energy_limit = user_input("DR/EnergyLimit", -1.0);   // Default no limit
    int min_continuum_l = user_input("DR/ContinuumLMin", 0);
    int max_continuum_l = user_input("DR/ContinuumLMax", 6);    // Default maximum L = 6
    double grid_min = user_input("DR/EnergyGrid/Min", 0.003675);  // Default 0.1eV
    bool use_single_particle_energy = user_input.search("DR/--single-particle-energy");

    auto math = MathConstant::Instance();

    ConfigGenerator gen(orbitals, user_input);
    pOrbitalMap valence(new OrbitalMap(*orbitals->valence));

    *outstream << "\nTarget occupancy:\n  " << target.Name() << std::endl;

    // Create relativistic target
    RelativisticConfiguration rtarget;
    for(const auto& pair: target)
    {
        rtarget.insert(std::make_pair(pair.first, int(std::floor(target.GetOccupancy(pair.first) + 0.1))));
    }

    // Select orbitals in target wavefunction
    std::set<OrbitalInfo> target_info_set;
    pOrbitalMap target_map(new OrbitalMap(lattice));
    for(const auto& pair: target)
    {
        target_info_set.insert(pair.first);
    }
    for(const auto& info: target_info_set)
    {
        target_map->AddState(valence->GetState(info));
    }
    // Add holes from compound states
    target_map->AddStates(*orbitals->hole);

    // Create operator for construction of continuum field
    pHFOperator hf_continuum;
    std::string config = user_input("DR/ContinuumResidue", "");
    if(config == "")
        hf_continuum = hf_open;
    else
    {   hf_continuum = hf_open->Clone();
        pCore continuum_core(new Core(lattice, config));
        for(auto& orbital: *continuum_core)
        {
            pOrbital basis_state = orbitals->all->GetState(orbital.first);
            continuum_core->AddState(basis_state);
        }

        hf_continuum->SetCore(continuum_core);
    }

    // HF solver for making continuum
    pODESolver ode_solver(new AdamsSolver(hf->GetIntegrator()));
    HartreeFocker HF_Solver(ode_solver);

    // Get target including core for occupancy
    OccupationMap target_with_core;
    for(const auto& orb: *orbitals->core)
    {
        target_with_core.insert(std::make_pair(orb.first, orb.first.MaxNumElectrons()));
    }
    target_with_core += target;

    *outstream << "\nAutoionization strength:"
               << "\n E(eV)   S(ns)   Configuration" << std::endl;

    double energy_unit_conversion = math->HartreeEnergyIneV();
    double rate_unit_conversion = math->AtomicFrequencySI() * 1.e-9;    // (ns-1)

    // Loop over all configurations
    for(auto& rconfig: *allconfigs)
    {
        // Subtract target from configurations to get differences
        auto rdiff = rconfig - rtarget;

        // Check that this nrconfig is a valid double-excitation of target
        if(rdiff.ParticleNumber() != 3)
            continue;

        Parity eps_parity = rconfig.GetParity() * rtarget.GetParity();

        // Build continuum
        double eps_energy = 0;
        if(use_single_particle_energy)
        {
            for(auto& pair: rdiff)
                eps_energy += valence->GetState(pair.first)->Energy() * pair.second;
        }
        else
        {   OccupationMap rcompound = target + rdiff;
            eps_energy = CalculateConfigurationAverageEnergy(rcompound, orbitals->valence, hf_electron, twobody_electron->GetIntegrals()) - ionization_energy;
        }

        if(energy_limit > 0.0 && (eps_energy > energy_limit))
            continue;
        double eps_energy_calculated = eps_energy;
        if(eps_energy_calculated <= grid_min)
            eps_energy_calculated = grid_min;

        // Get participating electrons:
        //     a -> h and b -> eps
        const OrbitalInfo* electrons[2] = {nullptr, nullptr};
        const OrbitalInfo* h = nullptr;

        bool double_occupancy_excitation = (rdiff.size() == 2);

        auto it = rdiff.begin();
        while(it != rdiff.end())
        {
            if(it->second < 0)
            {   h = &(it->first);
                it++;
            }
            else if(electrons[0] == nullptr)
            {   electrons[0] = &(it->first);
                if(!double_occupancy_excitation) // Do not increase iterator if it->second == 2.
                    it++;
            }
            else
            {   electrons[1] = &(it->first);
                it++;
            }
        }

        double rate = 0.0;

        for(int choose_electron_a = 0; choose_electron_a < (double_occupancy_excitation? 1: 2); choose_electron_a++)
        {
            const OrbitalInfo* a = electrons[choose_electron_a];
            const OrbitalInfo* b = electrons[1 - choose_electron_a];

            // Get largest k1 which inform eps_twoJ
            int max_k1_overall = a->L() + h->L();
            if(2 * max_k1_overall > a->TwoJ() + h->TwoJ())
                max_k1_overall -= 2;

            // Get epsilon
            // Minimum L and Parity -> Minimum J
            int min_eps_twoJ;
            if(eps_parity == (min_continuum_l%2? Parity::odd: Parity::even))
                min_eps_twoJ = mmax(2 * min_continuum_l - 1, 1);
            else
                min_eps_twoJ = 2 * min_continuum_l + 1;

            min_eps_twoJ = mmax(min_eps_twoJ, b->TwoJ() - 2 * max_k1_overall);

            // Maximum L and Parity -> Maximum J
            int max_eps_twoJ;
            if(eps_parity == (max_continuum_l%2? Parity::odd: Parity::even))
                max_eps_twoJ = 2 * max_continuum_l + 1;
            else
                max_eps_twoJ = 2 * max_continuum_l - 1;

            max_eps_twoJ = mmin(max_eps_twoJ, b->TwoJ() + 2 * max_k1_overall);

            for(int eps_twoJ = min_eps_twoJ; eps_twoJ <= max_eps_twoJ; eps_twoJ+=2)
            {
                // Kappa = (-1)^(j+1/2 + l) (j + 1/2) = (-1)^(j+1/2).P.(j+1/2)
                int eps_kappa = (eps_twoJ + 1)/2;
                eps_kappa = math->minus_one_to_the_power(eps_kappa) * Sign(eps_parity) * eps_kappa;

                OrbitalInfo eps_info(100, eps_kappa);
                pContinuumWave eps(new ContinuumWave(eps_kappa, 100, eps_energy_calculated));
                HF_Solver.CalculateContinuumWave(eps, hf_continuum);

                // Create two-body integrals
                pOrbitalMap continuum_map(new OrbitalMap(lattice));
                continuum_map->AddState(eps);

                pOrbitalManager all_orbitals(new OrbitalManager(*orbitals));
                all_orbitals->all->AddState(eps);
                all_orbitals->MakeStateIndexes();
                pSlaterIntegrals two_body_integrals(new SlaterIntegralsDenseHash(all_orbitals, hartreeY));
                two_body_integrals->clear();
                two_body_integrals->CalculateTwoElectronIntegrals(continuum_map, target_map, valence, valence);

                // Create operator for autoionization
                pTwoElectronCoulombOperator two_body_operator = std::make_shared<TwoElectronCoulombOperator>(two_body_integrals);

                // Sum_k1 < a, b || V^k1 || h, eps > * multiplier
                // multiplier = < a, b || V^k1 || h, eps > / (2k1 + 1)
                //              - Sum_k2 (-1)^(k1 + k2 + 1) * 6jSym * < a, b || V^k2 || eps, h >
                // 6jSym = { k1 jb jeps }
                //         { k2 ja  jh  }

                int min_k1 = mmax(abs(b->L() - eps_info.L()), abs(a->L() - h->L()));
                int max_k1 = mmin(b->L() + eps_info.L(), abs(a->L() + h->L()));

                for(int k1 = min_k1; k1 <= max_k1; k1 += 2)
                {
                    double Vah = two_body_operator->GetReducedMatrixElement(k1, *a, *b, *h, eps_info);

                    if(fabs(Vah) > 1.e-15)
                    {
                        double multiplier = Vah / (2*k1 + 1);

                        // Exchange part
                        int min_k2 = mmax(abs(a->L() - eps_info.L()), abs(b->L() - h->L()));
                        int max_k2 = mmin(a->L() + eps_info.L(), abs(b->L() + h->L()));

                        for(int k2 = min_k2; k2 <= max_k2; k2 += 2)
                        {
                            double partial_k2 = math->Wigner6j(k1, b->J(), eps->J(), k2, a->J(), h->J());
                            if(partial_k2)
                            {
                                double Vbh = two_body_operator->GetReducedMatrixElement(k2, *a, *b, eps_info, *h);
                                partial_k2 *= math->minus_one_to_the_power(k1 + k2) * Vbh;
                            }

                            multiplier += partial_k2;
                        }

                        rate += Vah * multiplier;
                    }
                }
            }
        }

        rate *= 2. * math->Pi() * target_with_core.GetOccupancy(*h)/h->MaxNumElectrons()
                * (1. - target_with_core.GetOccupancy(*electrons[0])/electrons[0]->MaxNumElectrons())
                * (1. - target_with_core.GetOccupancy(*electrons[1])/electrons[1]->MaxNumElectrons());

        if(ProcessorRank == 0)
        {
            char sep = ' ';
            *outstream << std::setprecision(8) << eps_energy * energy_unit_conversion << sep
                       << rate * rate_unit_conversion << sep;

            *outstream << rconfig.Name('_') << std::endl;
        }
    }
}
}

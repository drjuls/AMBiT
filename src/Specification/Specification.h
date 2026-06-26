#ifndef SPECIFICATION_CLASS_H
#define SPECIFICATION_CLASS_H

#include <string>
#include <vector>
#include <optional>
#include <memory>
#include "SpecificationMap.h"
#include "TransitionSpecifications.h"
#include "parapara/parapara.h"

namespace P = parapara;

namespace Ambit
{

/************************** CI ***********************/
class CIBaseSpecification
{
public:
    std::vector<std::string> leading_configurations;
    std::optional<std::vector<std::string>> leading_rel_configurations;
    std::optional<std::vector<std::string>> extra_configurations;
    std::optional<std::vector<std::string>> extra_rel_configurations;

    unsigned electron_excitations = 2;
    std::optional<std::vector<std::string>> electron_bounds;
    unsigned hole_excitations = 0;
    std::optional<std::vector<std::string>> hole_bounds;
    bool print_configurations = false;
    bool print_rel_configurations = false;
    std::optional<std::pair<double, double> > configuration_average_energy_range;
};

// This is used for the "full" CI options, as it includes extra stuff not in SmallSide (mainly
// around how many and which solutions to generate
class CISpecification: public CIBaseSpecification
{
public:
    std::vector<unsigned> even_parity_twoj;
    std::vector<unsigned> odd_parity_twoj;
    unsigned num_solutions = 6;
    bool all_symmetries = false;
    bool gfactors = false;
    bool no_gfactors = false;
    bool memory_saver = false;
    bool single_configuration_ci = false;
    bool scalapack = false;
    std::optional<double> max_energy;
    unsigned chunksize = 4;
    bool sort_matrix_by_configuration = false;
    bool output_print_hamiltonian = false;
    bool output_write_hamiltonian = false;
    std::optional<double> output_max_displayed_energy;
    std::optional<double> output_min_displayed_percentage;
    bool output_print_inline = false;
    std::optional<std::string> output_separator;
    bool output_print_relativistic_configurations = false;
    bool output_no_configurations = false;
};

/************************** Deprecated ***********************/
// Holds keywords that older input files may still set but that AMBiT no longer supports.
// Listing them here lets the parser accept the keys without raising "unknown key" errors;
// validateAndNormaliseSpecification() then rejects any that were actually used with a message
// pointing the user at the replacement.
//
// To deprecate a new keyword: add a field below (use std::optional<T> for value-bearing keys,
// plain bool for "--"/"-" flags), register it in MakeGlobalSpecificationsMappings(), and add
// the rejection branch in validateAndNormaliseSpecification().
class DeprecatedSpecification
{
public:
    // Old internal multi-run driver, now replaced by shell scripts (commit e283584).
    std::optional<std::string> multirun;
    bool r_flag = false;

    // HF mass-shift keys relocated under HF/MassShift/ (commit 54ab3c2).
    std::optional<double> nuclear_inverse_mass;
    bool hf_sms = false;
    bool hf_nms = false;
    bool hf_only_rel_nms = false;
    bool hf_nonrel_mass_shift = false;
    bool hf_include_lower_sms = false;

    // HF / Basis structural removals.
    std::optional<unsigned> hf_charge;
    std::optional<std::string> basis_size;
    bool basis_custom_basis = false;
    bool basis_r_basis = false;
    std::optional<std::vector<std::string>> basis_custom_orbitals;

    // MBPT renames.
    bool mbpt_use_subtraction = false;
    bool mbpt_brueckner_excited = false;

    // CI/Output renames.
    bool ci_output_no_configs = false;
    std::optional<double> ci_output_min_displayed_percent;
};

/************************** GlobalSpecification ***********************/

class GlobalSpecification {

public:
    static GlobalSpecification* Instance();

    // Perform global validation of specifications. Return non-empty error message on failure.
    std::string validateAndNormaliseSpecification();
    P::hopefully<void> PrintSpecification();

    // Reset all fields to their default values (useful for test isolation).
    void Reset();

    // RAII guard that redirects Instance() to a caller-supplied object for the duration of a scope.
    // Nesting is supported: the previous override (or nullptr) is restored on destruction.
    class ScopedOverride {
    public:
        explicit ScopedOverride(GlobalSpecification* gs) : previous_(s_override) { s_override = gs; }
        ~ScopedOverride() { s_override = previous_; }
        ScopedOverride(const ScopedOverride&) = delete;
        ScopedOverride& operator=(const ScopedOverride&) = delete;
    private:
        GlobalSpecification* previous_;
    };

    GlobalSpecification();

    // Helper functions for parsing
    //template <typename R, typename X>
    //void add_delegates_for(std::vector<P::specification<R>>& specs, X R::* field_ptr, std::string key_pfx, const std::vector<P::specification<X>>& x_specs, std::string desc_pfx = "");

    // CI
    template <typename T>
    std::vector<P::specification<T>> make_ci_base_specification();

    template <typename T>
    std::vector<P::specification<T>> make_ci_specification();

    // Transitions
    template <typename T>
    std::vector<P::specification<T>> make_transitions_base_specification();

    template <typename T>
    std::vector<P::specification<T>> make_transitions_ej_specification();
    template <typename T>
    std::vector<P::specification<T>> make_transitions_hfs1_specification();
    template <typename T>
    std::vector<P::specification<T>> make_transitions_hfs2_specification();
    template <typename T>
    std::vector<P::specification<T>> make_transitions_hfi_specification();
    template <typename T>
    std::vector<P::specification<T>> make_transitions_fs_specification();
    template <typename T>
    std::vector<P::specification<T>> make_transitions_qed_specification();
    template <typename T>
    std::vector<P::specification<T>> make_transitions_yukawa_specification();
    template <typename T>
    std::vector<P::specification<T>> make_transitions_ke_specification();
    template <typename T>
    std::vector<P::specification<T>> make_transitions_llit2_specification();
    template <typename T>
    std::vector<P::specification<T>> make_transitions_nms_specification();
    template <typename T>
    std::vector<P::specification<T>> make_ic_specification();
    template <typename T>
    std::vector<P::specification<T>> make_residual_coulomb_specification();
    template <typename T>
    std::vector<P::specification<T>> make_deprecated_specification();

public:
    // Make the vector containing pairs mapping from the CLI option to their corresponding
    // fields in a GlobalSpecification struct.
    void MakeGlobalSpecificationsMappings();

public:
    // Ungrouped options
    std::string ID;
    unsigned Z = 0;
    std::optional<std::string> level_directory;
    bool only_verify_input = false;
    bool check_sizes = false;
    bool clean_run = false;
    bool print_basis = false;
    bool no_new_mbpt = false;
    bool ci_complete = false;
    bool no_ci = false;
    bool configuration_average = false;
    bool no_diagonalisation = false;
    double alpha_squared_variation = 0;
    double nuclear_radius = 4.0;
    double nuclear_thickness = 2.3;
    // Lattice
    unsigned lattice_num_points = 1000;
    double lattice_start_point = 1.e-6;
    double lattice_end_point = 50.0;
    double lattice_beta = 4.0;
    bool lattice_exponential = false;
    double lattice_h = 0.05;
    // HF
    std::optional<unsigned> hf_N;
    std::optional<std::string> hf_configuration;
    bool hf_breit = false;
    bool hf_qed = false; // Alias for hf_qed_uehling and hf_qed_self_energy
    bool hf_local_exchange = false;
    double hf_xalpha = 1.0;
    // HF/QED
    bool hf_do_qed = false; // Checks if the [HF/QED] section is defined
    bool hf_qed_uehling = false;
    bool hf_qed_self_energy = false;
    bool hf_qed_use_nuclear_density = false;
    std::optional<double> hf_qed_nuclear_rms_radius;
    bool hf_qed_no_magnetic = false;
    bool hf_qed_no_electric = false;
    bool hf_qed_skip_offmass = false;
    bool hf_qed_use_electron_screening = false;
    // HF/MassShift
    double hf_ms_nuclear_inverse_mass = 0;
    bool hf_ms_sms = false;
    bool hf_ms_nms = false;
    bool hf_ms_only_rel_nms = false;
    bool hf_ms_nonrel_mass_shift = false;
    bool hf_ms_include_lower_sms = false;
    // HF/NuclearPolarisability
    bool hf_do_nuclear_polarisability = false;
    double hf_nuclear_polarisability_alpha_e = 0.0;
    double hf_nuclear_polarisability_ebar_mev = 8.0;
    // HF/Yukawa
    bool hf_do_yukawa = false;
    std::optional<double> hf_yukawa_mass;
    std::optional<double> hf_yukawa_massev;
    std::optional<double> hf_yukawa_rc;
    double hf_yukawa_scale = 1.0;
    // HF/AddLocalPotential
    bool hf_do_local_potential = false;
    std::string hf_addlocal_filename;
    double hf_addlocal_scale = 1.0;
    // Basis
    std::string basis_valence;
    std::optional<std::string> basis_frozen_core;
    std::vector<std::string> basis_include_valence;
    std::vector<std::string> basis_exclude_valence;
    std::optional<std::string> basis_residue;
    std::optional<std::vector<std::string>> basis_inject_orbitals;
    bool basis_reorthogonalise = false;
    bool basis_hf = false;
    bool basis_bspline = false;
    bool basis_xr = false;
    std::optional<std::string> basis_hf_orbitals;
    std::optional<std::vector<std::string> > basis_xr_custom_orbitals;
    // Basis/BSpline
    std::optional<double> basis_bspline_rmax;
    double basis_bspline_r0 = 0;
    unsigned basis_bspline_k = 7;
    unsigned basis_bspline_N = 40;
    std::string basis_bspline_splinetype = "Reno";
    // CI
    bool do_ci = false;
    CISpecification ci;
    // CI/Output is also in CISpecification
    // CI/SmallSide
    bool ci_do_smallside = false;
    CIBaseSpecification ci_smallside;
    // MBPT
    std::optional<std::string> mbpt_basis;
    std::optional<std::string> mbpt_energy_denom_orbitals;
    bool mbpt_use_valence = false;
    bool mbpt_no_core = false;
    bool mbpt_force_subtraction = false;
    bool mbpt_no_subtraction = false;
    bool mbpt_no_extra_box = false;
    double mbpt_energy_denom_floor = 0.01;
    double mbpt_delta = 0.0;
    std::optional<std::vector<unsigned> > mbpt_twobody_storage_limits;
    std::optional<std::vector<double> > mbpt_onebody_scaling;
    bool mbpt_brueckner = false;
    bool mbpt_use_brueckner = false;
    // NOTE: Always use these when checking which MBPT diagrams to include
    bool mbpt_one_body = false;
    bool mbpt_two_body = false;
    bool mbpt_three_body = false;
    // MBPT/Brueckner
    double mbpt_brueckner_startpoint = 4.35e-5;
    double mbpt_brueckner_endpoint = 8.0;
    unsigned mbpt_brueckner_stride = 4;
    std::vector<double> mbpt_brueckner_scaling;
    bool mbpt_brueckner_use_lower = false;
    bool mbpt_brueckner_use_lower_lower = false;

    // Bools for whether to do a given transition type
    bool transitions_do_e1 = false;
    bool transitions_do_e2 = false;
    bool transitions_do_e3 = false;
    bool transitions_do_e4 = false;
    bool transitions_do_m1 = false;
    bool transitions_do_m2 = false;
    bool transitions_do_m3 = false;
    bool transitions_do_m4 = false;
    bool transitions_do_hfs1 = false;
    bool transitions_do_hfs2 = false;
    bool transitions_do_hfi  = false;
    bool transitions_do_fs = false;
    bool transitions_do_qed = false;
    bool transitions_do_yukawa = false;
    bool transitions_do_ke = false;
    bool transitions_do_llit2 = false;
    bool transitions_do_nms = false;
    bool transitions_do_residual_coulomb = false;

    // Transitions specification objects
    // EJ transitions
    TransitionsEJSpecification transitions_e1;
    TransitionsEJSpecification transitions_e2;
    TransitionsEJSpecification transitions_e3;
    TransitionsEJSpecification transitions_e4;
    TransitionsEJSpecification transitions_m1;
    TransitionsEJSpecification transitions_m2;
    TransitionsEJSpecification transitions_m3;
    TransitionsEJSpecification transitions_m4;

    // Everything else
    TransitionsHFS1Specification transitions_hfs1;
    TransitionsHFS2Specification transitions_hfs2;
    TransitionsHFISpecification transitions_hfi;
    TransitionsFSSpecification transitions_fs;
    TransitionsQEDSpecification transitions_qed;
    TransitionsYukawaSpecification transitions_yukawa;
    TransitionsKESpecification transitions_ke;
    TransitionsLLIT2Specification transitions_llit2;
    TransitionsNMSSpecification transitions_nms;

    // Autoionization
    ICSpecification ic;
    ResidualCoulombSpecification residual_coulomb;

    // Recombination
    std::string dr_target_file;
    std::optional<int> dr_target_twoj;
    std::string dr_target_parity;
    int dr_target_index = 0;

    std::optional<double> dr_ionization_energy;
    std::optional<double> dr_energy_limit;
    int dr_continuum_l_min = 0;
    int dr_continuum_l_max = 6;
    std::optional<std::string> dr_continuum_residue;

    bool dr_use_energygrid = false;
    double dr_energygrid_min = 0.003675;
    double dr_energygrid_max = 3.675;
    double dr_energygrid_step = 0.03675;
    bool dr_use_single_particle_energy = false;

    // Deprecated keys (accepted by the parser, rejected during validation)
    DeprecatedSpecification deprecated;
private:
    static GlobalSpecification* s_override;

protected:
    bool _s1 = false;
    bool _s2 = false;
    bool _s3 = false;
    bool _s123 = false;
    bool _s12 = false;
    bool _s13 = false;
    bool _s23 = false;
public:
    std::vector<P::specification<GlobalSpecification>> global_specifications_mappings;
};


// Miscellaneous utility functions

// On success, return empty string.
// On failure, return (long) error message.
std::string importSpecificationFile(GlobalSpecification* gs, const std::string& fileName);
std::string importSpecificationINIString(GlobalSpecification* gs, const std::string& ini_string);
std::string importSpecificationKV(GlobalSpecification* gs, const std::string& assignment);

// Return a map-like interface to the configuration, which lets us do nice things like
// config_map["key"] = "value"
typedef keyed_record_view<GlobalSpecification> SpecificationMap;
typedef std::shared_ptr<SpecificationMap> pSpecificationMap;
SpecificationMap get_config_map_view(GlobalSpecification& gs);

} // namespace Ambit

#endif //ndef SPECIFICATION_CLASS_H

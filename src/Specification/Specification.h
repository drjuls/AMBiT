#ifndef SPECIFICATION_CLASS_H
#define SPECIFICATION_CLASS_H

#include <string>
#include <vector>
#include <optional>
#include <memory>
#include "SpecificationMap.h"
#include "parapara/parapara.h"

namespace P = parapara;

namespace Ambit
{

// Generic options for CI and CI/SmallSide
class CIBaseSpecification
{
public:
    std::string leading_configurations;
    std::string leading_rel_configurations;
    std::string extra_configurations;
    std::string extra_rel_configurations;

    /* TODO EVK: This is really annoying, because the spec says that it can be *either* an integer
     * (e.g. CI/ElectronExcitations=2) or a string (e.g. ElectronExcitations = '1,5spdf,2,5spd')
     * and this is really annoying to deal with. My current approach of splitting it into two fields
     * (one for the number of excitations, and one for the bounds) completely breaks backwards
     * compatibility with the old style of user input parsing.
     *
     * The alternative would be to use a std::variant, which would rule because
     * algebraic types rule, but std::variant's actual implementation in C++ is kind of yucky.
     * Also, the existing "list" form of arguments ('1, 5spdf, 2, 5spd') also has multiple types
     * within a single list (an integer followed by a string), which kind of inherently requires
     * weird type wrangling to get it to work.
    */
    unsigned electron_excitations;
    std::optional<std::vector<std::string>> electron_bounds;
    std::optional<std::vector<std::string>> hole_bounds;
    unsigned hole_excitations = 0;
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
    std::optional<double> output_min_displayed_percent;
    bool output_print_inline = false;
    std::optional<std::string> output_separator;
    bool output_print_relativistic_configurations = false;
    bool output_no_configs = false;
};

class GlobalSpecification {

public:
    static GlobalSpecification* Instance();

    /*
    // On success, return empty string.
    // On failure, return (long) error message.

    std::string importSpecificationFile(const std::string& fileName);
    std::string importSpecificationKV(const std::string& assignment);
    */

    // Perform global validation of specifications. Return non-empty error message on failure.
    std::string validateAndNormaliseSpecification();

protected:
    GlobalSpecification();

    // Helper functions for parsing
    //template <typename R, typename X>
    //void add_delegates_for(std::vector<P::specification<R>>& specs, X R::* field_ptr, std::string key_pfx, const std::vector<P::specification<X>>& x_specs, std::string desc_pfx = "");
    
    template <typename T>
    std::vector<P::specification<T>> make_ci_base_specification();

    template <typename T>
    std::vector<P::specification<T>> make_ci_specification();


public:
    // Make the vector containing pairs mapping from the CLI option to their corresponding
    // fields in a GlobalSpecification struct.
    void MakeGlobalSpecificationsMappings();

public: 
    // Ungrouped options
    std::string ID;
    unsigned Z = 0;
    std::optional<std::string> level_directory;
    bool no_new_mbpt = false;
    bool check_sizes = false;
    bool clean_run = false;
    bool print_basis = false;
    bool ci_complete = false;
    bool no_ci = false;
    bool configuration_average = false;
    double alpha_squared_variation = 0;
    double nuclear_inverse_mass = 0;
    double nuclear_radius = 0;
    double nuclear_thickness = 0;
    bool m = false;
    // Lattice
    unsigned lattice_num_points = 0;
    double lattice_start_point = 0;
    double lattice_end_point = 0;
    bool lattice_exponential = false;
    double lattice_H = 0.05;
    // HF
    std::optional<unsigned> hf_N;
    std::optional<int> hf_charge;
    std::string hf_configuration;
    bool hf_breit = false;
    bool hf_sms = false;
    bool hf_nms = false;
    bool hf_only_rel_nms = false;
    bool hf_nonrel_mass_shift = false;
    bool hf_include_lower_mass = false;
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
    // HF/NuclearPolarisability
    bool hf_do_nuclear_polarisability = false;
    double hf_nuclear_polarisability_alpha_e;
    double hf_nuclear_polarisability_ebar_mev;
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
    std::optional<std::string> basis_size;
    // TODO EVK: Should this be a std::optional?
    std::vector<std::string> basis_include_valence;
    std::vector<std::string> basis_exclude_valence;
    std::optional<std::string> basis_residue;
    std::optional<std::vector<std::string>> basis_inject_orbitals;
    bool basis_reorthogonalise = false;
    bool basis_hf = false;
    bool basis_bspline = true;
    bool basis_xr = false;
    std::optional<std::string> basis_hf_orbitals;
    std::optional<std::vector<std::string> > basis_xr_custom_orbitals;
    // Basis/BSpline
    std::optional<double> basis_bspline_rmax; 
    double basis_bspline_r0 = 0;
    double basis_bspline_k = 7;
    double basis_bspline_N = 40;
    std::string basis_bspline_splinetype;
    // CI
    CISpecification ci;
    // CI/Output is also in CISpecification
    // CI/SmallSide
    bool ci_do_smallside;
    CIBaseSpecification ci_smallside;
    // MBPT
    std::optional<std::string> mbpt_basis;
    std::optional<std::string> mbpt_energy_denom_orbitals;
    bool mbpt_use_valence;
    bool mbpt_no_core;
    bool mbpt_use_subtraction;
    bool mbpt_no_subtraction;
    bool mbpt_no_extra_box;
    double mbpt_energy_denom_floor = 0.01;
    double mbpt_delta = 0.0;
    std::optional<std::vector<unsigned> > mbpt_twobody_storage_limits;
    std::optional<std::vector<double> > mbpt_onebody_scaling;
    bool mbpt_brueckner = false;
    // NOTE: Always use these when checking which MBPT diagrams to include
    bool mbpt_one_body = false;
    bool mbpt_two_body = false;
    bool mbpt_three_body = false;
    // MBPT/Brueckner
    double mbpt_brueckner_startpoint = 4.35e-5;
    double mbpt_brueckner_endpoint = 8.0;
    unsigned mbpt_brueckner_stride = 4;
    std::vector<double> mbpt_brueckner_scaling;
    std::vector<double> mbpt_brueckner_energy_scaling;
    bool mbpt_brueckner_use_lower = false;
    bool mbpt_brueckner_use_lower_lower = false;
    bool mbpt_brueckner_excited = false;
    // We need to explicitly enumerate the options for MBPT diagrams: s1, s2, s3, s12, s123, s13
    // and s23. It's kind of clunky to use these directly in the main body of the code, so we
    // obfuscate them here and then, during the spec validation, set the more general attributes:
    // mbpt_one_body, mbpt_two_body, mbpt_three_body as appropriate
    // NOTE: Don't use these directly in AMBiT, use GlobalSpecification::mbpt_one_body and friends
    // TODO EVK: Really want to be able to make these private, but GlobalSpecification is declared
    // globally, so it doesn't work. Probably need to move the parsing stuff into a separate class
    // so I can make it a friend
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
    // Nice std::map-like interface for retrieving global_specifications fields based on
    // their CLI arguments
    //P::specification_map<GlobalSpecification> global_specifications_dict;
};

// On success, return empty string.
// On failure, return (long) error message.
std::string importSpecificationFile(GlobalSpecification* gs, const std::string& fileName);
std::string importSpecificationKV(GlobalSpecification* gs, const std::string& assignment);

// Perform global validation of specifications. Return non-empty error message on failure.
//std::string validateAndNormaliseSpecification(GlobalSpecification* gs);


// Return a map-like interface to the configuration, which lets us do nice things like 
// config_map["key"] = "value"
typedef keyed_record_view<GlobalSpecification> SpecificationMap;
typedef std::shared_ptr<SpecificationMap> pSpecificationMap;
SpecificationMap get_config_map_view(GlobalSpecification& gs);

} // namespace Ambit

#endif //ndef SPECIFICATION_CLASS_H

#include <algorithm>
#include <fstream>
#include <istream>
#include <sstream>
#include <string>
#include <string_view>
#include <vector>

#include "Specification.h"
#include "Specification/TransitionSpecifications.h"
#include "parapara/parapara.h"
#include "Include.h"

namespace P = parapara;

namespace Ambit
{

// Use a custom line-by-line ini importer to handle relative section headings
template <typename Record>
P::hopefully<void> custom_import_ini(Record& rec, const P::specification_map<Record>& specs, std::istream& in);

template <typename Value>
auto between(Value min, Value max, std::string constraint = "value between upper and lower bounds (x in (min, max] )")
{   return P::validator([min = std::move(min), max = std::move(max)](auto x) { return (x>min && x <= max); }, std::move(constraint));
}

inline auto positive = P::greater_than(0, "must be positive");
inline auto nonnegative = P::at_least(0, "must be positive");
inline auto zero_to_1 = between(0, 1, "must be positive");

/********************* GlobalSpecification methods *********************/
GlobalSpecification* GlobalSpecification::s_override = nullptr;

GlobalSpecification* GlobalSpecification::Instance()
{
    if (s_override) return s_override;
    static GlobalSpecification default_instance;
    return &default_instance;
}

GlobalSpecification::GlobalSpecification()
{
    // Create the mappings from CLI arguments to GlobalSpecification fields
    MakeGlobalSpecificationsMappings();
}

void GlobalSpecification::Reset()
{
    *this = GlobalSpecification();
}

/********************* Parsing methods *********************/
// Helper methods for the polymorphic "delegated" parsing (so we can e.g. re-use common options
// between [Transition] operators. This is lifted from "ex-delegate.cc" in the parapara
// examples: https://github.com/halfflat/parapara/blob/main/ex/ex-delegate.cc

// A bit of syntactic sugar for inheritance (even if it does involve some copies):
template <typename C>
C extend(C&& container, const std::initializer_list<typename C::value_type>& il) {
    std::copy(std::begin(il), std::end(il), std::back_inserter(container));
    return container;
}

/****************************** CI ****************************/
template <typename X>
std::vector<P::specification<X>> GlobalSpecification::make_ci_base_specification() {
    return {
        {"LeadingConfigurations",                &X::leading_configurations},
        {"LeadingRelativisticConfigurations",    &X::leading_rel_configurations},
        {"ExtraConfigurations",                  &X::extra_configurations},
        {"ExtraRelativisticConfigurations",      &X::extra_rel_configurations},
        {"ElectronExcitations",                  &X::electron_excitations},
        {"ElectronBounds",                       &X::electron_bounds},
        {"HoleExcitations",                      &X::hole_excitations},
        {"HoleBounds",                           &X::hole_bounds},
        {"--print-configurations",               &X::print_configurations},
        {"--print-relativistic-configurations",  &X::print_rel_configurations},
        {"ConfigurationAverageEnergyRange",      &X::configuration_average_energy_range},
    };
}

template <typename X>
std::vector<P::specification<X>> GlobalSpecification::make_ci_specification() {
    // Extend from the base specification class
    return extend(make_ci_base_specification<X>(), {
        // Now add the extras on top
        {"EvenParityTwoJ",                       &X::even_parity_twoj},
        {"OddParityTwoJ",                        &X::odd_parity_twoj},
        {"NumSolutions",                         &X::num_solutions},
        {"--all-symmetries",                     &X::all_symmetries},
        {"--gfactors",                           &X::gfactors},
        {"--no-gfactors",                        &X::no_gfactors},
        {"--memory-saver",                       &X::memory_saver},
        {"--single-configuration-CI",            &X::single_configuration_ci},
        {"--scalapack",                                     &X::scalapack},
        {"MaxEnergy",                                       &X::max_energy},
        {"ChunkSize",                                       &X::chunksize},
        {"--sort-matrix-by-configuration",                  &X::sort_matrix_by_configuration},
        {"Output/--write-hamiltonian",                      &X::output_write_hamiltonian},
        {"Output/--print-hamiltonian",                      &X::output_print_hamiltonian},
        {"Output/--print-inline",                           &X::output_print_inline},
        {"Output/--print-relativistic-configurations",      &X::output_print_relativistic_configurations},
        {"Output/--no-configurations",                      &X::output_no_configurations},
        {"Output/MaxDisplayedEnergy",                       &X::output_max_displayed_energy},
        {"Output/MinDisplayedPercentage",                   &X::output_min_displayed_percentage},
        {"Output/Separator",                                &X::output_separator},
     });
}

/****************************** Transitions ****************************/

template <typename X>
std::vector<P::specification<X>> GlobalSpecification::make_transitions_base_specification() {
    return {
        {"MatrixElements",                  &X::matrix_elements},
        {"AllBelow",                        &X::all_below},
        {"--reduced-elements",              &X::reduced_elements},
        {"--print-integrals",               &X::print_integrals},
        {"--rpa",                           &X::do_rpa},
        {"Scale",                           &X::scale},
        {"RPA/Scale",                       &X::rpa_scale},
        {"RPA/Weighting",                   &X::rpa_weighting, zero_to_1},
        {"RPA/MaxDeltaK",                   &X::rpa_max_delta_k, nonnegative},
        {"RPA/MaxIterations",               &X::rpa_max_iterations, positive},
        {"RPA/Configuration",               &X::rpa_configuration},
        {"RPA/--no-negative-states",        &X::rpa_no_negative_states},
        {"RPA/--print-field",               &X::rpa_print_field},
        {"RPA/BSpline/N",                   &X::rpa_bspline_n},
        {"RPA/BSpline/K",                   &X::rpa_bspline_k},
        {"RPA/BSpline/Rmax",                &X::rpa_bspline_rmax},
        {"RPA/BSpline/R0",                  &X::rpa_bspline_r0},
    };
}

template <typename X>
std::vector<P::specification<X>> GlobalSpecification::make_transitions_ej_specification() {
    // Start with the base transitions specification and add extras on top
    return extend(make_transitions_base_specification<X>(), {
        {"Frequency",                       &X::frequency},
        {"--velocity-gauge",                &X::velocity_gauge},
    });
}

template <typename X>
std::vector<P::specification<X>> GlobalSpecification::make_transitions_hfs1_specification() {
    // Start with the base transitions specification and add extras on top
    return extend(make_transitions_base_specification<X>(), {
        {"NuclearMagneticRadius",           &X::nuclear_magnetic_radius},
        {"gOnI",                            &X::g_on_i},
    });
}

template <typename X>
std::vector<P::specification<X>> GlobalSpecification::make_transitions_hfs2_specification() {
    // Start with the base transitions specification and add extras on top
    return extend(make_transitions_base_specification<X>(), {
        {"NuclearQuadrupoleRadius",         &X::nuclear_quadrupole_radius},
        {"Q",                               &X::q},
    });
}

template <typename X>
std::vector<P::specification<X>> GlobalSpecification::make_transitions_hfi_specification() {
    // Start with the base transitions specification and add extras on top
    return extend(make_transitions_base_specification<X>(), {
        {"NuclearRadius",           &X::nuclear_radius},
        {"Operator",                &X::op},
    });
}

template <typename X>
std::vector<P::specification<X>> GlobalSpecification::make_transitions_fs_specification() {
    // Start with the base transitions specification and add extras on top
    return extend(make_transitions_base_specification<X>(), {
        {"DeltaNuclearRadius",              &X::delta_nuclear_radius},
        {"DeltaNuclearThickness",           &X::delta_nuclear_thickness},
    });
}

template <typename X>
std::vector<P::specification<X>> GlobalSpecification::make_transitions_qed_specification() {
    // Start with the base transitions specification and add extras on top
    return extend(make_transitions_base_specification<X>(), {
    {"QED/--uehling",                   &X::uehling},
    {"QED/--self-energy",               &X::self_energy},
    {"QED/--use-nuclear-density",       &X::use_nuclear_density},
    {"QED/NuclearRMSRadius",            &X::nuclear_rms_radius, P::nonzero()},
    {"QED/--no-magnetic",               &X::no_magnetic},
    {"QED/--no-electric",               &X::no_electric},
    {"QED/--skip-offmass",              &X::skip_offmass},
    {"QED/--use-electron-screening",    &X::use_electron_screening},
    });
}

template <typename X>
std::vector<P::specification<X>> GlobalSpecification::make_transitions_yukawa_specification() {
    // Start with the base transitions specification and add extras on top
    return extend(make_transitions_base_specification<X>(), {
        {"Mass",            &X::mass},
        {"MassEV",          &X::massev},
        {"Rc",              &X::rc},
    });
}

// Don't need anything for KE or LLIT2 because they have no options beyond the base spec
template <typename X>
std::vector<P::specification<X>> GlobalSpecification::make_transitions_nms_specification() {
    // Start with the base transitions specification and add extras on top
    return extend(make_transitions_base_specification<X>(), {
        {"--only-relativistic-nms",         &X::only_relativistic_nms},
        {"--nonrelativistic-mass-shift",    &X::nonrelativistic_mass_shift},
    });
}

template <typename X>
std::vector<P::specification<X>> GlobalSpecification::make_ic_specification() {
    // Start with the base transitions specification and add extras on top
    return extend(make_transitions_hfi_specification<X>(), {
        {"Source/Filename",         &X::source_filename},
        {"Source/TwoJ",             &X::source_twoj, nonnegative},
        {"Source/Parity",           &X::source_parity},
        {"Source/Index",            &X::source_index, nonnegative},
        {"NuclearEnergyEV",         &X::nuclear_energy_ev, positive},
        {"ContinuumLMin",           &X::continuum_l_min},
        {"ContinuumLMax",           &X::continuum_l_max},
    });
}

template <typename X>
std::vector<P::specification<X>> GlobalSpecification::make_residual_coulomb_specification() {
    // Start with the base transitions specification and add extras on top
    return make_transitions_base_specification<X>();
}

template <typename X>
std::vector<P::specification<X>> GlobalSpecification::make_deprecated_specification() {
    // Deprecated keywords: parsed at the top level (empty key prefix), rejected during validation.
    return {
        {"Multirun",                             &X::multirun},
        {"-r",                                   &X::r_flag},
        {"NuclearInverseMass",                   &X::nuclear_inverse_mass},
        {"HF/--sms",                             &X::hf_sms},
        {"HF/--nms",                             &X::hf_nms},
        {"HF/--only-relativistic-nms",           &X::hf_only_rel_nms},
        {"HF/--nonrelativistic-mass-shift",      &X::hf_nonrel_mass_shift},
        {"HF/--include-lower-sms",               &X::hf_include_lower_sms},
        {"HF/Charge",                            &X::hf_charge},
        {"Basis/BasisSize",                      &X::basis_size},
        {"Basis/--custom-basis",                 &X::basis_custom_basis},
        {"Basis/--r-basis",                      &X::basis_r_basis},
        {"Basis/CustomOrbitals",                 &X::basis_custom_orbitals},
        {"MBPT/--use-subtraction",               &X::mbpt_use_subtraction},
        {"MBPT/Brueckner/--excited",             &X::mbpt_brueckner_excited},
        {"CI/Output/--no-configs",               &X::ci_output_no_configs},
        {"CI/Output/MinimumDisplayedPercentage", &X::ci_output_min_displayed_percent},
    };
}

/***************************** Delegated parsing *************************/
// Add specifications for members of record subobjects that have their own sets of specifications.
// Keys and field descriptions are given a prefix as supplied.
template <typename R, typename X>
void add_delegates_for(std::vector<P::specification<R>>& specs, X R::* field_ptr, std::string key_pfx, const std::vector<P::specification<X>>& x_specs, std::string desc_pfx = "")
{
    for (const auto& x_spec: x_specs) {
        specs.emplace_back(key_pfx+x_spec.key, field_ptr, x_spec, desc_pfx+x_spec.description);
    }
}

// Now for the Global options
void GlobalSpecification::MakeGlobalSpecificationsMappings()
{
    global_specifications_mappings = {
    // Ungrouped options
    {"ID",                      &GlobalSpecification::ID, P::nonempty()},
    {"Z",                       &GlobalSpecification::Z, P::at_least(1)},
    {"LevelDirectory",          &GlobalSpecification::level_directory},
    // NOTE: Don't check "GlobalSpecification::s1" and friends in the main body of the code, use
    // "GlobalSpecification::mbpt_one_body" etc instead; it's much cleaner to figure out which
    // integrals to include once when normalising the spec
    {"-s1",                     &GlobalSpecification::_s1},
    {"-s2",                     &GlobalSpecification::_s2},
    {"-s3",                     &GlobalSpecification::_s3},
    {"-s12",                    &GlobalSpecification::_s12},
    {"-s13",                    &GlobalSpecification::_s13},
    {"-s23",                    &GlobalSpecification::_s23},
    {"-s123",                   &GlobalSpecification::_s123},
    {"--no-new-mbpt",           &GlobalSpecification::no_new_mbpt},
    {"--check-sizes",           &GlobalSpecification::check_sizes},
    {"-c",                      &GlobalSpecification::clean_run},
    {"--clean",                 &GlobalSpecification::clean_run},
    {"-p",                      &GlobalSpecification::print_basis},
    {"--print-basis",           &GlobalSpecification::print_basis},
    {"--ci-complete",           &GlobalSpecification::ci_complete},
    {"--no-ci",                 &GlobalSpecification::no_ci},
    {"--configuration-average", &GlobalSpecification::configuration_average},
    {"--no-diagonalisation",    &GlobalSpecification::no_diagonalisation},
    {"NuclearRadius",           &GlobalSpecification::nuclear_radius},
    {"NuclearThickness",        &GlobalSpecification::nuclear_thickness},
    {"AlphaSquaredVariation",   &GlobalSpecification::alpha_squared_variation},
    // Lattice
    {"Lattice/NumPoints",       &GlobalSpecification::lattice_num_points,  positive},
    {"Lattice/StartPoint",      &GlobalSpecification::lattice_start_point, positive},
    {"Lattice/EndPoint",        &GlobalSpecification::lattice_end_point,   positive},
    {"Lattice/Beta",            &GlobalSpecification::lattice_beta,        positive},
    {"Lattice/--exp-lattice",   &GlobalSpecification::lattice_exponential},
    {"Lattice/H",               &GlobalSpecification::lattice_h,           positive},
    // HF
    {"HF/N",                            &GlobalSpecification::hf_N,        nonnegative},
    {"HF/Configuration",                &GlobalSpecification::hf_configuration},
    {"HF/--breit",                      &GlobalSpecification::hf_breit},
    {"HF/--qed",                        &GlobalSpecification::hf_qed},
    {"HF/--local-exchange",             &GlobalSpecification::hf_local_exchange},
    {"HF/Xalpha",                       &GlobalSpecification::hf_xalpha},
    // HF/QED
    {"HF/QED",                          &GlobalSpecification::hf_do_qed},
    {"HF/QED/--uehling",                &GlobalSpecification::hf_qed_uehling},
    {"HF/QED/--self-energy",            &GlobalSpecification::hf_qed_self_energy},
    {"HF/QED/--use-nuclear-density",    &GlobalSpecification::hf_qed_use_nuclear_density},
    {"HF/QED/NuclearRMSRadius",         &GlobalSpecification::hf_qed_nuclear_rms_radius, P::nonzero()},
    {"HF/QED/--no-magnetic",            &GlobalSpecification::hf_qed_no_magnetic},
    {"HF/QED/--no-electric",            &GlobalSpecification::hf_qed_no_electric},
    {"HF/QED/--skip-offmass",           &GlobalSpecification::hf_qed_skip_offmass},
    {"HF/QED/--use-electron-screening", &GlobalSpecification::hf_qed_use_electron_screening},
    // HF/MassShift
    {"HF/MassShift/NuclearInverseMass",           &GlobalSpecification::hf_ms_nuclear_inverse_mass},
    {"HF/MassShift/--sms",                        &GlobalSpecification::hf_ms_sms},
    {"HF/MassShift/--nms",                        &GlobalSpecification::hf_ms_nms},
    {"HF/MassShift/--only-relativistic-nms",      &GlobalSpecification::hf_ms_only_rel_nms},
    {"HF/MassShift/--nonrelativistic-mass-shift", &GlobalSpecification::hf_ms_nonrel_mass_shift},
    {"HF/MassShift/--include-lower-sms",          &GlobalSpecification::hf_ms_include_lower_sms},
    // HF/NuclearPolarisability
    {"HF/NuclearPolarisability",        &GlobalSpecification::hf_do_nuclear_polarisability},
    {"HF/NuclearPolarisability/AlphaE", &GlobalSpecification::hf_nuclear_polarisability_alpha_e, P::nonzero()},
    {"HF/NuclearPolarisability/EbarMeV",&GlobalSpecification::hf_nuclear_polarisability_ebar_mev, P::nonzero()},
    // HF/Yukawa
    {"HF/Yukawa",           &GlobalSpecification::hf_do_yukawa},
    {"HF/Yukawa/Mass",      &GlobalSpecification::hf_yukawa_mass},
    {"HF/Yukawa/MassEV",    &GlobalSpecification::hf_yukawa_massev},
    {"HF/Yukawa/Rc",        &GlobalSpecification::hf_yukawa_rc},
    {"HF/Yukawa/Scale",     &GlobalSpecification::hf_yukawa_scale},
    // HF/AddLocalPotential
    {"HF/AddLocalPotential",            &GlobalSpecification::hf_do_local_potential},
    {"HF/AddLocalPotential/Filename",   &GlobalSpecification::hf_addlocal_filename},
    {"HF/AddLocalPotential/Scale",      &GlobalSpecification::hf_addlocal_scale},
    // Basis
    {"Basis/ValenceBasis",      &GlobalSpecification::basis_valence, P::nonempty()},
    {"Basis/FrozenCore",        &GlobalSpecification::basis_frozen_core},
    {"Basis/IncludeValence",    &GlobalSpecification::basis_include_valence},
    {"Basis/ExcludeValence",    &GlobalSpecification::basis_exclude_valence},
    {"Basis/Residue",           &GlobalSpecification::basis_residue},
    {"Basis/InjectOrbitals",    &GlobalSpecification::basis_inject_orbitals},
    {"Basis/--reorthogonalise", &GlobalSpecification::basis_reorthogonalise},
    {"Basis/--hf-basis",        &GlobalSpecification::basis_hf},
    {"Basis/--bspline-basis",   &GlobalSpecification::basis_bspline},
    {"Basis/--xr-basis",        &GlobalSpecification::basis_xr},
    {"Basis/HFOrbitals",        &GlobalSpecification::basis_hf_orbitals},
    // Basis/BSpline
    {"Basis/BSpline/Rmax",          &GlobalSpecification::basis_bspline_rmax},
    {"Basis/BSpline/R0",            &GlobalSpecification::basis_bspline_r0},
    {"Basis/BSpline/K",             &GlobalSpecification::basis_bspline_k},
    {"Basis/BSpline/N",             &GlobalSpecification::basis_bspline_N},
    {"Basis/BSpline/SplineType",    &GlobalSpecification::basis_bspline_splinetype},
    // Basis/xR
    {"Basis/xR/CustomOrbitals", &GlobalSpecification::basis_xr_custom_orbitals},
    // CI
    {"CI",              &GlobalSpecification::do_ci},
    {"CI/SmallSide",    &GlobalSpecification::ci_do_smallside},
    // MBPT
    {"MBPT/Basis",                  &GlobalSpecification::mbpt_basis},
    {"MBPT/EnergyDenomOrbitals",    &GlobalSpecification::mbpt_energy_denom_orbitals},
    {"MBPT/--use-valence",          &GlobalSpecification::mbpt_use_valence},
    {"MBPT/--no-core",              &GlobalSpecification::mbpt_no_core},
    {"MBPT/--force-subtraction",    &GlobalSpecification::mbpt_force_subtraction},
    {"MBPT/--no-subtraction",       &GlobalSpecification::mbpt_no_subtraction},
    {"MBPT/--no-extra-box",         &GlobalSpecification::mbpt_no_extra_box},
    {"MBPT/EnergyDenomFloor",       &GlobalSpecification::mbpt_energy_denom_floor},
    {"MBPT/Delta",                  &GlobalSpecification::mbpt_delta},
    {"MBPT/TwoBody/StorageLimits",  &GlobalSpecification::mbpt_twobody_storage_limits},
    {"MBPT/OneBody/Scaling",        &GlobalSpecification::mbpt_onebody_scaling},
    {"MBPT/--brueckner",            &GlobalSpecification::mbpt_brueckner},
    {"MBPT/--mbpt-use-brueckner",   &GlobalSpecification::mbpt_use_brueckner},
    {"MBPT/Brueckner/StartPoint",   &GlobalSpecification::mbpt_brueckner_startpoint},
    {"MBPT/Brueckner/EndPoint",     &GlobalSpecification::mbpt_brueckner_endpoint},
    {"MBPT/Brueckner/Stride",       &GlobalSpecification::mbpt_brueckner_stride},
    {"MBPT/Brueckner/Scaling",      &GlobalSpecification::mbpt_brueckner_scaling},
    {"MBPT/Brueckner/--use-lower",  &GlobalSpecification::mbpt_brueckner_use_lower},
    {"MBPT/Brueckner/--use-lower-lower",&GlobalSpecification::mbpt_brueckner_use_lower_lower},
    {"Transitions/E1",              &GlobalSpecification::transitions_do_e1},
    {"Transitions/E2",              &GlobalSpecification::transitions_do_e2},
    {"Transitions/E3",              &GlobalSpecification::transitions_do_e3},
    {"Transitions/E4",              &GlobalSpecification::transitions_do_e4},
    {"Transitions/M1",              &GlobalSpecification::transitions_do_m1},
    {"Transitions/M2",              &GlobalSpecification::transitions_do_m2},
    {"Transitions/M3",              &GlobalSpecification::transitions_do_m3},
    {"Transitions/M4",              &GlobalSpecification::transitions_do_m4},
    {"Transitions/HFS1",            &GlobalSpecification::transitions_do_hfs1},
    {"Transitions/HFS2",            &GlobalSpecification::transitions_do_hfs2},
    {"Transitions/HFI",             &GlobalSpecification::transitions_do_hfi},
    {"Transitions/FS",              &GlobalSpecification::transitions_do_fs},
    {"Transitions/QED",             &GlobalSpecification::transitions_do_qed},
    {"Transitions/Yukawa",          &GlobalSpecification::transitions_do_yukawa},
    {"Transitions/KE",              &GlobalSpecification::transitions_do_ke},
    {"Transitions/LLIT2",           &GlobalSpecification::transitions_do_llit2},
    {"Transitions/NMS",             &GlobalSpecification::transitions_do_nms},
    {"Transitions/ResidualCoulomb", &GlobalSpecification::transitions_do_residual_coulomb},
    {"DR/Target/Filename",          &GlobalSpecification::dr_target_file, P::nonempty()},
    {"DR/Target/TwoJ",              &GlobalSpecification::dr_target_twoj, nonnegative},
    {"DR/Target/Parity",            &GlobalSpecification::dr_target_parity},
    {"DR/Target/Index",             &GlobalSpecification::dr_target_index, nonnegative},
    {"DR/IonizationEnergy",         &GlobalSpecification::dr_ionization_energy, nonnegative},
    {"DR/EnergyLimit",              &GlobalSpecification::dr_energy_limit},
    {"DR/ContinuumLMin",            &GlobalSpecification::dr_continuum_l_min},
    {"DR/ContinuumLMax",            &GlobalSpecification::dr_continuum_l_max},
    {"DR/ContinuumResidue",         &GlobalSpecification::dr_continuum_residue},
    {"DR/--energy-grid",            &GlobalSpecification::dr_use_energygrid},
    {"DR/EnergyGrid/Min",           &GlobalSpecification::dr_energygrid_min},
    {"DR/EnergyGrid/Max",           &GlobalSpecification::dr_energygrid_max},
    {"DR/EnergyGrid/Step",          &GlobalSpecification::dr_energygrid_step},
    {"DR/--single-particle-energy", &GlobalSpecification::dr_use_single_particle_energy},
};

/**************************** Delegated parsing (CI and transitions) ***********/
    // CI and CI/SmallSide specs
    auto ci_specs = make_ci_specification<CISpecification>();
    auto ci_smallside_specs = make_ci_base_specification<CIBaseSpecification>();

    // Delegates
    add_delegates_for(global_specifications_mappings, &GlobalSpecification::ci, "CI/", ci_specs, "CI full: ");
    add_delegates_for(global_specifications_mappings, &GlobalSpecification::ci_smallside, "CI/SmallSide/", ci_smallside_specs, "CI SmallSide: ");

    // Transition matrix elements
    // Base specs common to EJ/MJ operators
    auto ej_specs = make_transitions_ej_specification<TransitionsEJSpecification>();
    auto hfs1_specs = make_transitions_hfs1_specification<TransitionsHFS1Specification>();
    auto hfs2_specs = make_transitions_hfs2_specification<TransitionsHFS2Specification>();
    auto hfi_specs = make_transitions_hfi_specification<TransitionsHFISpecification>();
    auto fs_specs = make_transitions_fs_specification<TransitionsFSSpecification>();
    auto qed_specs = make_transitions_qed_specification<TransitionsQEDSpecification>();
    auto yukawa_specs = make_transitions_yukawa_specification<TransitionsYukawaSpecification>();
    auto ke_specs = make_transitions_base_specification<TransitionsKESpecification>();
    auto llit2_specs = make_transitions_base_specification<TransitionsLLIT2Specification>();
    auto nms_specs = make_transitions_nms_specification<TransitionsNMSSpecification>();

    auto ic_specs = make_ic_specification<ICSpecification>();
    auto residual_coulomb_specs = make_residual_coulomb_specification<ResidualCoulombSpecification>();

    // Now delegate the parsing for each EJ multipole operator
    // NOTE: The trailing "/"s in the key prefixes are signification (e.g. it's
    // "Transitions/E1/" *not "Transitions/E1")
    add_delegates_for(global_specifications_mappings, &GlobalSpecification::transitions_e1, "Transitions/E1/", ej_specs, "E1 transition operator: ");
    add_delegates_for(global_specifications_mappings, &GlobalSpecification::transitions_e2, "Transitions/E2/", ej_specs, "E2 transition operator: ");
    add_delegates_for(global_specifications_mappings, &GlobalSpecification::transitions_e3, "Transitions/E3/", ej_specs, "E3 transition operator: ");
    add_delegates_for(global_specifications_mappings, &GlobalSpecification::transitions_e4, "Transitions/E4/", ej_specs, "E4 transition operator: ");
    add_delegates_for(global_specifications_mappings, &GlobalSpecification::transitions_m1, "Transitions/M1/", ej_specs, "M1 transition operator: ");
    add_delegates_for(global_specifications_mappings, &GlobalSpecification::transitions_m2, "Transitions/M2/", ej_specs, "M2 transition operator: ");
    add_delegates_for(global_specifications_mappings, &GlobalSpecification::transitions_m3, "Transitions/M3/", ej_specs, "M3 transition operator: ");
    add_delegates_for(global_specifications_mappings, &GlobalSpecification::transitions_m4, "Transitions/M4/", ej_specs, "M4 transition operator: ");

    // And the rest...
    add_delegates_for(global_specifications_mappings, &GlobalSpecification::transitions_hfs1, "Transitions/HFS1/", hfs1_specs , "HFS dipole transitions: ");
    add_delegates_for(global_specifications_mappings, &GlobalSpecification::transitions_hfs2, "Transitions/HFS2/", hfs2_specs , "HFS quadrupole transitions: ");
    add_delegates_for(global_specifications_mappings, &GlobalSpecification::transitions_hfi, "Transitions/HFI/", hfi_specs , "Generalised HFS transitions: ");
    add_delegates_for(global_specifications_mappings, &GlobalSpecification::transitions_fs, "Transitions/FS/", fs_specs , "Field shift transitions: ");
    add_delegates_for(global_specifications_mappings, &GlobalSpecification::transitions_qed, "Transitions/QED/", qed_specs , "QED radiative potential transitions: ");
    add_delegates_for(global_specifications_mappings, &GlobalSpecification::transitions_yukawa, "Transitions/Yukawa/", yukawa_specs , "Yukawa potential transitions: ");
    add_delegates_for(global_specifications_mappings, &GlobalSpecification::transitions_ke, "Transitions/KE/", ke_specs , "Kinetic energy transitions: ");
    add_delegates_for(global_specifications_mappings, &GlobalSpecification::transitions_llit2, "Transitions/LLIT2/", llit2_specs , "Local-Lorentz Invariance T2 transitions: ");
    add_delegates_for(global_specifications_mappings, &GlobalSpecification::transitions_nms, "Transitions/NMS/", nms_specs , "Nuclear mass-shift (NMS) transitions: ");

    // Autoionization
    add_delegates_for(global_specifications_mappings, &GlobalSpecification::ic, "IC/", ic_specs, "Internal Conversion (IC):");
    add_delegates_for(global_specifications_mappings, &GlobalSpecification::residual_coulomb, "Transitions/ResidualCoulomb/", residual_coulomb_specs, "ResidualCoulomb:");

    // Deprecated keywords: registered at the top level (empty prefix) so the parser accepts them,
    // then rejected in validateAndNormaliseSpecification() with a pointer to the replacement.
    auto deprecated_specs = make_deprecated_specification<DeprecatedSpecification>();
    add_delegates_for(global_specifications_mappings, &GlobalSpecification::deprecated, "", deprecated_specs, "DEPRECATED: ");
}

// Note that this doesn't do any parsing, it just returns a view of the GlobalSpecification, with a
// nice key-value interface for usability. Therefore, it has to be called *after* we've done all
// the parsing and validation
SpecificationMap get_specification_map_view(GlobalSpecification& gs)
{
    keyed_record_view<GlobalSpecification> specification_map(gs.global_specifications_mappings, gs);
    return(specification_map);
}

// Print the specification key-value pairs to "outstream"
// Modified from export_ini in parapara, but modified to better reflect AMBiT's I/O needs
P::hopefully<void> GlobalSpecification::PrintSpecification()
{
    //P::export_ini(*this, global_specifications_mappings, P::default_writer(), *outstream);

    auto record = *this;
    auto specs = global_specifications_mappings;
    auto wtr = P::default_writer();
    std::string secsep = "/";

    constexpr auto npos = std::string::npos;
    std::map<std::string, std::stringstream> section_content;

    for (const auto& spec: specs) {
        std::string section;
        std::string_view key(spec.key);

        if (auto j = key.find(secsep); j!=npos) {
            section = key.substr(0, j);
            key.remove_prefix(j+1);
        }

        std::stringstream& content = section_content[section];
        if (content.tellp()) content << '\n';

        if (auto hs = spec.write(record, wtr)) {
            // Don't print the value if it's a "--" flag: the fact that it's present
            // implies it's true
            auto minus_idx = key.find("--");
            if(minus_idx != std::string::npos)
            {
                content << key;
            }
            else
            {
                content << key << " = " << hs.value();
            }
        }
        else if (hs.error().error == P::failure::empty_optional) {
            continue;
        }
        else {
            return P::unexpected(hs.error());
        }
    }

    //if (section_content.count("")) {
    //    *outstream << section_content[""].rdbuf(); //<< '\n';
    //}

    for (const auto& [section, content]: section_content) {
        if (section.empty()) continue;

        *outstream << '[' << section << "]\n\n" << content.rdbuf() << '\n';
    }

    return {};
}

/**********************************************************************
 * Actual parsing happens here
 **********************************************************************/
std::string importSpecificationFile(GlobalSpecification* gs, const std::string& fileName) {
    std::ifstream in(fileName);
    if (!in) return "unable to open input file '"+fileName+"'";

    P::specification_map<GlobalSpecification> global_specifications_dict(gs->global_specifications_mappings, P::keys_lc_nows);
    auto h = custom_import_ini(*gs, global_specifications_dict, in);
    if (!h) return P::explain(h.error(), true);

    return "";
}

// Read a full INI from a string, rather than a file. Mostly useful for unit testing
std::string importSpecificationINIString(GlobalSpecification* gs, const std::string& ini_string) {
    std::istringstream in(ini_string);
    P::specification_map<GlobalSpecification> global_specifications_dict(gs->global_specifications_mappings, P::keys_lc_nows);
    auto h = custom_import_ini(*gs, global_specifications_dict, in);
    if (!h) return P::explain(h.error(), true);

    return "";
}

std::string importSpecificationKV(GlobalSpecification* gs, const std::string& assignment) {
    *outstream << assignment << '\n';
    P::specification_map<GlobalSpecification> global_specifications_dict(gs->global_specifications_mappings, P::keys_lc_nows);
    auto h = P::import_k_eq_v(*gs, global_specifications_dict, assignment);
    if (!h) return P::explain(h.error(), true);
    else return "";
}

/**************************** Validation *******************************/
// Perform global validation of specifications. Return non-empty error message on failure.
std::string GlobalSpecification::validateAndNormaliseSpecification()
{
    // Reject deprecated keywords with a message pointing at the replacement.
    if(deprecated.multirun)
        return "'Multirun' is no longer supported; see the 'Multiple runs' section of the user guide.";
    if(deprecated.r_flag)
        return "The '-r' option is no longer supported; see the 'Multiple runs' section of the user guide.";
    if(deprecated.nuclear_inverse_mass)
        return "'NuclearInverseMass' at the top level is no longer supported; "
               "use 'HF/MassShift/NuclearInverseMass' instead.";
    if(deprecated.hf_sms)
        return "'HF/--sms' is no longer supported; use 'HF/MassShift/--sms' instead.";
    if(deprecated.hf_nms)
        return "'HF/--nms' is no longer supported; use 'HF/MassShift/--nms' instead.";
    if(deprecated.hf_only_rel_nms)
        return "'HF/--only-relativistic-nms' is no longer supported; "
               "use 'HF/MassShift/--only-relativistic-nms' instead.";
    if(deprecated.hf_nonrel_mass_shift)
        return "'HF/--nonrelativistic-mass-shift' is no longer supported; "
               "use 'HF/MassShift/--nonrelativistic-mass-shift' instead.";
    if(deprecated.hf_include_lower_sms)
        return "'HF/--include-lower-sms' is no longer supported; "
               "use 'HF/MassShift/--include-lower-sms' instead.";
    if(deprecated.hf_charge)
        return "'HF/Charge' is no longer supported; specify 'HF/N' instead "
               "(charge is derived as Z - HF/N).";
    if(deprecated.basis_size)
        return "'Basis/BasisSize' is no longer supported; "
               "use 'Basis/ValenceBasis' (and 'MBPT/Basis' for the MBPT basis) instead.";
    if(deprecated.basis_custom_basis)
        return "'Basis/--custom-basis' is no longer supported; "
               "use one of 'Basis/--hf-basis', 'Basis/--bspline-basis', or 'Basis/--xr-basis' instead.";
    if(deprecated.basis_r_basis)
        return "'Basis/--r-basis' is no longer supported; "
               "use one of 'Basis/--hf-basis', 'Basis/--bspline-basis', or 'Basis/--xr-basis' instead.";
    if(deprecated.basis_custom_orbitals)
        return "'Basis/CustomOrbitals' is no longer supported; use 'Basis/xR/CustomOrbitals' instead.";
    if(deprecated.mbpt_use_subtraction)
        return "'MBPT/--use-subtraction' is no longer supported; use 'MBPT/--force-subtraction' instead.";
    if(deprecated.mbpt_brueckner_excited)
        return "'MBPT/Brueckner/--excited' is no longer supported; use 'MBPT/--mbpt-use-brueckner' instead.";
    if(deprecated.ci_output_no_configs)
        return "'CI/Output/--no-configs' is no longer supported; use 'CI/Output/--no-configurations' instead.";
    if(deprecated.ci_output_min_displayed_percent)
        return "'CI/Output/MinimumDisplayedPercentage' is no longer supported; "
               "use 'CI/Output/MinDisplayedPercentage' instead.";

    // Check for critical stuff
    if(ID.empty())
        return "ID must be specified.";

    if(Z == 0)
        return "Z must be specified.";

    // Check for consistent Lattice settings:
    if(lattice_exponential && lattice_end_point>0)
        return "Lattice/EndPoint is ignored if Lattice/--exp-lattice is set. Use Lattice/H instead.";

    // Check for consistency in HF decorators (e.g. cannot request both nonrelativistic and
    // relativistic-only constraints)
    if (hf_ms_only_rel_nms && hf_ms_nonrel_mass_shift)
        return "HF/--only-relativistic-nms and HF/--nonrelativistic-mass-shift cannot both be set at the same time";
    // Must have at least one of HF/N and HF/Charge, so error out if these are missing
    if(!hf_N)
        return "Must specify HF/N, which must correspond to the electron count in HF/Configuration.";
    if(!hf_configuration)
        return "HF/Configuration must be set (can be set to '' or 0 for no-core case).";

    // Yukawa interaction has three possible parameters to determine its mass: Mass, MassEV and
    // rc, which are logically equivalent but differ in units. The input file must have exactly one
    // of these defined.
    //
    // Count the number of parameters which are set. There's only three so a dumb solution is
    // okay
    if(hf_do_yukawa)
    {
        unsigned count = 0;
        if(hf_yukawa_mass)
            count++;
        if(hf_yukawa_massev)
            count++;
        if(hf_yukawa_rc)
            count++;

        if(count != 1)
            return "Must specify exactly one of HF/Yukawa/Mass, HF/Yukawa/MassEV, or HF/Yukawa/Rc when using HF/Yukawa";
    }

    // One type of basis at a time
    unsigned basis_type_count = 0;
    if(basis_bspline)
        basis_type_count++;
    if(basis_xr)
        basis_type_count++;
    if(basis_hf)
        basis_type_count++;

    if(basis_type_count > 1)
        return "Must specify at most one of Basis/--bspline-basis, Basis/--xr-basis, or Basis/--hf-basis";

    // Default to bsplines if no options specified
    // Can't do this as a default value in the specification struct because there's no way
    // to "turn off" a basis type (so if we default to basis_bspline = true, then it will
    // *always* be true)
    if(basis_type_count == 0)
        basis_bspline = true;

    // Set the general MBPT (one|two|three)_body_integral variables
    if(_s1 || _s12 || _s13 || _s123)
    {
        mbpt_one_body = true;
    }
    if(_s2 || _s12 || _s23 || _s123)
    {
        mbpt_two_body = true;
    }
    if(_s3 || _s13 || _s23 || _s123)
    {
        mbpt_three_body = true;
    }

    // We must have at least one leading configuration for CI and CI/SmallSide
    if(!no_ci && do_ci && ci.leading_configurations.empty())
        return "Must have at least one config in CI/LeadingConfigurations";

    if(ci_do_smallside && ci_smallside.leading_configurations.empty())
        return "Must have at least one config in CI/SmallSide/LeadingConfigurations";

    // Must specify LevelDirectory if using CI/--memory-saver
    if(ci.memory_saver && !level_directory)
        return "Must specify LevelDirectory when using CI/--memory-saver";
    // Can't request both CI/--gfactors and CI/--no-gfactors
    if (ci.gfactors && ci.no_gfactors)
        return "CI/--gfactors and CI/--no-gfactors cannot both be set at the same time";

    // The number of electron/hole excitation bounds must be equal to the number of
    // electron/hole excitations
    if(ci.electron_bounds)
    {
        if(ci.electron_bounds.value().size() != ci.electron_excitations)
        {
            return "Incorrect number of arguments to CI/ElectronBounds";
        }
    }
    if(ci.hole_bounds)
    {
        if(ci.hole_bounds.value().size() != ci.hole_excitations)
        {
            return "Incorrect number of arguments to CI/HoleBounds";
        }
    }

    // Same for SmallSide
    if(ci_smallside.electron_bounds)
    {
        if(ci_smallside.electron_bounds.value().size() != ci_smallside.electron_excitations)
        {
            return "Incorrect number of arguments to CI/SmallSide/ElectronBounds";
        }
    }
    if(ci_smallside.hole_bounds)
    {
        if(ci_smallside.hole_bounds.value().size() != ci_smallside.hole_excitations)
        {
            return "Incorrect number of arguments to CI/SmallSide/HoleBounds";
        }
    }

    // BSpline order K must be strictly less than basis size N
    if(basis_bspline && basis_bspline_k >= basis_bspline_N)
        return "Basis/BSpline/K must be less than Basis/BSpline/N";

    // Oll Korrect
    return "";
}

/**********************************************************************
 * Custom INI import implementation:
 *********************************************************************/
using P::ini_record;
using P::ini_record_kind;

ini_record custom_ini_parser(std::string_view v) {
    using token = ini_record::token;
    using size_type = std::string_view::size_type;

    constexpr size_type npos = std::string_view::npos;
    constexpr std::string_view ws{" \t\f\v\r\n"};

    auto comment_at = [](std::string_view v, size_type p) { return v.substr(p, 2)=="//"; };
    size_type b = v.find_first_not_of(ws);

    // empty or comment?
    if (b==npos || comment_at(v, b)) return ini_record{ini_record_kind::empty};

    // section heading?
    if (v[b]=='[') {
        size_type e = v.find(']');

        // check for malformed heading
        if (e==npos) return {ini_record_kind::syntax_error, token("", b+1)};

        if (e+1<v.length()) {
            auto epilogue = v.find_first_not_of(ws, e+1);
            if (epilogue!=npos && !comment_at(v, epilogue)) {
                return {ini_record_kind::syntax_error, token("", epilogue)};
            }
        }

        b = v.find_first_not_of(ws, b+1);
        e = v.find_last_not_of(ws, e-1);

        return {ini_record_kind::section, token(v.substr(b, e>=b? e+1-b: 0), b+1)};
    }

    // expect key first, followed by ws and eol, =, or //.
    size_type j = std::min(v.find('=', b), v.find("//", b));
    token key_token{v.substr(b, j==b? 0: v.find_last_not_of(ws, j-1)+1-b), b+1};

    // key without value?
    if (j==npos || v[j]!='=') {
        return {ini_record_kind::key, key_token};
    }

    // skip to text after =, look for value
    size_type eq = j;
    size_type value_cindex = eq;

    if (j<v.length()) {
        j = v.find_first_not_of(ws, j+1);
        if (j!=npos && !comment_at(v, j)) {
            value_cindex = j+1;

            // if value is not quoted, take text up to eol or first eol, discarding trailing ws
            if (v[j]!='\'') {
                size_type end = v.find("//", j);
                if (end!=npos) --end;

                return {ini_record_kind::key_value, key_token,
                        token{v.substr(j, v.find_last_not_of(ws, end)-j+1), value_cindex}};
            }
            else {
                // quoted value; take until next unescaped '
                std::string value;
                size_type epilogue = npos;
                bool esc = false;

                for (size_type i = j+1; i<v.length(); ++i) {
                    if (esc) {
                        value += v[i];
                        esc = false;
                    }
                    else if (v[i]=='\'') {
                        epilogue = i+1;
                        break;
                    }
                    else if (v[i]=='\\') {
                        esc = true;
                    }
                    else {
                        value += v[i];
                    }
                }

                // unterminated quoted value or escaped eol?
                if (epilogue==npos || esc) {
                    return {ini_record_kind::syntax_error, token("", j)};
                }

                // extra stuff following value that is not a comment?
                if (epilogue<v.length()) {
                    epilogue = v.find_first_not_of(ws, epilogue);
                    if (epilogue!=npos && !comment_at(v, epilogue)) {
                        return {ini_record_kind::syntax_error, token("", epilogue)};
                    }
                }

                return {ini_record_kind::key_value, key_token, token{value, value_cindex}};
            }
        }
    }
    // key with empty value
    return {ini_record_kind::key_value, key_token, token{"", eq}};
}

// Use a custom line-by-line ini importer to handle relative section headings
template <typename Record>
P::hopefully<void> custom_import_ini(Record& rec, const P::specification_map<Record>& specs, std::istream& in)
{
    constexpr auto npos = std::string_view::npos;

    P::ini_style_importer importer(custom_ini_parser, in);
    while (importer) {
        std::string prev_sec{importer.section()};
        auto h = importer.run_one(rec, specs, P::default_reader(), "/");

        if (!h)
            return P::unexpected(h.error());

        if (h.value() == P::ini_record_kind::eof)
            break;

        if (h.value() != P::ini_record_kind::section)
        {   // Print key and value
            if(!prev_sec.empty())
                *outstream << prev_sec << '/';
            *outstream << importer.context().record << '\n';
            continue;
        }
        else
            *outstream << importer.context().record << '\n';

        std::string_view new_sec = importer.section();
        if (new_sec.substr(0, 2)=="./") {
            std::string section{prev_sec};
            section += new_sec.substr(1);
            importer.section(section);
        }
        else if (new_sec.substr(0, 3)=="../") {
            std::string_view prefix(prev_sec);
            do {
                new_sec.remove_prefix(3);
                if (auto tail = prefix.rfind('/'); tail != npos) {
                    prefix = prefix.substr(0, tail);
                }
            } while (new_sec.substr(0, 3)=="../");

            std::string section{prefix};
            section += "/";
            section += new_sec;
            importer.section(section);
        }
    }
    return {};
}

} // namespace Ambit

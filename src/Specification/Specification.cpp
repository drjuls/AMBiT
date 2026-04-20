#include <algorithm>
#include <fstream>
#include <istream>
#include <string>
#include <string_view>
#include <vector>

#include "Specification.h"
#include "parapara/parapara.h"

namespace P = parapara;

namespace Ambit
{

// Use a custom line-by-line ini importer to handle relative section headings
template <typename Record>
P::hopefully<void> custom_import_ini(Record& rec, const P::specification_map<Record>& specs, std::istream& in);

inline auto positive = P::greater_than(0, "must be positive");

/********************* GlobalSpecification methods *********************/
GlobalSpecification* GlobalSpecification::Instance()
{
    static GlobalSpecification instance;
    return &instance;
}

GlobalSpecification::GlobalSpecification()
{
    // Create the mappings from CLI arguments to GlobalSpecification fields
    MakeGlobalSpecificationsMappings();
}

/********************* Parsing methods *********************/
// Helper methods for the polymorphic "delegated" parsing (so we can e.g. re-use common options
// between [Transition] operators. Much of this is derived from "ex-delegate.cc" in the parapara
// examples: https://github.com/halfflat/parapara/blob/main/ex/ex-delegate.cc

// A bit of syntactic sugar for inheritance (even if it does involve some copies):
template <typename C>
C extend(C&& container, const std::initializer_list<typename C::value_type>& il) {
    std::copy(std::begin(il), std::end(il), std::back_inserter(container));
    return container;
}

template <typename X>
std::vector<P::specification<X>> GlobalSpecification::make_ci_base_specification() {
    return {
        {"LeadingConfigurations",                &X::leading_configurations},
        {"LeadingRelativisticConfigurations",    &X::leading_rel_configurations},
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
        {"Output/--no-configs",                             &X::output_no_configs}, 
        {"Output/MaxDisplayedEnergy",                       &X::output_max_displayed_energy}, 
        {"Output/MinDisplayedPercent",                      &X::output_min_displayed_percent}, 
        {"Output/Separator",                                &X::output_separator}, 
     });
}


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
    {"Z",                       &GlobalSpecification::Z, positive},
    {"LevelDirectory",          &GlobalSpecification::level_directory},
    // NOTE: Don't check "GlobalSpecification::s1" and friends in the main body of the code, use
    // "GlobalSpecification::mbpt_one_body" etc instead; it's much cleaner to figure out which
    // integrals to include once when normalising the spec
    {"-m",                       &GlobalSpecification::m},
    {"-s1",                      &GlobalSpecification::_s1},
    {"-s2",                      &GlobalSpecification::_s2},
    {"-s3",                      &GlobalSpecification::_s3},
    {"-s12",                     &GlobalSpecification::_s12},
    {"-s13",                     &GlobalSpecification::_s13},
    {"-s23",                     &GlobalSpecification::_s23},
    {"-s123",                   &GlobalSpecification::_s123},
    {"--no-new-mbpt",           &GlobalSpecification::no_new_mbpt},
    {"--check-sizes",           &GlobalSpecification::check_sizes},
    {"-c",                      &GlobalSpecification::clean_run},
    {"-p",                      &GlobalSpecification::print_basis},
    {"--ci-complete",           &GlobalSpecification::ci_complete},
    {"--no-ci",                 &GlobalSpecification::no_ci},
    {"--configuration-average", &GlobalSpecification::configuration_average},
    {"NuclearRadius",           &GlobalSpecification::nuclear_radius},
    {"NuclearThickness",        &GlobalSpecification::nuclear_thickness},
    {"NuclearInverseMass",      &GlobalSpecification::nuclear_inverse_mass},
    {"AlphaSquaredVariation",   &GlobalSpecification::alpha_squared_variation},
    // Lattice
    {"Lattice/NumPoints",           &GlobalSpecification::lattice_num_points,    P::nonzero()},
    {"Lattice/StartPoint",          &GlobalSpecification::lattice_start_point,   positive},
    {"Lattice/EndPoint",            &GlobalSpecification::lattice_end_point,     positive},
    {"Lattice/--exp-lattice",       &GlobalSpecification::lattice_exponential},
    {"Lattice/H",                   &GlobalSpecification::lattice_H,             positive},
    // HF
    {"HF/N",                            &GlobalSpecification::hf_N, positive},
    {"HF/Charge",                       &GlobalSpecification::hf_charge},
    {"HF/Configuration",                &GlobalSpecification::hf_configuration},
    {"HF/--breit",                      &GlobalSpecification::hf_breit},
    {"HF/--sms",                        &GlobalSpecification::hf_sms},
    {"HF/--nms",                        &GlobalSpecification::hf_nms},
    {"HF/--only-relativistic-nms",      &GlobalSpecification::hf_only_rel_nms},
    {"HF/--nonrelativistic-mass-shift", &GlobalSpecification::hf_nonrel_mass_shift},
    {"HF/--include-lower-mass",         &GlobalSpecification::hf_include_lower_mass},
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
    // HF/NuclearPolarisability
    {"HF/NuclearPolarisability",     &GlobalSpecification::hf_do_nuclear_polarisability},
    {"HF/NuclearPolarisability/AlphaE",     &GlobalSpecification::hf_nuclear_polarisability_alpha_e, P::nonzero()},
    {"HF/NuclearPolarisability/EbarMeV",    &GlobalSpecification::hf_nuclear_polarisability_ebar_mev, P::nonzero()},
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
    {"Basis/BasisSize",         &GlobalSpecification::basis_size},
    {"Basis/ValenceBasis",      &GlobalSpecification::basis_valence},
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
    {"Basis/XR/CustomOrbitals",    &GlobalSpecification::basis_xr_custom_orbitals},
    // Basis/BSpline
    {"Basis/BSpline/Rmax",          &GlobalSpecification::basis_bspline_rmax},
    {"Basis/BSpline/R0",            &GlobalSpecification::basis_bspline_r0},
    {"Basis/BSpline/K",             &GlobalSpecification::basis_bspline_k},
    {"Basis/BSpline/N",             &GlobalSpecification::basis_bspline_N},
    {"Basis/BSpline/SplineType",    &GlobalSpecification::basis_bspline_splinetype},
    // CI
// CI/SmallSide
    {"CI/SmallSide",    &GlobalSpecification::ci_do_smallside},
    // MBPT
    {"MBPT/Basis",                             &GlobalSpecification::mbpt_basis}, 
    {"MBPT/EnergyDenomOrbitals",               &GlobalSpecification::mbpt_energy_denom_orbitals}, 
    {"MBPT/--use-valence",                     &GlobalSpecification::mbpt_use_valence}, 
    {"MBPT/--no-core",                         &GlobalSpecification::mbpt_no_core}, 
    {"MBPT/--use-subtraction",                 &GlobalSpecification::mbpt_use_subtraction}, 
    {"MBPT/--no-subtraction",                  &GlobalSpecification::mbpt_no_subtraction}, 
    {"MBPT/--no-extra-box",                    &GlobalSpecification::mbpt_no_extra_box}, 
    {"MBPT/EnergyDenomFloor",                 &GlobalSpecification::mbpt_energy_denom_floor}, 
    {"MBPT/Delta",                            &GlobalSpecification::mbpt_delta}, 
    {"MBPT/TwoBody/StorageLimits",            &GlobalSpecification::mbpt_twobody_storage_limits}, 
    {"MBPT/OneBody/Scaling",                  &GlobalSpecification::mbpt_onebody_scaling}, 
    {"MBPT/--brueckner",                      &GlobalSpecification::mbpt_brueckner}
};

    // CI and CI/SmallSide specs
    auto ci_specs = make_ci_specification<CISpecification>();
    auto ci_smallside_specs = make_ci_base_specification<CIBaseSpecification>();

    // Delegates
    add_delegates_for(global_specifications_mappings, &GlobalSpecification::ci, "CI/", ci_specs, "CI full: ");
    add_delegates_for(global_specifications_mappings, &GlobalSpecification::ci_smallside, "CI/SmallSide/", ci_smallside_specs, "CI SmallSide: ");
}

// Note that this doesn't do any parsing, it just returns a view of the GlobalSpecification, with a
// nice key-value interface for usability. Therefore, it has to be called *after* we've done all
// the parsing and validation
SpecificationMap get_config_map_view(GlobalSpecification& gs)
{
    keyed_record_view<GlobalSpecification> config_map(gs.global_specifications_mappings, gs);
    return(config_map);
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

std::string importSpecificationKV(GlobalSpecification* gs, const std::string& assignment) {
    P::specification_map<GlobalSpecification> global_specifications_dict(gs->global_specifications_mappings, P::keys_lc_nows);
    auto h = P::import_k_eq_v(*gs, global_specifications_dict, assignment);
    if (!h) return P::explain(h.error(), true);
    else return "";
}

/**************************** Validation *******************************/
// Perform global validation of specifications. Return non-empty error message on failure.
std::string GlobalSpecification::validateAndNormaliseSpecification() {
    // Check for consistent Lattice settings:

    if (lattice_exponential && lattice_end_point>0)
        return "Lattice/EndPoint is ignored if Lattice/--exp-lattice is set";

    // Check for consistency in HF decorators (e.g. cannot request both nonrelativistic and
    // relativistic-only constraints)
    if (hf_only_rel_nms && hf_nonrel_mass_shift)
        return "HF/--only-relativistic-nms and HF/--nonrelativistic-mass-shift cannot both be set at the same time";
    // Must have at least one of HF/N and HF/Charge, so error out if these are missing
    if((!hf_charge) && (!hf_N))
        return "Must specify at least one of HF/N or HF/Charge";
    
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
    if(ci.leading_configurations.empty())
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

        // Ignore unrecognized keys until specification coverage is complete
        if (!h) {
             if (h.error().error!=P::failure::unrecognized_key) return P::unexpected(h.error());
             else continue; 
        }

        if (h.value() != P::ini_record_kind::section) continue;

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

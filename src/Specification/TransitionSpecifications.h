#ifndef TRANSITIONS_SPECIFICATIONS_H
#define TRANSITIONS_SPECIFICATIONS_H

#include <optional>
#include <vector>
#include <string>

namespace Ambit
{

/************************** Transitions ***********************/
// Base class from which all transition operators define their CLI spec
class TransitionsBaseSpecification
{
public:
    std::optional<std::vector<std::string>> matrix_elements;
    std::optional<double> all_below;
    bool reduced_elements = false;
    bool print_integrals = false;

    double scale = 1.0;

    // RPA options
    bool do_rpa = false;
    std::optional<double> rpa_scale;
    std::optional<double> rpa_weighting;
    std::optional<int> rpa_max_delta_k;
    int rpa_max_iterations = 100;
    std::string rpa_configuration;
    bool rpa_no_negative_states = false;
    bool rpa_print_field = false;
    int rpa_bspline_n = 40;
    int rpa_bspline_k = 7;
    std::optional<double> rpa_bspline_rmax;
    double rpa_bspline_r0 = 0.0;
};

class TransitionsEJSpecification: public TransitionsBaseSpecification
{
public:
    double frequency = 0.1;
    bool velocity_gauge = false;
};

class TransitionsHFS1Specification: public TransitionsBaseSpecification
{
public:
    std::optional<double> nuclear_magnetic_radius;
    double g_on_i = 1.0;
};

class TransitionsHFS2Specification: public TransitionsBaseSpecification
{
public:
    std::optional<double> nuclear_quadrupole_radius;
    double q = 1.0;
};

class TransitionsHFSISpecification: public TransitionsBaseSpecification
{
public:
    std::optional<double> nuclear_radius;
    std::optional<std::string> op; // Needs to be sanitized/validated to have the form
                                   // "EJ" or "MJ" (where J is some multipolarity)
};

class TransitionsFSSpecification: public TransitionsBaseSpecification
{
public:
    double delta_nuclear_radius = 0.1;
    double delta_nuclear_thickness = 0.0;
};

class TransitionsQEDSpecification: public TransitionsBaseSpecification
{
public:
    bool uehling = false;
    bool self_energy = false;
    bool use_nuclear_density = false;
    std::optional<double> nuclear_rms_radius;
    bool no_magnetic = false;
    bool no_electric = false;
    bool skip_offmass = false;
    bool use_electron_screening = false;
};

class TransitionsYukawaSpecification: public TransitionsBaseSpecification
{
public:
    std::optional<double> mass;
    std::optional<double> massev;
    std::optional<double> rc;
};

class TransitionsKESpecification: public TransitionsBaseSpecification
{};

class TransitionsLLIT2Specification: public TransitionsBaseSpecification
{};

class TransitionsNMSSpecification: public TransitionsBaseSpecification
{
public:
    bool only_relativistic_nms = false;
    bool nonrelativistic_mass_shift = false;
};

// Atom_Autoionization stuff
// IC options should be a subclass of the HFSI specification, since it gets used to
// initialise a GeneralisedHyperfine calculator
class ICSpecification : public TransitionsHFSISpecification
{
public:
    std::string source_filename;
    std::optional<int> source_twoj;
    std::string source_parity;
    int source_index;
    double nuclear_energy_ev;
    unsigned continuum_l_min = 0;
    unsigned continuum_l_max = 6;
    std::optional<std::string> residue;
};

class ResidualCoulombSpecification: public TransitionsBaseSpecification
{};

}
#endif

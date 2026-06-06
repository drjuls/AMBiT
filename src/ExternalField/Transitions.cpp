#include "Transitions.h"
#include "Include.h"
#include "Atom/HamiltonianTypes.h"
#include "HartreeFock/ConfigurationParser.h"
#include "TimeDependentSpinorOperator.h"
#include "ExternalField/RPAOperator.h"
#include "ExternalField/RPASolver.h"

namespace Ambit
{
std::string Name(const LevelID& levelid)
{
    return (levelid.first->Name() + ":" + itoa(levelid.second));
}

void TransitionCalculator::CalculateAndPrint()
{
    // Get all named matrix elements
    int num_transitions = 0;
    if(specification.matrix_elements)
        num_transitions = specification.matrix_elements.value().size();
    bool all_below = specification.all_below.has_value();
    bool print_integrals = specification.print_integrals;

    *outstream << "\n";
    PrintHeader();

    if(print_integrals)
    {   PrintIntegrals();
        *outstream << std::endl;
        PrintHeader();
    }

    for(int i = 0; i < num_transitions; i++)
        CalculateTransition(specification.matrix_elements.value()[i]);

    if(!all_below)
    {
        if(!num_transitions)
            *outstream << "  No transitions requested." << std::endl;

        return;
    }

    // Calculate all transitions of a certain type below a given energy
    double max_energy = specification.all_below.value_or(0.0);
    bool found_one = false;

    auto left_it = levels->begin();
    while(left_it != levels->end())
    {
        auto right_it = left_it;
        while(right_it != levels->end())
        {
            if(TransitionExists((*left_it)->GetSymmetry(), (*right_it)->GetSymmetry()))
            {
                LevelVector left_vec = levels->GetLevels(*left_it);
                for(int i = 0; i < left_vec.NumLevels(); i++)
                {
                    if(left_vec.eigenvalues[i] > max_energy)
                        break;

                    LevelVector right_vec = levels->GetLevels(*right_it);
                    for(int j = 0; j < right_vec.NumLevels(); j++)
                    {
                        if(right_vec.eigenvalues[j] > max_energy)
                            break;

                        found_one = true;
                        CalculateTransition(std::make_pair(*left_it, i), std::make_pair(*right_it, j));
                    }
                }
            }
            right_it++;
        }
        left_it++;
    }

    if(!found_one)
    {   *outstream << "  No transitions below E = " << max_energy << std::endl;
        return;
    }
}

void TransitionCalculator::PrintAll() const
{
    if(matrix_elements.size())
    {
        *outstream << "\n";
        PrintHeader();
        for(auto& pair: matrix_elements)
            PrintTransition(pair.first.first, pair.first.second, pair.second);
    }
}

void TransitionCalculator::PrintIntegrals()
{
    if(integrals == nullptr)
    {
        // Create new TransitionIntegrals object and calculate integrals
        integrals = std::make_shared<TransitionIntegrals>(orbitals, op);
        integrals->CalculateOneElectronIntegrals(orbitals->valence, orbitals->valence);
    }

    *outstream << "One-body transition reduced matrix elements (a.u.):" << std::endl;
    *outstream << std::setprecision(8);
    auto it1 = orbitals->valence->begin();
    while(it1 != orbitals->valence->end())
    {
        const auto& orb1 = it1->first;

        auto it2 = it1;
        while(it2 != orbitals->valence->end())
        {
            const auto& orb2 = it2->first;

            if(op->IsNonZero(orb1, orb2))
                *outstream << "  " << orb1.Name() << " -> " << orb2.Name()
                           << " = " << integrals->GetReducedMatrixElement(orb1, orb2) << std::endl;

            ++it2;
        }
        ++it1;
    }
}

pRPAOperator TransitionCalculator::MakeRPA(pSpinorOperator external, pHFOperatorConst hf, pHartreeY hartreeY)
{
    // Get BSplineBasis parameters and make RPA solver
    int N = specification.rpa_bspline_n;
    int K = specification.rpa_bspline_n;
    double Rmax;
    if(specification.rpa_bspline_rmax)
    {
        Rmax = specification.rpa_bspline_rmax.value();
    }
    else
    {
        Rmax = hf->GetLattice()->R(hf->GetCore()->LargestOrbitalSize());
    }

    double dR0 = specification.rpa_bspline_r0;

    // TODO EVK: This uses the default SplineType (Reno), but should probably be upgraded to get
    // this value from Specification. Would need to change the spec, since currently there's no
    // RPA/BSpline/SplineType
    pBSplineBasis basis_maker = std::make_shared<BSplineBasis>(hf->GetLattice(), N, K, Rmax, dR0, SplineType::Reno);

    bool use_negative_states = true;
    if(specification.rpa_no_negative_states)
        use_negative_states = false;

    pRPASolver rpa_solver = std::make_shared<RPASolver>(basis_maker, use_negative_states);

    double propnew;
    if(specification.rpa_weighting)
    {
        propnew = specification.rpa_weighting.value();
        rpa_solver->SetTDHFWeighting(propnew);
    }

    if(specification.rpa_max_delta_k)
    {
        int deltaKlimit = specification.rpa_max_delta_k.value();
        rpa_solver->SetLimitDeltaK(deltaKlimit);
    }

    int max_iterations = specification.rpa_max_iterations;

    // Limit RPA orbitals if requested
    // Make RPA operator
    std::string rpa_configuration = specification.rpa_configuration;
    pRPAOperator rpa;
    if(rpa_configuration.empty())
    {
        rpa = std::make_shared<RPAOperator>(external, hf, hartreeY, rpa_solver);
    }
    else
    {
        OccupationMap rpa_occupancies = ConfigurationParser::ParseFractionalConfiguration(rpa_configuration);
        rpa = std::make_shared<RPAOperator>(external, hf, hartreeY, rpa_occupancies, rpa_solver);
    }

    // Scale operator
    if(!specification.rpa_scale)
    {
        if(rpa->IsStaticRPA())
            scale = 1.0;
        else    // Scale EM transitions by default
            scale = 0.01;
    }
    else
    {
        scale = specification.rpa_scale.value();
    }
    rpa->SetScale(scale);

    DebugOptions.LogHFIterations(true);
    rpa->SolveRPA();

    return rpa;
}

double TransitionCalculator::CalculateTransition(const LevelID& left, const LevelID& right)
{
    double return_value = 0.0;

    // Check transition is allowed
    if(TransitionExists(left, right))
    {
        // Check it doesn't exist already
        TransitionID id = make_transitionID(left, right);

        auto found_it = matrix_elements.find(id);
        if(found_it != matrix_elements.end())
        {
            return_value = found_it->second;
        }
        else
        {   // Get levels
            LevelVector left_levels = levels->GetLevels(left.first);
            if(left_levels.NumLevels() <= left.second)
            {   *errstream << "TransitionCalculator: Level " << Name(left) << " not found." << std::endl;
                return 0.;
            }

            LevelVector right_levels = levels->GetLevels(right.first);
            if(right_levels.NumLevels() <= right.second)
            {   *errstream << "TransitionCalculator: Level " << Name(right) << " not found." << std::endl;
                return 0.;
            }

            if(integrals == nullptr)
            {
                // Create new TransitionIntegrals object and calculate integrals
                integrals = std::make_shared<TransitionIntegrals>(orbitals, op);
                integrals->CalculateOneElectronIntegrals(orbitals->valence, orbitals->valence);
            }

            ManyBodyOperator<pTransitionIntegrals> many_body_operator(integrals);

            // Get matrix elements for all transitions with same HamiltonianIDs
            std::vector<double> values = many_body_operator.GetMatrixElement(left_levels, right_levels);

            // Add to matrix_elements map
            auto value_iterator = values.begin();
            for(int i = 0; i < left_levels.NumLevels(); i++)
            {
                for(int j = 0; j < right_levels.NumLevels(); j++)
                {
                    TransitionID current_id = make_transitionID(std::make_pair(left.first, i), std::make_pair(right.first, j));
                    matrix_elements.insert(std::make_pair(current_id, *value_iterator/scale));
                    value_iterator++;
                }
            }

            return_value = matrix_elements[id];
        }
    }

    // Return transition specifically requested
    PrintTransition(left, right, return_value);
    return return_value;
}

double TransitionCalculator::CalculateTransition(const std::string& transition)
{
    int pos = transition.find("->");
    if(pos == std::string::npos)
    {
        LevelID only(make_LevelID(transition.substr(0, pos)));

        if(!only.first)
        {   *errstream << "TransitionCalculator: " << transition << " incorrectly formed." << std::endl;
            return 0.;
        }
        return CalculateTransition(only, only);
    }
    else
    {   LevelID left(make_LevelID(transition.substr(0, pos)));
        LevelID right(make_LevelID(transition.substr(pos+2)));

        if(!left.first || !right.first)
        {   *errstream << "TransitionCalculator: " << transition << " incorrectly formed." << std::endl;
            return 0.;
        }

        return CalculateTransition(left, right);
    }
}

LevelID TransitionCalculator::make_LevelID(const std::string& name)
{
    LevelID ret(nullptr, 0);

    // Decide on type and then parse
    int index_break = name.find(":");
    int fullstop_break = name.find(".");
    bool noParityFound = (name.find('e') == std::string::npos && name.find('o') == std::string::npos);

    // Single particle type
    if(index_break == std::string::npos && fullstop_break == std::string::npos && noParityFound)
    {
        OrbitalInfo info = ConfigurationParser::ParseOrbital<OrbitalInfo>(name);

        ret.first = std::make_shared<SingleOrbitalID>(info);
        ret.second = 0;
    }
    // Non-relativistic configuration type
    else if(fullstop_break != std::string::npos)
    {
        if(index_break == std::string::npos || fullstop_break > index_break)
            return ret;

        NonRelConfiguration config = ConfigurationParser::ParseConfiguration<NonRelInfo, int>(name.substr(0, fullstop_break));

        // Get symmetry
        int p = fullstop_break+1;
        int twoJ = atoi(name.c_str() + p);
        if((twoJ < 0) || (twoJ > 100))
            return ret;
        while(name[p] && (isdigit(name[p]) || isblank(name[p])))
            p++;

        Parity P;
        if(name[p] == 'o')
            P = Parity::odd;
        else if(name[p] == 'e')
            P = Parity::even;
        else
            return ret;

        int index = atoi(name.c_str() + index_break + 1);
        if(index < 0 || P != config.GetParity())
            return ret;

        ret.first = std::make_shared<NonRelID>(config, twoJ);
        ret.second = index;
    }
    // Usual HamiltonianID type
    else
    {   int twoJ = -1;
        Parity P;
        int index = -1;

        std::stringstream ss(name);
        ss >> twoJ;
        if((twoJ < 0) || (twoJ > 100) || !ss.good())
            return ret;

        char c = ss.get();
        if(c == 'o' && ss.good())
            P = Parity::odd;
        else if(c == 'e' && ss.good())
            P = Parity::even;
        else
            return ret;

        // Colon is optional for usual type
        if(ss.peek() == ':')
            ss.get();

        ss >> index;
        if(index < 0)
            return ret;

        ret.first = std::make_shared<HamiltonianID>(twoJ, P);
        ret.second = index;
    }

    return ret;
}
}

#include "Transitions.h"
#include "ExternalField/EJOperator.h"
#include "Include.h"

std::string Name(const LevelID& levelid)
{
    return (levelid.first->Name() + itoa(levelid.second));
}

LevelID make_LevelID(const std::string& name)
{
    int twoJ;
    Parity P;
    int index;

    std::stringstream ss(name);
    ss >> twoJ;
    if(ss.get() == 'o')
        P = Parity::odd;
    else
        P = Parity::even;

    ss >> index;

    return std::make_pair(std::make_shared<HamiltonianID>(twoJ, P), index);
}

TransitionMap::TransitionMap(Atom& atom, TransitionGauge gauge, TransitionType max):
    atom(atom), preferred_gauge(gauge), max_type(max)
{
    orbitals = atom.GetBasis();
    levels = atom.GetLevels();
}

TransitionID TransitionMap::make_transitionID(const LevelID& left, const LevelID& right, TransitionType type) const
{
    if(left < right)
        return std::make_tuple(left, right, type);
    else
        return std::make_tuple(right, left, type);
}

bool TransitionMap::TransitionExists(const Symmetry& left, const Symmetry& right, TransitionType type) const
{
    int deltaJ = abs(left.GetTwoJ() - right.GetTwoJ())/2;
    if(((left.GetTwoJ() + right.GetTwoJ())/2 < type.second)
       || (deltaJ > type.second))
    {
        return false;
    }

    if((type.first == MultipolarityType::E && type.second%2 == 1)       // E1, E3, ...
       || (type.first == MultipolarityType::M && type.second%2 == 0))   // M2, M4, ...
    {
        if(left.GetParity() == right.GetParity())
        {
            return false;
        }
    }
    else if(left.GetParity() != right.GetParity())
    {
        return false;
    }

    return true;
}

double TransitionMap::CalculateTransition(const LevelID& left, const LevelID& right)
{
    // Get minimum transition type
    int deltaJ = abs(left.first->GetTwoJ() - right.first->GetTwoJ())/2;
    if(deltaJ == 0)
    {
        if(left.first->GetTwoJ() == 0) // 0 -> 0 transition
            return 0.;
        deltaJ++;
    }

    TransitionType type;
    type.second = deltaJ;

    if((deltaJ%2 == 1 && left.first->GetParity() != right.first->GetParity())
        || (deltaJ%2 == 0 && left.first->GetParity() == right.first->GetParity()))
    {
        type.first = MultipolarityType::E;
    }
    else
        type.first = MultipolarityType::M;

    if(max_type < type)
        return 0.;

    return CalculateTransition(left, right, type);
}

double TransitionMap::CalculateTransition(const LevelID& left, const LevelID& right, TransitionType type)
{
    // Check it doesn't exist already
    TransitionID id = make_transitionID(left, right, type);

    if(this->find(id) != this->end())
        return (*this)[id];

    // Check transition is allowed
    if(!TransitionExists(left, right, type))
    {   *errstream << "Transitions::AddTransition(): Transition not allowed: "
                   << Name(left) << " -> " << Name(right) << " (" << type << ")" << std::endl;
        return 0.;
    }

    // Get levels
    auto level_it = levels->find(left.first);
    if(level_it == levels->end() || level_it->second.size() <= left.second)
    {   *errstream << "Transitions::AddTransition(): Level " << Name(left) << " not found." << std::endl;
        return 0.;
    }
    const LevelVector& left_levels = level_it->second;

    level_it = levels->find(right.first);
    if(level_it == levels->end() || level_it->second.size() <= right.second)
    {   *errstream << "Transitions::AddTransition(): Level " << Name(right) << " not found." << std::endl;
        return 0.;
    }
    const LevelVector& right_levels = level_it->second;

    // Get transition integrals
    pTransitionIntegrals transition_integral = integrals[type];
    if(transition_integral == nullptr)
    {
        // Create new TransitionIntegrals object and calculate integrals
        pOPIntegrator integrator(new SimpsonsIntegrator(atom.GetLattice()));
        pSpinorMatrixElementConst pOperator;
        if(type.first == MultipolarityType::E)
            pOperator = std::make_shared<EJOperator>(atom.GetPhysicalConstants(), type.second, integrator, preferred_gauge);
        else
            pOperator = std::make_shared<MJOperator>(atom.GetPhysicalConstants(), type.second, integrator);

        transition_integral.reset(new TransitionIntegrals(orbitals, pOperator));

        integrals[type] = transition_integral;
        transition_integral->CalculateOneElectronIntegrals(orbitals->valence, orbitals->valence);
    }

    ManyBodyOperator<pTransitionIntegrals> many_body_operator(transition_integral);

    // Get matrix elements for all transitions with same HamiltonianIDs
    std::vector<double> values;
    values = many_body_operator.GetMatrixElement(left_levels, right_levels);

    // Convert to strength and add to TransitionMap
    double angular_projection_factor = 1./MathConstant::Instance()->Electron3j(right.first->GetTwoJ(), left.first->GetTwoJ(), type.second, right.first->GetTwoJ(), -left.first->GetTwoJ());

    auto value_iterator = values.begin();
    for(int i = 0; i < left_levels.size(); i++)
    {
        for(int j = 0; j < right_levels.size(); j++)
        {
            double value = (*value_iterator) * angular_projection_factor;
            value = value * value;

            TransitionID current_id = make_transitionID(std::make_pair(left.first, i), std::make_pair(right.first, j), type);
            (*this)[current_id] = value;
            value_iterator++;
        }
    }

    // Get transition specifically requested, print and return
    double value = (*this)[id];
    *outstream << "  " << Name(left) << " -> " << Name(right) << ": S(" << type << ") = " << value << std::endl;

    return value;
}

void TransitionMap::Print() const
{
    *outstream << "\nTransition strengths (S):\n";
    for(auto& pair: *this)
    {
        const TransitionID& id = pair.first;
        *outstream << "  " << Name(std::get<0>(id)) << " -> " << Name(std::get<1>(id))
                   << " (" << std::get<2>(id) << ") = " << pair.second << std::endl;
    }
}

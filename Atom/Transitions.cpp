#include "Transitions.h"
#include "ExternalField/EJOperator.h"
#include "Include.h"

TransitionMap::TransitionMap(Atom& atom, TransitionGauge gauge, TransitionType max):
    atom(atom), preferred_gauge(gauge), max_type(max)
{
    orbitals = atom.GetBasis();
    levels = atom.GetLevels();
}

TransitionID TransitionMap::make_transitionID(const LevelID& left, const LevelID& right, TransitionType type) const
{
    TransitionID id;
    if(left < right)
    {   std::get<0>(id) = left;
        std::get<1>(id) = right;
    }
    else
    {   std::get<0>(id) = right;
        std::get<1>(id) = left;
    }
    std::get<2>(id) = type;

    return id;
}

bool TransitionMap::TransitionExists(const LevelID& left, const LevelID& right, TransitionType type) const
{
    int deltaJ = abs(left.GetTwoJ() - right.GetTwoJ())/2;
    if(((left.GetTwoJ() + right.GetTwoJ())/2 < type.second)
       || (deltaJ > type.second))
    {
        return false;
    }

    if((type.first == MultipolarityType::E && type.second%2 == 1)       // E1, E3, ...
       || (type.first == MultipolarityType::M && type.second%2 == 0))  // M2, M4, ...
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
    int deltaJ = abs(left.GetTwoJ() - right.GetTwoJ())/2;
    if(deltaJ == 0)
    {
        if(left.GetTwoJ() == 0) // 0 -> 0 transition
            return 0.;
        deltaJ++;
    }

    TransitionType type;
    type.second = deltaJ;

    if((deltaJ%2 == 1 && left.GetParity() != right.GetParity())
        || (deltaJ%2 == 0 && left.GetParity() == right.GetParity()))
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
                   << left.Name() << " -> " << right.Name() << " (" << type << ")" << std::endl;
        return 0.;
    }

    // Get levels
    pLevel left_level, right_level;
    auto level_it = levels->find(left);
    if(level_it == levels->end())
    {   *errstream << "Transitions::AddTransition(): Level " << left.Name() << " not found." << std::endl;
        return 0.;
    }
    left_level = level_it->second;

    if(left != right)
    {
        level_it = levels->find(right);
        if(level_it == levels->end())
        {   *errstream << "Transitions::AddTransition(): Level " << right.Name() << " not found." << std::endl;
            return 0.;
        }
        right_level = level_it->second;
    }

    // Get transition integrals
    pTransitionIntegrals transition_integral = integrals[type];
    if(transition_integral == nullptr)
    {
        // Create new TransitionIntegrals object and calculate integrals
        pOPIntegrator integrator(new SimpsonsIntegrator(atom.GetLattice()));
        pSpinorMatrixElementConst pOperator;
        if(type.first == MultipolarityType::E)
            pOperator.reset(new EJOperator(atom.GetPhysicalConstants(), type.second, integrator, preferred_gauge));
        else
            pOperator.reset(new MJOperator(atom.GetPhysicalConstants(), type.second, integrator));

        transition_integral.reset(new TransitionIntegrals(orbitals, pOperator));

        integrals[type] = transition_integral;
        transition_integral->CalculateOneElectronIntegrals(orbitals->valence, orbitals->valence);
    }

    ManyBodyOperator<pTransitionIntegrals> many_body_operator(transition_integral);

    // Get matrix elements for all transitions with same symmetry
    std::vector<double> values;
    Symmetry leftsym = left.GetSymmetry();
    Symmetry rightsym = right.GetSymmetry();
    values = many_body_operator.GetMatrixElement(levels->begin(leftsym), levels->end(leftsym), levels->begin(rightsym), levels->end(rightsym));

    // Convert to strength and add to TransitionMap
    double angular_projection_factor = 1./MathConstant::Instance()->Electron3j(right.GetTwoJ(), left.GetTwoJ(), type.second, right.GetTwoJ(), -left.GetTwoJ());

    auto value_iterator = values.begin();
    for(auto left_it = levels->begin(leftsym); left_it != levels->end(leftsym); left_it++)
    {
        for(auto right_it = levels->begin(rightsym); right_it != levels->end(rightsym); right_it++)
        {
            double value = (*value_iterator) * angular_projection_factor;
            value = value * value;

            TransitionID current_id = make_transitionID(left_it->first, right_it->first, type);
            (*this)[current_id] = value;
            value_iterator++;
        }
    }

    // Get transition specifically requested, print and return
    double value = (*this)[id];
    *outstream << "  " << left.Name() << " -> " << right.Name() << ": S(" << type << ") = " << value << std::endl;

    return value;
}

void TransitionMap::Print() const
{
    *outstream << "\nTransition strengths (S):\n";
    for(auto& pair: *this)
    {
        const TransitionID& id = pair.first;
        *outstream << "  " << std::get<0>(id).Name() << " -> " << std::get<1>(id).Name()
                   << " (" << std::get<2>(id) << ") = " << pair.second << std::endl;
    }
}

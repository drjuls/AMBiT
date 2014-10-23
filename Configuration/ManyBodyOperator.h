#ifndef MANY_BODY_OPERATOR_H
#define MANY_BODY_OPERATOR_H

#include "Level.h"
#include "Projection.h"
#include <tuple>
#include <boost/iterator/counting_iterator.hpp>
#include <boost/iterator/indirect_iterator.hpp>
#include <boost/utility/enable_if.hpp>

class ZeroOperator
{
public:
    template<typename... Args>
    inline double GetMatrixElement(Args&... args) { return 0.; }
};

/** ManyBodyOperator is instantiated with a set of pointers to operators,
    each of which must supply GetMatrixElement(const ElectronInfo& e1, ...).
    The number of operators supplied must be between one and three.
    Ordering of electrons in GetMatrixElement (la->ra, lb->rb, etc):
        One-body:   (la, ra)
        Two-body:   (la, lb, ra, rb)
        Three-body: (la, lb, lc, ra, rb, rc)
 */
template <typename... pElectronOperators>
class ManyBodyOperator
{
public:
    ManyBodyOperator(pElectronOperators... operators): pOperators(operators...)
    {   static_assert(sizeof...(pElectronOperators) < 4, "ManyBodyOperator<> instantiated with too many operators (template arguments).");
    }

    typedef std::vector<const ElectronInfo*> IndirectProjection;
    inline void make_indirect_projection(const Projection& proj, IndirectProjection& indirect_proj) const;

    /** Rearrange p1 and p2 so that differences are at the beginning (up to max differences).
        abs(return value) is number of differences found
        sign(return value) indicates number of permutations (+ve is even number, -ve is odd).
     */
    template<int max_diffs>
    inline int GetProjectionDifferences(IndirectProjection& p1, IndirectProjection& p2) const;

    /** Matrix element between two projections (abstract function). */
    inline double GetMatrixElement(const Projection& proj_left, const Projection& proj_right) const;

    /** Returns <left | O | right> */
    inline double GetMatrixElement(const Level& left, const Level& right) const;

    /** Equivalent to GetMatrixElement(level, level). */
    inline double GetMatrixElement(const Level& level) const;

    /** Equivalent to calculating GetMatrixElement(level) for each level in range [begin, end).
        Return vector of matrix elements.
     */
    inline std::vector<double> GetMatrixElement(LevelMap::const_symmetry_iterator begin, LevelMap::const_symmetry_iterator end) const;

    /** Equivalent to calculating GetMatrixElement(left, right) for each level in their respective ranges [begin, end).
        Return array of matrix elements indexed by left_index * (num_right) + right_index.
     */
    inline std::vector<double> GetMatrixElement(LevelMap::const_symmetry_iterator begin_left, LevelMap::const_symmetry_iterator end_left,
                                                LevelMap::const_symmetry_iterator begin_right, LevelMap::const_symmetry_iterator end_right) const;

protected:
    std::tuple<pElectronOperators...> pOperators;

    // There is always a one-body operator
    inline double OneBodyMatrixElements(const ElectronInfo& la, const ElectronInfo& ra) const
    {
        return std::get<0>(pOperators)->GetMatrixElement(la, ra);
    }

    // Sometimes there is no two-body operator, therefore need to protect std::get<1>().
    // In this version, std::get<1>(pOperators) will always compile
    template<int S = sizeof...(pElectronOperators)>
    inline double TwoBodyMatrixElements(typename boost::enable_if_c< S >= 2, const ElectronInfo>::type& la,
        const ElectronInfo& lb, const ElectronInfo& ra, const ElectronInfo& rb) const
    {
        auto& p_operator = std::get<1>(pOperators);
        return p_operator->GetMatrixElement(la, lb, ra, rb)
                - p_operator->GetMatrixElement(la, lb, rb, ra);
    }

    // This function is never actually called, but must be valid.
    template<int S = sizeof...(pElectronOperators)>
    inline double TwoBodyMatrixElements(typename boost::enable_if_c< S <= 1, const ElectronInfo>::type& la,
        const ElectronInfo& lb, const ElectronInfo& ra, const ElectronInfo& rb) const
    {   return 0.;
    }

    // As with TwoBodyMatrixElements(), we need to protect std::get here.
    template<int S = sizeof...(pElectronOperators)>
    inline double ThreeBodyMatrixElements(typename boost::enable_if_c< S >= 3, const ElectronInfo>::type& la,
        const ElectronInfo& lb, const ElectronInfo& lc, const ElectronInfo& ra, const ElectronInfo& rb, const ElectronInfo& rc) const
    {
        auto& p_operator = std::get<2>(pOperators);
        return p_operator->GetMatrixElement(la, lb, lc, ra, rb, rc)
                + p_operator->GetMatrixElement(la, lb, lc, rb, rc, ra)
                + p_operator->GetMatrixElement(la, lb, lc, rc, ra, rb)
                - p_operator->GetMatrixElement(la, lb, lc, rc, rb, ra)
                - p_operator->GetMatrixElement(la, lb, lc, rb, ra, rc)
                - p_operator->GetMatrixElement(la, lb, lc, ra, rc, rb);
    }

    // This function is never actually called, but must be valid.
    template<int S = sizeof...(pElectronOperators)>
    inline double ThreeBodyMatrixElements(typename boost::enable_if_c< S <= 2, const ElectronInfo>::type& la,
        const ElectronInfo& lb, const ElectronInfo& lc, const ElectronInfo& ra, const ElectronInfo& rb, const ElectronInfo& rc) const
    {   return 0.;
    }
};

/** ManyBodyOperator: one body operator specialization. */
template <typename... pElectronOperators>
double ManyBodyOperator<pElectronOperators...>::GetMatrixElement(const Projection& proj_left, const Projection& proj_right) const
{
    int num_diffs = 0;

    // NB: These are static for speed: to prevent memory (de)allocation
    static ManyBodyOperator<>::IndirectProjection left, right;
    make_indirect_projection(proj_left, left);

    // Skip this for same projection
    if(&proj_left != &proj_right)
    {
        make_indirect_projection(proj_right, right);
        num_diffs = GetProjectionDifferences<sizeof...(pElectronOperators)>(left, right);
    }

    double matrix_element = 0.0;

    switch(sizeof...(pElectronOperators))
    {
        case 1:
            if(num_diffs == 0)
            {
                for(auto& e: left)
                {
                    int sign = e->IsHole()? -1 : 1;
                    matrix_element += OneBodyMatrixElements(*e, *e) * sign;
                }
            }
            else if(abs(num_diffs) == 1)
            {
                matrix_element += OneBodyMatrixElements(*left[0], *right[0]);
            }
            break;
        case 2:
            if(num_diffs == 0)
            {
                auto i = boost::make_indirect_iterator(left.begin());
                auto end = boost::make_indirect_iterator(left.end());
                while(i != end)
                {
                    int sign = i->IsHole()? -1 : 1;
                    matrix_element += OneBodyMatrixElements(*i, *i) * sign;

                    auto j = i;
                    while(j != end)
                    {
                        int total_sign = sign * (j->IsHole()? -1 : 1);
                        matrix_element += TwoBodyMatrixElements(*i, *j, *i, *j) * total_sign;
                        j++;
                    }
                    i++;
                }
            }
            else if(abs(num_diffs) == 1)
            {
                matrix_element += OneBodyMatrixElements(*left[0], *right[0]);

                auto i = boost::make_indirect_iterator(left.begin());
                auto end = boost::make_indirect_iterator(left.end());
                i++;    // Skip first element
                while(i != end)
                {
                    int sign = i->IsHole()? -1 : 1;
                    matrix_element += TwoBodyMatrixElements(*left.front(), *i, *right.front(), *i) * sign;
                    i++;
                }
            }
            else if(abs(num_diffs) == 2)
            {
                matrix_element += TwoBodyMatrixElements(*left[0], *left[1], *right[0], *right[1]);
            }
            break;
        default:
            break;
    }

    if(num_diffs < 0)
        return (-matrix_element);
    else
        return matrix_element;
}

template<typename... pElectronOperators>
void ManyBodyOperator<pElectronOperators...>::make_indirect_projection(const Projection& proj, IndirectProjection& indirect_proj) const
{
    if(indirect_proj.size() != proj.size())
        indirect_proj.resize(proj.size());

    auto indirect_it = indirect_proj.begin();
    for(auto& e: proj)
        *indirect_it++ = &e;
}

template<typename... pElectronOperators>
template<int max_diffs>
int ManyBodyOperator<pElectronOperators...>::GetProjectionDifferences(ManyBodyOperator<>::IndirectProjection& p1, ManyBodyOperator<>::IndirectProjection& p2) const
{
    // NB: These are static for speed: to prevent memory (de)allocation
    static IndirectProjection diff1, diff2;
    diff1.clear();
    diff2.clear();
    static IndirectProjection sorted_p1, sorted_p2;
    sorted_p1.clear();
    sorted_p2.clear();

    auto it1 = p1.begin();
    auto it2 = p2.begin();

    // Number of permutations required to get correct ordering
    int permutations = 0;

    // Cumulative number of operators that are the same on both sides
    int num_same = 0;

    // Electrons and hole operators that have not yet been absorbed into a "difference"
    int num_holes1 = 0;
    int num_electrons1 = 0;
    int num_holes2 = 0;
    int num_electrons2 = 0;

    while(diff1.size() <= max_diffs && diff2.size() <= max_diffs && it1 != p1.end() && it2 != p2.end())
    {
        const ElectronInfo& e1(**it1);
        const ElectronInfo& e2(**it2);

        if(e1 == e2)
        {   sorted_p1.push_back(*it1++);
            sorted_p2.push_back(*it2++);
            num_same++;
        }
        else if(e1 < e2)
        {
            permutations += num_same;
            if(e1.IsHole())
            {
                diff2.push_back(*it1++);

                if(num_electrons1)  // There is an electron to combine already on this side
                {   permutations += num_electrons1;
                    num_electrons1--;
                }
                else if(num_holes2) // Combine with a hole on the other side
                {   num_holes2--;
                }
                else
                    num_holes1++;
            }
            else
            {
                diff1.push_back(*it1++);

                if(num_holes1)  // There is a hole to combine already on this side
                {   permutations += (num_holes1 - 1);
                    num_holes1--;
                }
                else if(num_electrons2) // Combine with electron on the other side
                {   num_electrons2--;
                }
                else
                    num_electrons1++;
            }
        }
        else
        {
            permutations += num_same;
            if(e2.IsHole())
            {
                diff1.push_back(*it2++);

                if(num_electrons2)  // There is an electron to combine already on this side
                {   permutations += num_electrons2;
                    num_electrons2--;
                }
                else if(num_holes1) // Combine with a hole on the other side
                {   num_holes1--;
                }
                else
                    num_holes2++;
            }
            else
            {
                diff2.push_back(*it2++);

                if(num_holes2)  // There is a hole to combine already on this side
                {   permutations += (num_holes2 - 1);
                    num_holes2--;
                }
                else if(num_electrons1) // Combine with electron on the other side
                {   num_electrons1--;
                }
                else
                    num_electrons2++;
            }
        }
    }
    while(it1 != p1.end() && (diff1.size() <= max_diffs))
    {
        permutations += num_same;
        if((*it1)->IsHole())
        {
            diff2.push_back(*it1++);

            if(num_electrons1)  // There is an electron to combine already on this side
            {   permutations += num_electrons1;
                num_electrons1--;
            }
            else if(num_holes2) // Combine with a hole on the other side
            {   num_holes2--;
            }
            else
                num_holes1++;
        }
        else
        {
            diff1.push_back(*it1++);

            if(num_holes1)  // There is a hole to combine already on this side
            {   permutations += (num_holes1 - 1);
                num_holes1--;
            }
            else if(num_electrons2) // Combine with electron on the other side
            {   num_electrons2--;
            }
            else
                num_electrons1++;
        }
    }
    while(it2 != p2.end() && (diff2.size() <= max_diffs))
    {
        permutations += num_same;
        if((*it2)->IsHole())
        {
            diff1.push_back(*it2++);

            if(num_electrons2)  // There is an electron to combine already on this side
            {   permutations += num_electrons2;
                num_electrons2--;
            }
            else if(num_holes1) // Combine with a hole on the other side
            {   num_holes1--;
            }
            else
                num_holes2++;
        }
        else
        {
            diff2.push_back(*it2++);

            if(num_holes2)  // There is a hole to combine already on this side
            {   permutations += (num_holes2 - 1);
                num_holes2--;
            }
            else if(num_electrons1) // Combine with electron on the other side
            {   num_electrons1--;
            }
            else
                num_electrons2++;
        }
    }

    int num_diffs = diff1.size();

    if((diff1.size() > max_diffs) || (diff2.size() > max_diffs))
        return max_diffs+1;

    else if(diff1.empty() && diff2.empty())
        return 0;

    // Additional sign changes for holes; add to permutations
    it1 = diff1.begin();
    it2 = diff2.begin();
    while(it1 != diff1.end())
    {
        if((*it1)->IsHole() || (*it2)->IsHole())
            permutations++;

        it1++;
        it2++;
    }

    //Copy differences
    diff1.insert(diff1.end(), sorted_p1.begin(), sorted_p1.end());
    diff2.insert(diff2.end(), sorted_p2.begin(), sorted_p2.end());

    p1 = diff1;
    p2 = diff2;

    if(permutations%2 == 0)
        return num_diffs;
    else
        return -num_diffs;
}

template<typename... pElectronOperators>
double ManyBodyOperator<pElectronOperators...>::GetMatrixElement(const Level& level) const
{
    double total = 0.;
    const auto& configs = level.GetRelativisticConfigList();
    const auto& eigenvector = level.GetEigenvector();

    auto proj_it = configs->projection_begin();
    while(proj_it != configs->projection_end())
    {
        auto proj_jt = proj_it;

        while(proj_jt != configs->projection_end())
        {
            double matrix_element = GetMatrixElement(*proj_it, *proj_jt);

            if(matrix_element)
            {
                // Summation over CSFs
                double coeff = 0.;

                for(auto coeff_i = proj_it.CSF_begin(); coeff_i != proj_it.CSF_end(); coeff_i++)
                {
                    RelativisticConfigList::const_CSF_iterator start_j = proj_jt.CSF_begin();

                    for(auto coeff_j = start_j; coeff_j != proj_jt.CSF_end(); coeff_j++)
                    {
                        coeff += (*coeff_i) * (*coeff_j)
                                * eigenvector[coeff_i.index()]
                                * eigenvector[coeff_j.index()];
                    }
                }

                // If the projections are different, count twice
                if(proj_it != proj_jt)
                    coeff *= 2.;

                total += coeff * matrix_element;
            }
            proj_jt++;
        }
        proj_it++;
    }

    return total;
}

template<typename... pElectronOperators>
std::vector<double> ManyBodyOperator<pElectronOperators...>::GetMatrixElement(LevelMap::const_symmetry_iterator begin, LevelMap::const_symmetry_iterator end) const
{
    std::vector<const double*> eigenvector;

    for(auto it = begin; it != end; it++)
        eigenvector.push_back(it->second->GetEigenvector());

    if(eigenvector.empty())
        return std::vector<double>();

    pRelativisticConfigListConst configs = begin->second->GetRelativisticConfigList();
    std::vector<double> total(eigenvector.size(), 0.);
    std::vector<double> coeff(eigenvector.size(), 0.);

    unsigned int solution = 0;

    // Iterate over projections
    auto proj_it = configs->projection_begin();
    while(proj_it != configs->projection_end())
    {
        auto proj_jt = proj_it;

        while(proj_jt != configs->projection_end())
        {
            double matrix_element = GetMatrixElement(*proj_it, *proj_jt);

            // coefficients
            if(matrix_element)
            {
                // Summation over CSFs
                for(double& e: coeff)
                    e = 0.;

                for(auto coeff_i = proj_it.CSF_begin(); coeff_i != proj_it.CSF_end(); coeff_i++)
                {
                    RelativisticConfigList::const_CSF_iterator start_j = proj_jt.CSF_begin();

                    for(auto coeff_j = start_j; coeff_j != proj_jt.CSF_end(); coeff_j++)
                    {
                        for(solution = 0; solution < coeff.size(); solution++)
                        {
                            coeff[solution] += (*coeff_i) * (*coeff_j)
                                                * eigenvector[solution][coeff_i.index()]
                                                * eigenvector[solution][coeff_j.index()];
                        }
                    }
                }

                // If the projections are different, count twice
                if(proj_it != proj_jt)
                {   for(double& element: coeff)
                        element *= 2.;
                }
                
                for(solution = 0; solution < coeff.size(); solution++)
                    total[solution] += coeff[solution] * matrix_element;
            }
            proj_jt++;
        }
        proj_it++;
    }

    return total;
}

template<typename... pElectronOperators>
double ManyBodyOperator<pElectronOperators...>::GetMatrixElement(const Level& left, const Level& right) const
{
    double total = 0.;
    const auto& configs_left = left.GetRelativisticConfigList();
    const auto& configs_right = right.GetRelativisticConfigList();
    const auto& eigenvector_left = left.GetEigenvector();
    const auto& eigenvector_right = right.GetEigenvector();

    auto proj_it = configs_left->projection_begin();
    while(proj_it != configs_left->projection_end())
    {
        auto proj_jt = configs_right->projection_begin();
        while(proj_jt != configs_right->projection_end())
        {
            double matrix_element = GetMatrixElement(*proj_it, *proj_jt);

            if(matrix_element)
            {
                // Summation over CSFs
                double coeff = 0.;

                for(auto coeff_i = proj_it.CSF_begin(); coeff_i != proj_it.CSF_end(); coeff_i++)
                {
                    for(auto coeff_j = proj_jt.CSF_begin(); coeff_j != proj_jt.CSF_end(); coeff_j++)
                    {
                        coeff += (*coeff_i) * (*coeff_j)
                                * eigenvector_left[coeff_i.index()]
                                * eigenvector_right[coeff_j.index()];
                    }
                }

                total += coeff * matrix_element;
            }
            proj_jt++;
        }
        proj_it++;
    }
    
    return total;
}

template<typename... pElectronOperators>
std::vector<double> ManyBodyOperator<pElectronOperators...>::GetMatrixElement(LevelMap::const_symmetry_iterator begin_left, LevelMap::const_symmetry_iterator end_left,
                                                                              LevelMap::const_symmetry_iterator begin_right, LevelMap::const_symmetry_iterator end_right) const
{
    std::vector<const double*> left_eigenvector;
    std::vector<const double*> right_eigenvector;
    
    for(auto it = begin_left; it != end_left; it++)
        left_eigenvector.push_back(it->second->GetEigenvector());

    for(auto it = begin_right; it != end_right; it++)
        right_eigenvector.push_back(it->second->GetEigenvector());

    unsigned int return_size = left_eigenvector.size() * right_eigenvector.size();
    if(return_size == 0)
        return std::vector<double>();
    
    pRelativisticConfigListConst configs_left = begin_left->second->GetRelativisticConfigList();
    pRelativisticConfigListConst configs_right = begin_right->second->GetRelativisticConfigList();
    std::vector<double> total(return_size, 0.);
    std::vector<double> coeff(return_size, 0.);
    
    unsigned int solution = 0;
    
    // Iterate over projections
    auto proj_it = configs_left->projection_begin();
    while(proj_it != configs_left->projection_end())
    {
        auto proj_jt = configs_right->projection_begin();
        while(proj_jt != configs_right->projection_end())
        {
            double matrix_element = GetMatrixElement(*proj_it, *proj_jt);
            
            // coefficients
            if(matrix_element)
            {
                // Summation over CSFs
                for(double& e: coeff)
                    e = 0.;
                
                for(auto coeff_i = proj_it.CSF_begin(); coeff_i != proj_it.CSF_end(); coeff_i++)
                {
                    for(auto coeff_j = proj_jt.CSF_begin(); coeff_j != proj_jt.CSF_end(); coeff_j++)
                    {
                        solution = 0;   // Note: solution = left_index * right_eigenvector.size() + right_index
                        for(unsigned int left_index = 0; left_index < left_eigenvector.size(); left_index++)
                        {
                            for(unsigned int right_index = 0; right_index < right_eigenvector.size(); right_index++)
                            {
                                coeff[solution] += (*coeff_i) * (*coeff_j)
                                    * left_eigenvector[left_index][coeff_i.index()]
                                    * right_eigenvector[right_index][coeff_j.index()];

                                solution++;
                            }
                        }
                    }
                }

                for(solution = 0; solution < return_size; solution++)
                    total[solution] += coeff[solution] * matrix_element;
            }
            proj_jt++;
        }
        proj_it++;
    }
    
    return total;
}

#endif

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
    {   static_assert(sizeof...(pElectronOperators) > 0, "ManyBodyOperator<> instantiated with too few operators (template arguments).");
        static_assert(sizeof...(pElectronOperators) < 4, "ManyBodyOperator<> instantiated with too few operators (template arguments).");
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
                    matrix_element += OneBodyMatrixElements(*e, *e);
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
                    matrix_element += OneBodyMatrixElements(*i, *i);

                    auto j = i;
                    while(j != end)
                    {
                        matrix_element += TwoBodyMatrixElements(*i, *j, *i, *j);
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
                    matrix_element += TwoBodyMatrixElements(*left.front(), *i, *right.front(), *i);
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
    indirect_proj.clear();
    indirect_proj.reserve(proj.size());

    std::copy(boost::make_counting_iterator(proj.data()),
              boost::make_counting_iterator(proj.data_end()),
              std::back_inserter(indirect_proj));
}

template<typename... pElectronOperators>
template<int max_diffs>
int ManyBodyOperator<pElectronOperators...>::GetProjectionDifferences(ManyBodyOperator<>::IndirectProjection& p1, ManyBodyOperator<>::IndirectProjection& p2) const
{
    // NB: These are static for speed: to prevent memory (de)allocation
    static IndirectProjection diff1, diff2;
    diff1.clear(); diff1.reserve(p1.size() + max_diffs);
    diff2.clear(); diff2.reserve(p2.size() + max_diffs);
    static IndirectProjection sorted_p1, sorted_p2;
    sorted_p1.clear(); sorted_p1.reserve(p1.size());
    sorted_p2.clear(); sorted_p2.reserve(p2.size());

    auto it1 = p1.begin();
    auto it2 = p2.begin();

    int permutations = 0;

    while(it1 != p1.end() && it2 != p2.end() && diff1.size() <= max_diffs && diff2.size() <= max_diffs)
    {
        const ElectronInfo& e1(**it1);
        const ElectronInfo& e2(**it2);

        if(e1 == e2)
        {   sorted_p1.push_back(*it1++);
            sorted_p2.push_back(*it2++);
        }
        else if(e1 < e2)
        {   permutations += *it1 - p1.front();
            diff1.push_back(*it1++);
        }
        else
        {   permutations += *it2 - p2.front();
            diff2.push_back(*it2++);
        }
    }
    while(it1 != p1.end() && (diff1.size() <= max_diffs))
    {   permutations += *it1 - p1.front();
        diff1.push_back(*it1++);
    }
    while(it2 != p2.end() && (diff2.size() <= max_diffs))
    {   permutations += *it2 - p2.front();
        diff2.push_back(*it2++);
    }

    int num_diffs = diff1.size();

    // TODO: Hole shifties
    if((diff1.size() > max_diffs) || (diff2.size() > max_diffs))
        return max_diffs+1;

    else if(diff1.empty() && diff2.empty())
        return 0;

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

#endif

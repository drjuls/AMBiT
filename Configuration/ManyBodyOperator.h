#ifndef MANY_BODY_OPERATOR_H
#define MANY_BODY_OPERATOR_H

#include "Projection.h"
#include "LevelMap.h"
#include <tuple>
#include <boost/iterator/counting_iterator.hpp>
#include <boost/iterator/indirect_iterator.hpp>
#include <boost/utility/enable_if.hpp>
#ifdef AMBIT_USE_MPI
#include <mpi.h>
#endif

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
#ifdef AMBIT_USE_OPENMP
        indirects_list.resize(omp_get_max_threads());
#endif
    }

    typedef std::vector<const ElectronInfo*> IndirectProjection;

    // Bundle all the indirect projections into a struct to allow for one indirect projection per thread
    struct IndirectProjections {
        IndirectProjection left, right;
        IndirectProjection diff1, diff2;
        IndirectProjection sorted_p1, sorted_p2;
    };

    /** Rearrange indirects.left and indirects.right so that differences are at the beginning (up to 
        max differences). abs(return value) is number of differences found sign(return value) 
        indicates number of permutations (+ve is even number, -ve is odd). 
        If skipped_right_electron != nullptr, assume indirects.right is one electron smaller than 
        indirects.left; the first element of indirects.right is the skipped electron which we can 
        assume is not present in indirects.left.
        
        N.B. This must be public to allow unit tests to call it
     */
    template<int max_diffs>
    inline int GetProjectionDifferences(IndirectProjections& indirects, const ElectronInfo* skipped_p2_electron = nullptr) const;
 
    inline void make_indirect_projection(const Projection& proj, IndirectProjection& indirect_proj) const;

    /** Matrix element between two projections.
        If eplsion != nullptr, we want
            < proj_left | O | {proj_right, epsilon} >
        where proj_left is one electron larger than proj_right.
     */
    inline double GetMatrixElement(const Projection& proj_left, const Projection& proj_right, const ElectronInfo* epsilon = nullptr) const;

    /** Returns <left | O | right>.
        If eplsion != nullptr, we want
            < proj_left | O | {proj_right, epsilon} >
        where proj_left is one electron larger than proj_right.
     */
    inline double GetMatrixElement(const Level& left, const Level& right, const ElectronInfo* epsilon = nullptr) const;

    /** Equivalent to GetMatrixElement(level, level). */
    inline double GetMatrixElement(const Level& level) const;

    /** Equivalent to calculating GetMatrixElement(level) for each level in vector.
        Return vector of matrix elements.
     */
    inline std::vector<double> GetMatrixElement(const LevelVector& levels) const;

    /** Equivalent to calculating GetMatrixElement(left, right, epsilon) for each pair of levels in their respective vectors.
        Return array of matrix elements indexed by left_index * (num_right) + right_index.
     */
    inline std::vector<double> GetMatrixElement(const LevelVector& left_levels, const LevelVector& right_levels, const ElectronInfo* epsilon = nullptr) const;

protected:
    std::tuple<pElectronOperators...> pOperators;

    
    // NB: Indirect projections are class members to prevent expensive memory (de)allocations
#ifdef AMBIT_USE_OPENMP
    mutable std::vector<IndirectProjections> indirects_list;
#else
    mutable IndirectProjections my_projections;
#endif
   
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

    bool IsMyJob(int index) const
    {   return index%NumProcessors == ProcessorRank;
    }
};

/** ManyBodyOperator: one body operator specialization. */
template <typename... pElectronOperators>
double ManyBodyOperator<pElectronOperators...>::GetMatrixElement(const Projection& proj_left, const Projection& proj_right, const ElectronInfo* epsilon) const
{
    int num_diffs = 0;


#ifdef AMBIT_USE_OPENMP
    // Grab the indirect projections for this thread
    IndirectProjections& my_projections = indirects_list[omp_get_thread_num()];
#endif
    make_indirect_projection(proj_left, my_projections.left); 

    // Skip this for same projection
    if(&proj_left != &proj_right)
    {
        make_indirect_projection(proj_right, my_projections.right);
        num_diffs = GetProjectionDifferences<sizeof...(pElectronOperators)>(my_projections, epsilon);
    }

    double matrix_element = 0.0;

    switch(sizeof...(pElectronOperators))
    {
        case 1:
            if(num_diffs == 0)
            {
                for(auto& e: my_projections.left)
                {
                    int sign = e->IsHole()? -1 : 1;
                    matrix_element += OneBodyMatrixElements(*e, *e) * sign;
                }
            }
            else if(abs(num_diffs) == 1)
            {
                matrix_element += OneBodyMatrixElements(*my_projections.left[0], *my_projections.right[0]);
            }
            break;
        case 2:
            if(num_diffs == 0)
            {
                auto i = boost::make_indirect_iterator(my_projections.left.begin());
                auto end = boost::make_indirect_iterator(my_projections.left.end());
                while(i != end)
                {
                    int sign = i->IsHole()? -1 : 1;
                    matrix_element += OneBodyMatrixElements(*i, *i) * sign;

                    auto j = std::next(i);
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
                matrix_element += OneBodyMatrixElements(*my_projections.left[0], *my_projections.right[0]);

                auto i = boost::make_indirect_iterator(my_projections.left.begin());
                auto end = boost::make_indirect_iterator(my_projections.left.end());
                i++;    // Skip first element
                while(i != end)
                {
                    int sign = i->IsHole()? -1 : 1;
                    matrix_element += TwoBodyMatrixElements(*my_projections.left.front(), *i, *my_projections.right.front(), *i) * sign;
                    i++;
                }
            }
            else if(abs(num_diffs) == 2)
            {
                matrix_element += TwoBodyMatrixElements(*my_projections.left[0], *my_projections.left[1], *my_projections.right[0], *my_projections.right[1]);
            }
            break;
        case 3:
            if(num_diffs == 0)
            {
                auto i = boost::make_indirect_iterator(my_projections.left.begin());
                auto end = boost::make_indirect_iterator(my_projections.left.end());
                while(i != end)
                {
                    int sign = i->IsHole()? -1 : 1;
                    matrix_element += OneBodyMatrixElements(*i, *i) * sign;

                    auto j = std::next(i);
                    while(j != end)
                    {
                        int total_sign = sign * (j->IsHole()? -1 : 1);
                        matrix_element += TwoBodyMatrixElements(*i, *j, *i, *j) * total_sign;

                        auto k = std::next(j);
                        while(k != end)
                        {
                            int threebody_sign = total_sign * (k->IsHole()? -1 : 1);
                            matrix_element += ThreeBodyMatrixElements(*i, *j, *k, *i, *j, *k) * threebody_sign;
                            k++;
                        }
                        j++;
                    }
                    i++;
                }
            }
            else if(abs(num_diffs) == 1)
            {
                matrix_element += OneBodyMatrixElements(*my_projections.left[0], *my_projections.right[0]);

                auto i = boost::make_indirect_iterator(my_projections.left.begin());
                auto end = boost::make_indirect_iterator(my_projections.left.end());
                i++;    // Skip first element
                while(i != end)
                {
                    int sign = i->IsHole()? -1 : 1;
                    matrix_element += TwoBodyMatrixElements(*my_projections.left.front(), *i, *my_projections.right.front(), *i) * sign;

                    auto j = std::next(i);
                    while(j != end)
                    {
                        int total_sign = sign * (j->IsHole()? -1 : 1);
                        matrix_element += ThreeBodyMatrixElements(*my_projections.left.front(), *i, *j, *my_projections.right.front(), *i, *j) * total_sign;
                        j++;
                    }
                    i++;
                }
            }
            else if(abs(num_diffs) == 2)
            {
                matrix_element += TwoBodyMatrixElements(*my_projections.left[0], *my_projections.left[1], *my_projections.right[0], *my_projections.right[1]);

                auto i = boost::make_indirect_iterator(my_projections.left.begin());
                auto end = boost::make_indirect_iterator(my_projections.left.end());
                std::advance(i, 2);    // Skip first two elements
                while(i != end)
                {
                    int sign = i->IsHole()? -1 : 1;
                    matrix_element += ThreeBodyMatrixElements(*my_projections.left[0], *my_projections.left[1], *i, *my_projections.right[0], *my_projections.right[1], *i) * sign;

                    i++;
                }
            }
            else if(abs(num_diffs) == 3)
            {
                matrix_element += ThreeBodyMatrixElements(*my_projections.left[0], *my_projections.left[1], *my_projections.left[2], *my_projections.right[0], *my_projections.right[1], *my_projections.right[2]);
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
int ManyBodyOperator<pElectronOperators...>::GetProjectionDifferences(IndirectProjections& indirects, const ElectronInfo* skipped_p2_electron) const
{
    indirects.diff1.clear();
    indirects.diff2.clear();
    indirects.sorted_p1.clear();
    indirects.sorted_p2.clear();

    auto it1 = indirects.left.begin();
    auto it2 = indirects.right.begin();

    // Number of permutations required to get correct ordering
    int permutations = 0;

    // Cumulative number of operators that are the same on both sides
    int num_same = 0;

    // Electrons and hole operators that have not yet been absorbed into a "difference"
    int num_holes1 = 0;
    int num_electrons1 = 0;
    int num_holes2 = 0;
    int num_electrons2 = 0;

    if(skipped_p2_electron != nullptr)
    {
        // Assume first element of left was a skipped electron (not hole)
        indirects.diff2.push_back(skipped_p2_electron);
        num_electrons2++;
    }

    while(indirects.diff1.size() <= max_diffs && indirects.diff2.size() <= max_diffs && it1 != indirects.left.end() && it2 != indirects.right.end())
    {
        const ElectronInfo& e1(**it1);
        const ElectronInfo& e2(**it2);

        if(e1 == e2)
        {   indirects.sorted_p1.push_back(*it1++);
            indirects.sorted_p2.push_back(*it2++);
            num_same++;
        }
        else if(e1 < e2)
        {
            permutations += num_same;
            if(e1.IsHole())
            {
                indirects.diff2.push_back(*it1++);

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
                indirects.diff1.push_back(*it1++);

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
                indirects.diff1.push_back(*it2++);

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
                indirects.diff2.push_back(*it2++);

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
    while(it1 != indirects.left.end() && (indirects.diff1.size() <= max_diffs))
    {
        permutations += num_same;
        if((*it1)->IsHole())
        {
            indirects.diff2.push_back(*it1++);

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
            indirects.diff1.push_back(*it1++);

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
    while(it2 != indirects.right.end() && (indirects.diff2.size() <= max_diffs))
    {
        permutations += num_same;
        if((*it2)->IsHole())
        {
            indirects.diff1.push_back(*it2++);

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
            indirects.diff2.push_back(*it2++);

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

    int num_diffs = indirects.diff1.size();

    if((indirects.diff1.size() > max_diffs) || (indirects.diff2.size() > max_diffs))
        return max_diffs+1;

    else if(indirects.diff1.empty() && indirects.diff2.empty())
        return 0;

    // Additional sign changes for holes; add to permutations
    it1 = indirects.diff1.begin();
    it2 = indirects.diff2.begin();
    while(it1 != indirects.diff1.end())
    {
        if((*it1)->IsHole() || (*it2)->IsHole())
            permutations++;

        it1++;
        it2++;
    }

    //Copy differences
    indirects.diff1.insert(indirects.diff1.end(), indirects.sorted_p1.begin(), indirects.sorted_p1.end());
    indirects.diff2.insert(indirects.diff2.end(), indirects.sorted_p2.begin(), indirects.sorted_p2.end());

    indirects.left = indirects.diff1;
    indirects.right = indirects.diff2;

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

    auto config_it = configs->begin();
    int config_index = 0;
    while(config_it != configs->end())
    {
        if(IsMyJob(config_index))
        {
            auto config_jt = config_it;
            while(config_jt != configs->end())
            {
                if(config_it->GetConfigDifferencesCount(*config_jt) <= sizeof...(pElectronOperators))
                {
                    // Iterate over projections
                    auto proj_it = config_it.projection_begin();
                    while(proj_it != config_it.projection_end())
                    {
                        RelativisticConfiguration::const_projection_iterator proj_jt;
                        if(config_it == config_jt)
                            proj_jt = proj_it;
                        else
                            proj_jt = config_jt.projection_begin();

                        while(proj_jt != config_jt.projection_end())
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
                }
                config_jt++;
            }
        }
        config_it++;
        config_index++;
    }

#ifdef AMBIT_USE_MPI
    double reduced_total = 0.;
    MPI_Allreduce(&total, &reduced_total, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return reduced_total;
#else
    return total;
#endif
}

template<typename... pElectronOperators>
std::vector<double> ManyBodyOperator<pElectronOperators...>::GetMatrixElement(const LevelVector& levels) const
{
    std::vector<const double*> eigenvector;

    eigenvector.reserve(levels.size());
    for(auto it = levels.begin(); it != levels.end(); it++)
        eigenvector.push_back((*it)->GetEigenvector().data());

    if(eigenvector.empty())
        return std::vector<double>();

    pRelativisticConfigListConst configs = levels.front()->GetRelativisticConfigList();
    std::vector<double> total(eigenvector.size(), 0.);

    unsigned int solution = 0;

    auto config_it = configs->begin();
    int config_index = 0;
    while(config_it != configs->end())
    {
        auto config_jt = config_it;
        while(config_jt != configs->end())
        {
            if(config_it->GetConfigDifferencesCount(*config_jt) <= sizeof...(pElectronOperators))
            {
                if(IsMyJob(config_index))
                {
                    // Iterate over projections
                    auto proj_it = config_it.projection_begin();
                    while(proj_it != config_it.projection_end())
                    {
                        RelativisticConfiguration::const_projection_iterator proj_jt;
                        if(config_it == config_jt)
                            proj_jt = proj_it;
                        else
                            proj_jt = config_jt.projection_begin();

                        while(proj_jt != config_jt.projection_end())
                        {
                            double matrix_element = GetMatrixElement(*proj_it, *proj_jt);

                            // coefficients
                            if(matrix_element)
                            {
                                // If the projections are different, count twice
                                if(proj_it != proj_jt)
                                {   matrix_element *= 2.;
                                }

                                for(solution = 0; solution < eigenvector.size(); solution++)
                                {
                                    for(auto coeff_i = proj_it.CSF_begin(); coeff_i != proj_it.CSF_end(); coeff_i++)
                                    {
                                        double left_coeff_and_matrix_element = matrix_element * (*coeff_i) * eigenvector[solution][coeff_i.index()];

                                        RelativisticConfigList::const_CSF_iterator start_j = proj_jt.CSF_begin();

                                        const double* pright = &eigenvector[solution][start_j.index()];
                                        for(auto coeff_j = start_j; coeff_j != proj_jt.CSF_end(); coeff_j++)
                                        {
                                            total[solution] += left_coeff_and_matrix_element
                                                                * (*coeff_j) * (*pright);
                                            pright++;
                                        }
                                    }
                                }
                            }
                            proj_jt++;
                        }
                        proj_it++;
                    }
                }
                config_index++;
            }
            config_jt++;
        }
        config_it++;
    }

#ifdef AMBIT_USE_MPI
    std::vector<double> reduced_total(total.size(), 0.);
    MPI_Allreduce(total.data(), reduced_total.data(), total.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return reduced_total;
#else
    return total;
#endif
}

template<typename... pElectronOperators>
double ManyBodyOperator<pElectronOperators...>::GetMatrixElement(const Level& left, const Level& right, const ElectronInfo* epsilon) const
{
    double total = 0.;
    const auto& configs_left = left.GetRelativisticConfigList();
    const auto& configs_right = right.GetRelativisticConfigList();
    const auto& eigenvector_left = left.GetEigenvector();
    const auto& eigenvector_right = right.GetEigenvector();

    auto config_it = configs_left->begin();
    int config_index = 0;
    while(config_it != configs_left->end())
    {
        if(IsMyJob(config_index))
        {
            auto config_jt = configs_right->begin();
            while(config_jt != configs_right->end())
            {
                if(config_it->GetConfigDifferencesCount(*config_jt) <= sizeof...(pElectronOperators))
                {
                    // Iterate over projections
                    auto proj_it = config_it.projection_begin();
                    while(proj_it != config_it.projection_end())
                    {
                        auto proj_jt = config_jt.projection_begin();
                        while(proj_jt != config_jt.projection_end())
                        {
                            double matrix_element = GetMatrixElement(*proj_it, *proj_jt, epsilon);

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
                }
                config_jt++;
            }
        }

        config_it++;
        config_index++;
    }

#ifdef AMBIT_USE_MPI
    double reduced_total = 0.;
    MPI_Allreduce(&total, &reduced_total, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return reduced_total;
#else
    return total;
#endif
}

template<typename... pElectronOperators>
std::vector<double> ManyBodyOperator<pElectronOperators...>::GetMatrixElement(const LevelVector& left_levels, const LevelVector& right_levels, const ElectronInfo* epsilon) const
{
    std::vector<const double*> left_eigenvector;
    std::vector<const double*> right_eigenvector;

    left_eigenvector.reserve(left_levels.size());
    for(auto it = left_levels.begin(); it != left_levels.end(); it++)
        left_eigenvector.push_back((*it)->GetEigenvector().data());

    right_eigenvector.reserve(right_levels.size());
    for(auto it = right_levels.begin(); it != right_levels.end(); it++)
        right_eigenvector.push_back((*it)->GetEigenvector().data());

    unsigned int return_size = left_eigenvector.size() * right_eigenvector.size();
    if(return_size == 0)
        return std::vector<double>();

    pRelativisticConfigListConst configs_left = left_levels.front()->GetRelativisticConfigList();
    pRelativisticConfigListConst configs_right = right_levels.front()->GetRelativisticConfigList();
    std::vector<double> total(return_size, 0.);

#ifdef AMBIT_USE_OPENMP
    /* Make a vector to hold the total for each thread. This needs to be static and thread private 
       so its contents persist across different OpenMP tasks (N.B. this is only here because gcc 
       and clang differ in how they handle data persistence between tasks)
    */
    static std::vector<double> my_total;
    #pragma omp threadprivate(my_total)
    
    // Clear the running totals *before* we enter the parallel region
    my_total.assign(return_size, 0.);
#endif

    unsigned int solution = 0;

    auto config_it = configs_left->begin();
    int config_index = 0;
#ifdef AMBIT_USE_OPENMP
    #pragma omp parallel copyin(my_total)
    {
    #pragma omp single 
#endif
    while(config_it != configs_left->end())
    {
        auto config_jt = configs_right->begin();
        while(config_jt != configs_right->end())
        {
            if(config_it->GetConfigDifferencesCount(*config_jt) <= sizeof...(pElectronOperators))
            {
                if(IsMyJob(config_index))
                {
#ifdef AMBIT_USE_OPENMP
                    #pragma omp task default(shared) firstprivate(config_it, config_jt, config_index, solution) shared(total)
                    {
#endif
                    // Iterate over projections
                    auto proj_it = config_it.projection_begin();
                    while(proj_it != config_it.projection_end())
                    {
                        int left_start_CSF_index = proj_it.CSF_begin().index();
                        int left_end_CSF_index = proj_it.CSF_end().index();

                        auto proj_jt = config_jt.projection_begin();
                        while(proj_jt != config_jt.projection_end())
                        {
                            double matrix_element = GetMatrixElement(*proj_it, *proj_jt, epsilon);

                            // coefficients
                            if(matrix_element)
                            {
                                int right_start_CSF_index = proj_jt.CSF_begin().index();
                                int right_end_CSF_index = proj_jt.CSF_end().index();

                                solution = 0;
                                for(unsigned int left_index = 0; left_index < left_eigenvector.size(); left_index++)
                                {
                                    for(unsigned int right_index = 0; right_index < right_eigenvector.size(); right_index++)
                                    {
                                        auto coeff_i = proj_it.CSF_begin();
                                        for(const double* pleft = &left_eigenvector[left_index][left_start_CSF_index];
                                            pleft != &left_eigenvector[left_index][left_end_CSF_index]; pleft++)
                                        {
                                            double left_coeff_and_matrix_element = matrix_element * (*coeff_i) * (*pleft);

                                            auto coeff_j = proj_jt.CSF_begin();
                                            for(const double* pright = &right_eigenvector[right_index][right_start_CSF_index];
                                                pright != &right_eigenvector[right_index][right_end_CSF_index]; pright++)
                                            {
#ifdef AMBIT_USE_OPENMP
                                                my_total[solution] += left_coeff_and_matrix_element * (*coeff_j) * (*pright);

#else
                                                total[solution] += left_coeff_and_matrix_element * (*coeff_j) * (*pright);
#endif
                                                coeff_j++;
                                            }
                                            coeff_i++;
                                        }

                                        solution++;
                                    }
                                }
                            }
                            proj_jt++;
                        }
                        proj_it++;
                    }
#ifdef AMBIT_USE_OPENMP
                    } // OpenMP task 
#endif
                } // MPI work distribution
                config_index++;
            }
            config_jt++;
        }
        config_it++;
    } // config_it + OpenMP single construct. There is an implicit OpenMP barrier here
#ifdef AMBIT_USE_OPENMP
    #pragma omp taskwait
    // Finally, gather all the thread totals into the final results
    #pragma omp critical(MANY_BODY_OP)
    for(int ii = 0; ii < return_size; ++ii)
    {
        total[ii] += my_total[ii];
    }
    } // OpenMP parallel region
#endif

#ifdef AMBIT_USE_MPI
    std::vector<double> reduced_total(return_size, 0.);
    MPI_Allreduce(total.data(), reduced_total.data(), return_size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return reduced_total;
#else
    return total;
#endif
}

#endif

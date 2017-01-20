#ifndef RELATIVISTIC_CONFIGURATION_H
#define RELATIVISTIC_CONFIGURATION_H

#include "HartreeFock/OrbitalInfo.h"
#include "HartreeFock/OrbitalMap.h"
#include "HartreeFock/Configuration.h"
#include "MBPT/OneElectronIntegrals.h"
#include "MBPT/SlaterIntegrals.h"
#include "Projection.h"
#include "AngularData.h"
#include "SortedList.h"
#include <boost/iterator_adaptors.hpp>
#include "IndexedIterator.h"
#include <string>

class RelativisticConfigList;

/** RelativisticConfiguration extends configuration by adding a set of projections
    and corresponding coefficients for a particular |J, M>.
 */
class RelativisticConfiguration: public Configuration<OrbitalInfo, int>
{
protected:
    friend class RelativisticConfigList;

public:
    RelativisticConfiguration() {}
    RelativisticConfiguration(const BaseConfiguration& other): BaseConfiguration(other) {}
    RelativisticConfiguration(BaseConfiguration&& other): BaseConfiguration(other) {}
    RelativisticConfiguration(const RelativisticConfiguration& other);
    RelativisticConfiguration(RelativisticConfiguration&& other);
    virtual ~RelativisticConfiguration() {}

    const RelativisticConfiguration& operator=(const RelativisticConfiguration& other);
    RelativisticConfiguration& operator=(RelativisticConfiguration&& other);

    typedef indexed_iterator<const double*> const_CSF_iterator;

    /** Return whether a suitable projection with (J, M) was found. */
    bool GetProjections(pAngularDataLibrary data, const Symmetry& sym, int two_m);

    /** Get the number of CSFs that have been calculated. */
    unsigned int NumCSFs() const;

    /** Calculate the largest projection possible for this configuration. */
    int GetTwiceMaxProjection() const;

    /** Calculate number of levels (of all J and M) corresponding to this configuration:
        \f[ \prod_a \frac{g_a!}{n_a!(g_a - n_a)!} \f]
     */
    int GetNumberOfLevels() const;

    /** Calculate configuration average energy. */
    double CalculateConfigurationAverageEnergy(pOrbitalMapConst orbitals, pHFIntegrals one_body, pSlaterIntegrals two_body) const;

    /** Return 2J if angular data has been set using GetProjections(), otherwise return -1. */
    int GetTwoJ() const
    {   if(angular_data)
            return angular_data->GetTwoJ();
        else return -1;
    }

    /** Return 2M if angular data has been set using GetProjections(). */
    int GetTwoM() const
    {   if(angular_data)
            return angular_data->GetTwoM();
        else return 0;
    }

public:
    /** Iterator over projections and CSFs. */
    class const_projection_iterator : public boost::iterator_adaptor<
        const_projection_iterator,
        ProjectionList::const_iterator,
        boost::use_default,
        boost::forward_traversal_tag>
    {
    public:
        const_projection_iterator():
            const_projection_iterator::iterator_adaptor_(), current_CSFs_start(nullptr), index_offset(0)
        {}

        const_projection_iterator(ProjectionList::const_iterator p):
            const_projection_iterator::iterator_adaptor_(p), current_CSFs_start(nullptr), index_offset(0)
        {}

        const_projection_iterator(ProjectionList::const_iterator p, const pAngularDataConst& angular_data, int offset = 0):
            const_projection_iterator::iterator_adaptor_(p), current_CSFs_start(angular_data->GetCSFs()), num_CSFs(angular_data->NumCSFs()), index_offset(offset)
        {}

        const_projection_iterator(const RelativisticConfiguration& rconfig, int offset = 0):
            const_projection_iterator::iterator_adaptor_(rconfig.projections.begin()), current_CSFs_start(rconfig.angular_data->GetCSFs()), num_CSFs(rconfig.angular_data->NumCSFs()), index_offset(offset)
        {}

        const_projection_iterator(const const_projection_iterator& other):
            const_projection_iterator::iterator_adaptor_(other.base_reference()), current_CSFs_start(other.current_CSFs_start), num_CSFs(other.num_CSFs), index_offset(other.index_offset)
        {}

        const_CSF_iterator CSF_begin() const { return const_CSF_iterator(current_CSFs_start, index_offset); }
        const_CSF_iterator CSF_end() const { return const_CSF_iterator(current_CSFs_start + num_CSFs, index_offset + num_CSFs); }
        unsigned int NumCSFs() const { return num_CSFs; }

    protected:
        friend class boost::iterator_core_access;
        void increment() { this->base_reference()++; current_CSFs_start += num_CSFs; }

    protected:
        const double* current_CSFs_start;
        unsigned int num_CSFs;
        int index_offset;
    };

    /** Usual begin() and end() are inherited from Configuration and iterate over pair(OrbitalInfo, occupancy).
        projection_begin(), projection_end(), etc, provide iterators over the projection list.
        csf_offset is an additional start index for the const_CSF_iterator to know where the CSF lies in
        the RelativisticConfigList. If using only one RelativisticConfiguration, use csf_offset = 0.
        NB: no default is provided for csf_offset to prevent accidental misuse from confusion with
            RelativisticConfigList::iterator.projection_begin().
     */
    const_projection_iterator projection_begin(int csf_offset) const
    {
        return const_projection_iterator(projections.begin(), angular_data, csf_offset);
    }

    const_projection_iterator projection_end(int csf_offset) const
    {
        return const_projection_iterator(projections.end(), angular_data, csf_offset);
    }

    unsigned int projection_size() const
    {
        if(!angular_data)
            return 0;

        return angular_data->projection_size();
    }

    void Read(FILE* fp);        //!< Read configuration only (angular data can be recovered from library)
    void Write(FILE* fp) const; //!< Write configuration only (not projections or angular data)

    /** Print Configuration::Name() to outstream with projections and, optionally, CSFs.
        Note that Configuration::operator<< prints only Name().
     */
    void Print(bool include_CSFs = false) const;

protected:
    /** Pointer to angular data (CSFs) for this RelativisticConfiguration. */
    pAngularData angular_data;

    /** Complete list of projections with required M, ordering taken directly from AngularData. */
    ProjectionList projections;
};

/** Calculate configuration average energy. */
template<class OccupancyType>
double CalculateConfigurationAverageEnergy(const Configuration<OrbitalInfo, OccupancyType>& config, pOrbitalMapConst orbitals, pHFIntegrals one_body, pSlaterIntegrals two_body)
{
    MathConstant* math = MathConstant::Instance();
    double energy = 0.;

    auto it_a = config.begin();
    while(it_a != config.end())
    {
        // one body: n_a E_a
        energy += it_a->second * one_body->GetMatrixElement(it_a->first, it_a->first);

        auto it_b = it_a;
        while(it_b != config.end())
        {
            double weight = it_a->second;
            if(it_a == it_b)
            {   if(it_a->second > 0)
                weight *= (it_b->second - 1)/2. * it_a->first.MaxNumElectrons()/double(it_a->first.MaxNumElectrons()-1.);
            else
                weight *= (it_b->second + 1)/2. * it_a->first.MaxNumElectrons()/double(it_a->first.MaxNumElectrons()-1.);
            }
            else
                weight *= it_b->second;

            if(weight)
            {
                // R^0_abab
                double U_ab = two_body->GetTwoElectronIntegral(0, it_a->first, it_b->first, it_a->first, it_b->first);

                // Sum_l R^k_abba (j_a  j_b k)^2 \xi(l_a + l_b + k)
                //                (1/2 -1/2 0)
                int k = (it_a->first.L() + it_b->first.L())%2;
                int kmax = (it_a->first.TwoJ() + it_b->first.TwoJ())/2;
                while(k <= kmax)
                {
                    double threej = math->Electron3j(it_a->first.TwoJ(), it_b->first.TwoJ(), k);
                    U_ab -= threej * threej * two_body->GetTwoElectronIntegral(k, it_a->first, it_b->first, it_b->first, it_a->first);

                    k += 2;
                }

                energy += weight * U_ab;
            }
            ++it_b;
        }
        ++it_a;
    }

    return energy;
}

#endif

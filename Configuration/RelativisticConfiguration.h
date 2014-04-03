#ifndef RELATIVISTIC_CONFIGURATION_H
#define RELATIVISTIC_CONFIGURATION_H

#include "HartreeFock/OrbitalInfo.h"
#include "HartreeFock/Configuration.h"
#include "Projection.h"
#include "AngularData.h"
#include "SortedList.h"
#include <boost/iterator_adaptors.hpp>
#include <string>

#ifndef PARITY_ENUM
#define PARITY_ENUM
    enum Parity { even, odd };
#endif

class RelativisticConfiguration: public Configuration<OrbitalInfo, int>
{
    /** RelativisticConfiguration extends configuration by adding a set of projections
        and corresponding coefficients for a particular |J, M>.
     */
public:
    RelativisticConfiguration() {}
    RelativisticConfiguration(const BaseConfiguration& other): BaseConfiguration(other) {}
    RelativisticConfiguration(BaseConfiguration&& other): BaseConfiguration(other) {}
    RelativisticConfiguration(const RelativisticConfiguration& other);
    RelativisticConfiguration(RelativisticConfiguration&& other);
    virtual ~RelativisticConfiguration() {}

    RelativisticConfiguration& operator=(const RelativisticConfiguration& other);
    RelativisticConfiguration& operator=(RelativisticConfiguration&& other);

    typedef const double* const_CSF_iterator;

    /** Return whether a suitable projection with J = M = two_m/2 was found. */
    bool GetProjections(pAngularDataLibrary data);

    /** Get the number of CSFs that have been calculated. */
    unsigned int NumCSFs() const;

    /** Calculate the largest projection possible for this configuration. */
    int GetTwiceMaxProjection() const;

public:
    class const_projection_iterator : public boost::iterator_facade<
        const_projection_iterator,
        const Projection,
        boost::bidirectional_traversal_tag
        >
    {
        /** Iterator over projections and CSFs. */
    public:
        const_projection_iterator():
        m_base(), current_CSFs(0)
        {}

        explicit const_projection_iterator(ProjectionList::const_iterator p, const_CSF_iterator csf, unsigned int num_CSFs):
            m_base(p), current_CSFs(csf), num_CSFs(num_CSFs)
        {}

        const_projection_iterator(const const_projection_iterator& other):
            m_base(other.m_base), current_CSFs(other.current_CSFs), num_CSFs(other.num_CSFs)
        {}

        const_CSF_iterator CSF_begin() const { return current_CSFs; }
        const_CSF_iterator CSF_end() const { return current_CSFs + num_CSFs; }
        unsigned int NumCSFs() const { return num_CSFs; }

    private:
        friend class boost::iterator_core_access;

        value_type dereference() const { return *m_base; }
        bool equal(const const_projection_iterator& other) const { return this->m_base == other.m_base; }
        void increment() { m_base++; current_CSFs += num_CSFs; }
        void decrement() { m_base--; current_CSFs -= num_CSFs; }

    private:
        ProjectionList::const_iterator m_base;
        AngularData::const_CSF_iterator current_CSFs;
        unsigned int num_CSFs;
    };

    const_projection_iterator projection_begin() const { return const_projection_iterator(projections.begin(), angular_data->CSF_begin(0), angular_data->NumCSFs()); }
    const_projection_iterator projection_end() const { return const_projection_iterator(projections.end(), angular_data->CSF_begin(angular_data->NumCSFs()), angular_data->NumCSFs()); }
    unsigned int projection_size() const { return angular_data->projection_size(); }

protected:
    /** Pointer to angular data (CSFs) for this RelativisticConfiguration. */
    pAngularData angular_data;

    /** Complete list of projections with required M, ordering taken directly from AngularData. */
    ProjectionList projections;
};

class MostCSFsFirstComparator;

class RelativisticConfigList: public SortedList<RelativisticConfiguration, MostCSFsFirstComparator>
{
public:
    RelativisticConfigList() {}
    RelativisticConfigList(const RelativisticConfigList& other): BaseSortedList(other), num_electrons(other.num_electrons) {}
    RelativisticConfigList(RelativisticConfigList&& other): BaseSortedList(other), num_electrons(other.num_electrons) {}
    RelativisticConfigList(const RelativisticConfiguration& val): BaseSortedList(val) { num_electrons = front().ElectronNumber(); }
    RelativisticConfigList(RelativisticConfiguration&& val): BaseSortedList(val) { num_electrons = front().ElectronNumber(); }
    virtual ~RelativisticConfigList() {}

    RelativisticConfigList& operator=(const RelativisticConfigList& other)
    {   BaseSortedList::operator=(other);
        num_electrons = other.num_electrons;
        return *this;
    }

    RelativisticConfigList& operator=(RelativisticConfigList&& other)
    {   BaseSortedList::operator=(other);
        num_electrons = other.num_electrons;
        return *this;
    }

protected:
    int num_electrons;
};

class MostCSFsFirstComparator
{
public:
    bool operator()(const RelativisticConfiguration& first, const RelativisticConfiguration& second) const
    {
        if(first.NumCSFs() > second.NumCSFs())
            return true;
        else if(first.NumCSFs() < second.NumCSFs())
            return false;
        else
            return(first < second);
    }
};

typedef boost::shared_ptr<RelativisticConfigList> pRelativisticConfigList;
typedef boost::shared_ptr<const RelativisticConfigList> pRelativisticConfigListConst;

#endif

#ifndef RELATIVISTIC_CONFIGURATION_H
#define RELATIVISTIC_CONFIGURATION_H

#include "HartreeFock/OrbitalInfo.h"
#include "Projection.h"
#include "AngularData.h"
#include <boost/iterator_adaptors.hpp>
#include <string>

#ifndef PARITY_ENUM
#define PARITY_ENUM
    enum Parity { even, odd };
#endif

class RelativisticConfiguration
{
//    friend struct RelConfProjectionSizeRanking;
//    friend struct RelConfNumJStatesRanking;

    /** RelativisticConfiguration extends configuration by adding a set of projections
        and corresponding coefficients for a particular |J, M>.
     */
public:
    RelativisticConfiguration() {}
    RelativisticConfiguration(const RelativisticConfiguration& other);
    virtual ~RelativisticConfiguration() {}

    typedef std::map<OrbitalInfo, int>::iterator iterator;
    typedef std::map<OrbitalInfo, int>::const_iterator const_iterator;
    typedef const double* const_CSF_iterator;

    iterator begin() { return config.begin(); }
    const_iterator begin() const { return config.begin(); }
    iterator end() { return config.end(); }
    const_iterator end() const { return config.end(); }
    iterator find(const OrbitalInfo& info);
    const_iterator find(const OrbitalInfo& info) const;

    void clear();
    bool empty() const;
    int size() const { return config.size(); }

    /** Different to STL map insert, in that it always changes the occupancy associated
        with OrbitalInfo, regardless of whether it existed previously.
        Returns iterator to new value, or end() if occupancy = 0.
     */
    iterator insert(const std::pair<OrbitalInfo, int>& val);

    iterator erase(const_iterator position);
    int erase(const OrbitalInfo& info);
    iterator erase(const_iterator first, const_iterator last);

    /** Get occupancy of a particular single particle state (zero if absent). */
    int GetOccupancy(const OrbitalInfo& info) const;
    int& operator[](const OrbitalInfo& info);

    /** Number of particles = number of electrons - number of holes. */
    virtual int ParticleNumber() const;

    /** Excitation number = number of electrons + number of holes. */
    virtual int ExcitationNumber() const;
    Parity GetParity() const;

    bool operator<(const RelativisticConfiguration& other) const;
    bool operator==(const RelativisticConfiguration& other) const;

    /** Return whether a suitable projection with J = M = two_m/2 was found. */
    bool GetProjections(pAngularDataLibrary data);

    /** Get the number of CSFs that have been calculated. */
    unsigned int NumCSFs() const;

    /** Calculate the largest projection possible for this configuration. */
    int GetTwiceMaxProjection() const;

    virtual std::string Name() const;
    
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
    std::map<OrbitalInfo, int> config;

    /** Pointer to angular data (CSFs) for this RelativisticConfiguration. */
    pAngularData angular_data;

    /** Complete list of projections with required M, ordering taken directly from AngularData. */
    ProjectionList projections;
};

typedef std::list<RelativisticConfiguration> RelativisticConfigList;
typedef boost::shared_ptr<RelativisticConfigList> pRelativisticConfigList;
typedef boost::shared_ptr<const RelativisticConfigList> pRelativisticConfigListConst;

//struct RelConfProjectionSizeRanking
//{
//    // Sort in descending order of projection
//    inline bool operator()(const RelativisticConfiguration& first, const RelativisticConfiguration& second) const
//    {
//        if(first.projections.size() > second.projections.size())
//            return true;
//        else if(first.projections.size() < second.projections.size())
//            return false;
//        else
//            return (first < second);
//    }
//};
//
//struct RelConfNumJStatesRanking
//{
//    // Sort in descending order of number of Jstates
//    inline bool operator()(const RelativisticConfiguration& first, const RelativisticConfiguration& second) const
//    {
//        if(first.num_states > second.num_states)
//            return true;
//        else if(first.num_states < second.num_states)
//            return false;
//        else
//            return (first < second);
//    }
//};

#endif

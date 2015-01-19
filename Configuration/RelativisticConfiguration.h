#ifndef RELATIVISTIC_CONFIGURATION_H
#define RELATIVISTIC_CONFIGURATION_H

#include "HartreeFock/OrbitalInfo.h"
#include "HartreeFock/Configuration.h"
#include "Projection.h"
#include "AngularData.h"
#include "SortedList.h"
#include <boost/iterator_adaptors.hpp>
#include "IndexedIterator.h"
#include <string>

/** RelativisticConfiguration extends configuration by adding a set of projections
    and corresponding coefficients for a particular |J, M>.
 */
class RelativisticConfiguration: public Configuration<OrbitalInfo, int>
{
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

    /** Return whether a suitable projection with J = M = two_m/2 was found. */
    bool GetProjections(pAngularDataLibrary data);

    /** Get the number of CSFs that have been calculated. */
    unsigned int NumCSFs() const;

    /** Calculate the largest projection possible for this configuration. */
    int GetTwiceMaxProjection() const;

public:
    /** Iterator over projections and CSFs. */
    class const_projection_iterator : public boost::iterator_facade<
        const_projection_iterator,
        const Projection,
        boost::forward_traversal_tag
        >
    {
    public:
        const_projection_iterator():
            m_base(), current_CSFs_start(0)
        {}

        explicit const_projection_iterator(ProjectionList::const_iterator p, const double* CSFs, unsigned int num_CSFs):
            m_base(p), current_CSFs_start(CSFs), num_CSFs(num_CSFs)
        {}

        const_projection_iterator(const RelativisticConfiguration& rconfig):
            m_base(rconfig.projections.begin()), current_CSFs_start(rconfig.angular_data->GetCSFs()), num_CSFs(rconfig.angular_data->NumCSFs())
        {}

        const_projection_iterator(const const_projection_iterator& other):
            m_base(other.m_base), current_CSFs_start(other.current_CSFs_start), num_CSFs(other.num_CSFs)
        {}

        const_CSF_iterator CSF_begin(int index_offset = 0) const { return const_CSF_iterator(current_CSFs_start, index_offset); }
        const_CSF_iterator CSF_end(int index_offset = 0) const { return const_CSF_iterator(current_CSFs_start + num_CSFs, index_offset + num_CSFs); }
        unsigned int NumCSFs() const { return num_CSFs; }

    protected:
        friend class boost::iterator_core_access;

        const value_type& dereference() const { return *m_base; }
        bool equal(const const_projection_iterator& other) const { return this->m_base == other.m_base; }
        void increment() { m_base++; current_CSFs_start += num_CSFs; }

    protected:
        ProjectionList::const_iterator m_base;
        const double* current_CSFs_start;
        unsigned int num_CSFs;
    };

    /** Usual begin() and end() are inherited from Configuration and iterate over pair(OrbitalInfo, occupancy).
        projection_begin(), projection_end(), etc, provide iterators over the projection list.
     */
    const_projection_iterator projection_begin() const
    {
        if(!angular_data)
            return const_projection_iterator(projections.begin(), nullptr, 0);

        return const_projection_iterator(projections.begin(), angular_data->GetCSFs(), angular_data->NumCSFs());
    }

    const_projection_iterator projection_end() const
    {
        if(!angular_data)
            return const_projection_iterator(projections.end(), nullptr, 0);

        return const_projection_iterator(projections.end(), angular_data->CSF_begin(angular_data->NumCSFs()), angular_data->NumCSFs());
    }

    unsigned int projection_size() const
    {
        if(!angular_data)
            return 0;

        return angular_data->projection_size();
    }

    void Read(FILE* fp);        //!< Read configuration only (angular data can be recovered from library)
    void Write(FILE* fp) const; //!< Write configuration only (not projections or angular data)

protected:
    /** Pointer to angular data (CSFs) for this RelativisticConfiguration. */
    pAngularData angular_data;

    /** Complete list of projections with required M, ordering taken directly from AngularData. */
    ProjectionList projections;
};

class MostCSFsFirstComparator;
class FewestProjectionsFirstComparator;

/** RelativisticConfigList extends SortedList to give projection, CSF information over the whole list.
    In particular a projection iterator just like RelativisticConfig::const_projection_iterator is
    provided which iterates over all projections in all configurations in the list, and from which CSF
    coefficients can be accessed.
    To facilitate multithreading over individual RelativisticConfigurations, a "random access"
    operator[] is provided that returns the RelativisticConfiguration iterator and the CSF offset.
 */
class RelativisticConfigList: public SortedList<RelativisticConfiguration, FewestProjectionsFirstComparator>
{
public:
    RelativisticConfigList() {}
    RelativisticConfigList(const RelativisticConfigList& other): BaseSortedList(other) {}
    RelativisticConfigList(RelativisticConfigList&& other): BaseSortedList(other) {}
    RelativisticConfigList(const RelativisticConfiguration& val): BaseSortedList(val) {}
    RelativisticConfigList(RelativisticConfiguration&& val): BaseSortedList(val) {}
    virtual ~RelativisticConfigList() {}

    const RelativisticConfigList& operator=(const RelativisticConfigList& other)
    {   BaseSortedList::operator=(other);
        return *this;
    }

    RelativisticConfigList& operator=(RelativisticConfigList&& other)
    {   BaseSortedList::operator=(other);
        return *this;
    }

    /** Return total number of CSFs stored in entire list. */
    unsigned int NumCSFs() const;

    /** Get the ith RelativisticConfiguration, and the CSF offset.
        PRE: i < size().
     */
    std::pair<iterator, int> operator[](unsigned int i);
    std::pair<const_iterator, int> operator[](unsigned int i) const;

public:
    typedef RelativisticConfiguration::const_CSF_iterator const_CSF_iterator;

    /** Iterator over projections and CSFs. */
    class const_projection_iterator : public boost::iterator_facade<
        const_projection_iterator,
        const Projection,
        boost::forward_traversal_tag
    >
    {
    public:
        const_projection_iterator():
            m_base(), m_csf_index(0), m_configlist_it(), m_configlist_end()
        {}

        explicit const_projection_iterator(const RelativisticConfigList* list, RelativisticConfiguration::const_projection_iterator proj_it, RelativisticConfigList::const_iterator list_it, int start_csf_index):
            m_base(proj_it), m_configlist_it(list_it), m_configlist_end(list->end()), m_csf_index(start_csf_index)
        {}

        const_projection_iterator(const RelativisticConfigList* list):
            m_base(), m_configlist_it(list->begin()), m_configlist_end(list->end()), m_csf_index(0)
        {
            // Check that there is a relativistic configuration to point to, otherwise m_base is null
            if(m_configlist_it != m_configlist_end)
            {
                m_base = m_configlist_it->projection_begin();
                // Advance to first projection (in case m_configlist_it->projection_size() == 0)
                advance_to_next_projection();
            }
        }

        const_projection_iterator(const const_projection_iterator& other):
            m_base(other.m_base), m_configlist_it(other.m_configlist_it), m_configlist_end(other.m_configlist_end), m_csf_index(other.m_csf_index)
        {}

        const_CSF_iterator CSF_begin() const { return m_base.CSF_begin(m_csf_index); }
        const_CSF_iterator CSF_end() const { return m_base.CSF_end(m_csf_index); }
        unsigned int NumCSFs() const { return m_base.NumCSFs(); }

    protected:
        friend class boost::iterator_core_access;

        const value_type& dereference() const { return *m_base; }

        bool equal(const const_projection_iterator& other) const { return this->m_base == other.m_base; }

        void increment()
        {
            m_base++;
            advance_to_next_projection();
        }

        void advance_to_next_projection()
        {
            // Keep going while we are at projection_end
            while(m_configlist_it != m_configlist_end &&
                  m_base == m_configlist_it->projection_end())
            {
                m_csf_index += m_configlist_it->NumCSFs();
                m_configlist_it++;

                if(m_configlist_it != m_configlist_end)
                    m_base = m_configlist_it->projection_begin();
            }
        }

    protected:
        RelativisticConfiguration::const_projection_iterator m_base;    //!< Base projection iterator
        RelativisticConfigList::const_iterator m_configlist_it;         //!< Current relativistic configuration
        RelativisticConfigList::const_iterator m_configlist_end;        //!< End of relativistic config list
        int m_csf_index;                                                //!< CSF index of current relativistic configuration
    };

    const_projection_iterator projection_begin() const;
    const_projection_iterator projection_end() const;

    /** Return total number of projections stored in entire list. */
    unsigned int projection_size() const;

    void Read(FILE* fp);        //!< Read configurations
    void Write(FILE* fp) const; //!< Write configurations
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

class FewestProjectionsFirstComparator
{
public:
    bool operator()(const RelativisticConfiguration& first, const RelativisticConfiguration& second) const
    {
        if(first.projection_size() < second.projection_size())
            return true;
        else if(first.projection_size() > second.projection_size())
            return false;
        else
            return(first < second);
    }
};

typedef std::shared_ptr<RelativisticConfigList> pRelativisticConfigList;
typedef std::shared_ptr<const RelativisticConfigList> pRelativisticConfigListConst;

#endif

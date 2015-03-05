#ifndef RELATIVISTIC_CONFIG_LIST_H
#define RELATIVISTIC_CONFIG_LIST_H

#include "RelativisticConfiguration.h"

class MostCSFsFirstComparator;
class FewestProjectionsFirstComparator;

/** RelativisticConfigList extends SortedList to give projection, CSF information over the whole list.
    There are two types of iterators provided:
    (const_)iterator: an iterator over the RelativisticConfigurations that has member functions
            projection_begin() and projection_end() that give a
            RelativisticConfiguration::const_projection_iterator over the current configuration.
    const_projection_iterator: which is exactly the same as RelativisticConfiguration::const_projection_iterator
            except it iterates over all projections in all configurations in the list and from which
            CSF coefficients can be accessed in exactly the same manner.

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

public: // Functions on whole list
    unsigned int NumCSFs() const;           //!< Total number of CSFs stored in entire list
    unsigned int projection_size() const;   //!< Total number of projections stored in entire list

    void Read(FILE* fp);            //!< Read configurations
    void Write(FILE* fp) const;     //!< Write configurations

public: // Functions for iterator and const_iterator
    /** Iterator over RelativisticConfigurations in list. Templating just for const/non-const iterators. */
    template <class Base>
    class relconfiglist_iterator : public boost::iterator_adaptor<
        relconfiglist_iterator<Base>,
        Base,
        boost::use_default,
        boost::bidirectional_traversal_tag>
    {
    private:
        struct enabler {};

    public:
        relconfiglist_iterator(): relconfiglist_iterator::iterator_adaptor_(), m_csf_index(0) {}

        explicit relconfiglist_iterator(Base p, int start_index = 0):
        relconfiglist_iterator::iterator_adaptor_(p), m_csf_index(start_index) {}

        template <class OtherBase>
        relconfiglist_iterator(relconfiglist_iterator<OtherBase> const& other,
                               typename boost::enable_if<boost::is_convertible<OtherBase,Base>, enabler>::type = enabler()):
        relconfiglist_iterator::iterator_adaptor_(other.base()), m_csf_index(other.m_csf_index) {}

        /** Get iterator over the projection list of the current RelativisticConfiguration. */
        RelativisticConfiguration::const_projection_iterator projection_begin() const
        {
            const Base& base = this->base_reference();
            return RelativisticConfiguration::const_projection_iterator(base->projections.begin(), base->angular_data, m_csf_index);
        }

        /** Get end() of the projection list of the current RelativisticConfiguration. */
        RelativisticConfiguration::const_projection_iterator projection_end() const
        {
            const Base& base = this->base_reference();
            return RelativisticConfiguration::const_projection_iterator(base->projections.end(), base->angular_data, m_csf_index);
        }

        /** Number of projections in the current RelativisticConfiguration. */
        unsigned int projection_size() const
        {
            pAngularDataConst angular_data = this->base_reference()->angular_data;
            if(!angular_data)
                return 0;

            return angular_data->projection_size();
        }

        /** Start index of current CSFs (i.e. those from current RelativisticConfiguration) in the complete RelativisticConfigList. */
        int csf_offset() const { return m_csf_index; }

    private:
        friend class boost::iterator_core_access;
        int m_csf_index;

        void increment()
        {   m_csf_index += this->base_reference()->NumCSFs();
            this->base_reference()++;
        }

        void decrement()
        {   this->base_reference()--;
            m_csf_index += this->base_reference()->NumCSFs();
        }
    };

    typedef relconfiglist_iterator<BaseSortedList::const_iterator> const_iterator;
    typedef relconfiglist_iterator<BaseSortedList::iterator> iterator;

    iterator begin() { return iterator(BaseSortedList::begin(), 0); }
    const_iterator begin() const { return const_iterator(BaseSortedList::begin(), 0); }

    iterator end() { return iterator(BaseSortedList::end()); }
    const_iterator end() const { return const_iterator(BaseSortedList::end()); }

    iterator erase(iterator position) { return iterator(m_list.erase(position.base())); }

    /** Get iterator for the ith RelativisticConfiguration (with correct CSF offset).
        PRE: i < size().
     */
    iterator operator[](unsigned int i);
    const_iterator operator[](unsigned int i) const;

public: // Functions for const_projection_iterator
    typedef RelativisticConfiguration::const_CSF_iterator const_CSF_iterator;

    class const_projection_iterator : public boost::iterator_adaptor<
        const_projection_iterator,
        RelativisticConfiguration::const_projection_iterator,
        boost::use_default,
        boost::forward_traversal_tag>
    {
    public:
        const_projection_iterator():
        const_projection_iterator::iterator_adaptor_(), m_csf_index(0), m_configlist_it(), m_configlist_end()
        {}

        explicit const_projection_iterator(const RelativisticConfigList* list, RelativisticConfiguration::const_projection_iterator proj_it, RelativisticConfigList::const_iterator list_it, int start_csf_index):
        const_projection_iterator::iterator_adaptor_(proj_it), m_configlist_it(list_it), m_configlist_end(list->end()), m_csf_index(start_csf_index)
        {}

        /** Iterator points to start of list. */
        const_projection_iterator(const RelativisticConfigList* list):
            const_projection_iterator::iterator_adaptor_(), m_configlist_it(list->begin()), m_configlist_end(list->end()), m_csf_index(0)
        {
            // Check that there is a relativistic configuration to point to, otherwise m_base is null
            if(m_configlist_it != m_configlist_end)
            {
                this->base_reference() = m_configlist_it.projection_begin();
                // Advance to first projection (in case m_configlist_it->projection_size() == 0)
                advance_to_next_projection();
            }
        }

        const_projection_iterator(const const_projection_iterator& other):
            const_projection_iterator::iterator_adaptor_(other.base_reference()), m_configlist_it(other.m_configlist_it), m_configlist_end(other.m_configlist_end), m_csf_index(other.m_csf_index)
        {}

        const_CSF_iterator CSF_begin() const { return this->base_reference().CSF_begin(); }
        const_CSF_iterator CSF_end() const { return this->base_reference().CSF_end(); }
        unsigned int NumCSFs() const { return this->base_reference().NumCSFs(); }

    protected:
        void increment() { this->base_reference()++; advance_to_next_projection(); }
        void advance_to_next_projection()
        {
            // Keep going while we are at projection_end
            while(m_configlist_it != m_configlist_end &&
                  this->base_reference() == m_configlist_it.projection_end())
            {
                m_csf_index += m_configlist_it->NumCSFs();
                m_configlist_it++;

                if(m_configlist_it != m_configlist_end)
                    this->base_reference() = m_configlist_it.projection_begin();
            }
        }

    private:
        friend class boost::iterator_core_access;
        RelativisticConfigList::const_iterator m_configlist_it;         //!< Current relativistic configuration
        RelativisticConfigList::const_iterator m_configlist_end;        //!< End of relativistic config list
        int m_csf_index;                                                //!< CSF index of current relativistic configuration
    };

    const_projection_iterator projection_begin() const;
    const_projection_iterator projection_end() const;
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

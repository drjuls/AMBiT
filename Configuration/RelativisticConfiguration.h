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

    /** Iterator over projections and CSFs. */
//    class projection_iterator : public boost::iterator_adaptor<
//        projection_iterator,
//    typename ProjectionList::iterator,   //Base
//    typename ConfigType::value_type  // Value
//    >
//    {
//    private:
//        struct enabler {};  // a private type avoids misuse
//
//    public:
//        config_iterator():
//        config_iterator::iterator_adaptor_(typename ConfigType::iterator()) {}
//
//        explicit config_iterator(typename ConfigType::iterator p):
//        config_iterator::iterator_adaptor_(p) {}
//
//        template <class OtherConfigType>
//        config_iterator(config_iterator<OtherConfigType> const& other,
//                        typename boost::enable_if< boost::is_convertible<OtherConfigType,ConfigType>, enabler
//                        >::type = enabler()
//                        ):
//        config_iterator::iterator_adaptor_(other.base()) {}
//    };
//
//    typedef config_iterator<std::map<OrbitalInfo, int> > iterator;
//    typedef config_iterator<const std::map<OrbitalInfo, int> > const_iterator;

//    bool GenerateJCoefficients(double J);

    /** PRE: num_Jstates > 0
             coefficients = double[num_Jstates * projections.size()]
        This class assumes responsibility for freeing the memory in coefficients.
     */
//    void SetJCoefficients(unsigned int num_Jstates, double* coefficients);
//
//    inline const ProjectionSet& GetProjections() const
//    {   return projections; }
//
//    inline unsigned int NumJStates() const
//    {   return num_states;  }
//
//    inline const double* GetJCoefficients() const
//    {   return j_coefficients; }
//
//    inline double* GetJCoefficients()
//    {   return j_coefficients; }

    int GetTwiceMaxProjection() const;

    virtual std::string Name() const;

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

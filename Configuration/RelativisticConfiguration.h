#ifndef RELATIVISTIC_CONFIGURATION_H
#define RELATIVISTIC_CONFIGURATION_H

#include "Configuration.h"
#include "HartreeFock/OrbitalInfo.h"
#include "Projection.h"
#include "AngularData.h"

#include <string>

class RelativisticConfiguration
{
//    friend struct RelConfProjectionSizeRanking;
//    friend struct RelConfNumJStatesRanking;

    /** RelativisticConfiguration extends configuration by adding a set of projections
        and corresponding coefficients for a particular |J, M>.
        map-like member functions are provided, as though RelativisticConfiguration
        is a map<OrbitalInfo, int>, however internally it is stored as
        list< pair<OrbitalInfo, int> > so that we can sort on occupancy as well.
        This will make some operations (e.g. find, insert) inefficient.
        Because we sort on occupancy also, operator[] is not implemented as changing the
        associated value could break the ordering.
     */
public:
    RelativisticConfiguration():
        num_states(0), j_coefficients(NULL)
    {}
    RelativisticConfiguration(const RelativisticConfiguration& other);
    virtual ~RelativisticConfiguration()
    {   if(num_states)
            delete[] j_coefficients;
    }

    typedef std::list< std::pair<OrbitalInfo, int> >::iterator iterator;
    typedef std::list< std::pair<OrbitalInfo, int> >::const_iterator const_iterator;

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

    /** Number of particles = number of electrons - number of holes. */
    virtual int ParticleNumber() const;

    /** Excitation number = number of electrons + number of holes. */
    virtual int ExcitationNumber() const;
    Parity GetParity() const;

    bool operator<(const RelativisticConfiguration& other) const;
    bool operator==(const RelativisticConfiguration& other) const;

    /** Return whether a suitable projection with J = M = two_m/2 was found.
     */
//    bool GenerateProjections(int two_m);
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
//
//    int GetTwiceMaxProjection() const;

//    virtual std::string Name() const;
    
    Configuration GetNonRelConfiguration() const;

    // File storage (binary)
//    virtual void Write(FILE* fp) const;
//    virtual void Read(FILE* fp);

protected:
    static bool compare(std::pair<OrbitalInfo, int>& one, std::pair<OrbitalInfo, int>& two);
//    double GetJSquared(const Projection& first, const Projection& second) const;

    /** A complete set of projections with M = J */
    ProjectionSet projections;

    /** A superposition of projections that form an eigenstate of J^2 is called
        a "JState". The set of coefficients for each JState is j_coefficients.
        The number of Jstates is num_jstates.
        j_coefficients = double[num_states * projections.size()]
      */
    unsigned int num_states;
    double* j_coefficients;

    std::list< std::pair<OrbitalInfo, int> > config;
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

#ifndef RELATIVISTIC_CONFIGURATION_H
#define RELATIVISTIC_CONFIGURATION_H

#include "Configuration.h"
#include "HartreeFock/StateInfo.h"
#include "Projection.h"
#include <string>

class RelativisticConfiguration : public Configuration
{
    friend struct RelConfProjectionSizeRanking;

    /** RelativisticConfiguration extends configuration by adding a set of projections
        and corresponding coefficients for a particular |J, M>. Thus it should not be used
        with non-relativistic configurations.
     */
public:
    RelativisticConfiguration():
        Configuration(), num_states(0), j_coefficients(NULL)
    {}
    RelativisticConfiguration(const RelativisticConfiguration& other);
    RelativisticConfiguration(const Configuration& other):
        Configuration(other), num_states(0), j_coefficients(NULL)
    {}
    virtual ~RelativisticConfiguration(void)
    {   if(num_states)
            delete[] j_coefficients;
    }

    StateInfo GetInfo() const;
    virtual bool AddSingleParticle(const StateInfo& info);

    /** Return whether a suitable projection with J = M = two_m/2 was found.
     */
    bool GenerateProjections(int two_m);
    bool GenerateJCoefficients(double J);

    /** PRE: num_Jstates > 0
             coefficients = double[num_Jstates * projections.size()]
        This class assumes responsibility for freeing the memory in coefficients.
     */
    void SetJCoefficients(unsigned int num_Jstates, double* coefficients);

    inline const ProjectionSet& GetProjections() const
    {   return projections; }
    inline unsigned int NumJStates() const
    {   return num_states;  }
    inline const double* GetJCoefficients() const
    {   return j_coefficients; }

    inline double* GetJCoefficients()
    {   return j_coefficients; }

    int GetTwiceMaxProjection() const;

    virtual std::string Name() const;
    
    Configuration GetNonRelConfiguration() const;

protected:
    void DoElectron(std::vector<ElectronInfo>& electrons, unsigned int index);

    double GetJSquared(const Projection& first, const Projection& second) const;

    /** A complete set of projections with M = J */
    ProjectionSet projections;

    /** A superposition of projections that form an eigenstate of J^2 is called
        a "JState". The set of coefficients for each JState is j_coefficients.
        The number of Jstates is num_jstates.
        j_coefficients = double[num_states * projections.size()]
      */
    unsigned int num_states;
    double* j_coefficients;
};

typedef std::list<RelativisticConfiguration> RelativisticConfigList;

struct RelConfProjectionSizeRanking
{
    // Sort in descending order of projection
    inline bool operator()(const RelativisticConfiguration& first, const RelativisticConfiguration& second) const
    {
        if(first.projections.size() > second.projections.size())
            return true;
        else if(first.projections.size() < second.projections.size())
            return false;
        else
            return (first < second);
    }
};

#endif

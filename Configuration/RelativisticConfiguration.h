#ifndef RELATIVISTIC_CONFIGURATION_H
#define RELATIVISTIC_CONFIGURATION_H

#include "Configuration.h"
#include "RelativisticInfo.h"
#include "Projection.h"
#include <string>

class RelativisticConfiguration : public Configuration
{
public:
    RelativisticConfiguration():
        Configuration()
    {}
    RelativisticConfiguration(const RelativisticConfiguration& other):
        Configuration(other), projections(other.projections), j_coefficients(other.j_coefficients)
    {}
    RelativisticConfiguration(const Configuration& other):
        Configuration(other)
    {}
    virtual ~RelativisticConfiguration(void) {}

    RelativisticInfo GetInfo() const;
    virtual bool AddSingleParticle(const SingleParticleInfo& info);

    /** Return whether a suitable projection with J = M = two_m/2 was found
     */
    bool GenerateProjections(int two_m);
    inline const ProjectionSet& GetProjections() const
    {   return projections; }
    inline const std::vector< std::vector<double> >& GetJCoefficients() const
    {   return j_coefficients; }

    int GetTwiceMaxProjection() const;

    virtual std::string Name() const;
    
    Configuration GetNonRelConfiguration() const;

protected:
    void DoElectron(std::vector<ElectronInfo>& electrons, unsigned int index);
    bool GetProjectionCoefficients(double J);

    double GetJSquared(const Projection& first, const Projection& second) const;

    /** A complete set of projections with M = J */
    ProjectionSet projections;

    /** A superposition of projections that form an eigenstate of J^2 is called
        a "JState". The set of coefficients for each JState is j_coefficients.
      */
    std::vector< std::vector<double> > j_coefficients;
};

typedef std::list<RelativisticConfiguration> RelativisticConfigList;

#endif
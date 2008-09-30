#ifndef PROJECTION_H
#define PROJECTION_H

#include "Configuration.h"
#include "ElectronInfo.h"
#include <vector>
#include <list>

#ifndef PARITY_ENUM
#define PARITY_ENUM
    enum Parity { even, odd };
#endif

class Projection
{
    /** A projection is really a kind of configuration, however there
        is no occupancy: either a state has one electron or it isn't
        part of the configuration.
     */
public:
    Projection() {}
    Projection(const Projection& other):
        Config(other.Config), NonRelConfiguration(other.NonRelConfiguration)
    {}
    virtual ~Projection(void) {}

    void Add(const ElectronInfo& info);
    void Remove(const ElectronInfo& info);

    unsigned int Size() const { return Config.size(); }
    ElectronInfo& operator[](unsigned int i);
    const ElectronInfo& operator[](unsigned int i) const;

    Parity GetParity() const;
    int GetTwoM() const;

    /** Sort the projection. Return true if there were an odd number of swaps. */
    bool Sort();

    bool operator<(const Projection& other) const;
    bool operator==(const Projection& other) const;
    std::string Name() const;

    /** This function finds differences between p1 and p2 (up to two differences).
        The absolute value of the return is the number of differences. The sign of the return
            indicates whether getting the differences to align required an odd or even number
            of swaps. If there are more than two differences the function returns 3.
        diff[4] is an array which will store the positions of the differences,
            diff[0] = projection p1, difference 1.
            diff[1] = projection p2, difference 1.
            diff[2] = projection p1, difference 2.
            diff[3] = projection p2, difference 2.
      */
    static int GetProjectionDifferences(const Projection& p1, const Projection& p2, unsigned int* diff);

    /** This function is the same as GetProjectionDifferences(), but returns up
            to three differences.
            If there are more than three differences the function returns 4.
        diff[6] stores the positions of the three differences (see above).
     */
    static int GetProjectionDifferences3(const Projection& p1, const Projection& p2, unsigned int* diff);

    /** Get equivalent non-relativistic configuration. */
    inline const Configuration& GetNonRelConfiguration() const;

protected:
    bool Sort(std::vector<ElectronInfo>& config) const;

    std::vector<ElectronInfo> Config;
    Configuration NonRelConfiguration;
};

const Configuration& Projection::GetNonRelConfiguration() const
{
    if(NonRelConfiguration.Empty())
        exit(44);
    return NonRelConfiguration;
}

typedef std::list<Projection> ProjectionSet;

#endif

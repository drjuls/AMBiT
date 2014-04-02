#ifndef PROJECTION_H
#define PROJECTION_H

#include "ElectronInfo.h"
#include <vector>
#include <list>

#ifndef PARITY_ENUM
#define PARITY_ENUM
    enum Parity { even, odd };
#endif

class RelativisticConfiguration;

class Projection
{
    /** A projection is really a kind of configuration, however there
        is no occupancy: either a state has one electron or it isn't
        part of the configuration.
     */
public:
    Projection(const RelativisticConfiguration& relconfig, const std::vector<int>& twoMs);
    Projection(const Projection& other);
    Projection(Projection&& other);
    virtual ~Projection(void) {}

    Projection& operator=(const Projection& other) { config = other.config; return *this; }
    Projection& operator=(Projection&& other) { config.swap(other.config); return *this; }

    typedef std::vector<ElectronInfo>::iterator iterator;
    typedef std::vector<ElectronInfo>::const_iterator const_iterator;

    iterator begin() { return config.begin(); }
    const_iterator begin() const { return config.begin(); }
    iterator end() { return config.end(); }
    const_iterator end() const { return config.end(); }

    unsigned int size() const { return config.size(); }
    ElectronInfo& operator[](unsigned int i);
    const ElectronInfo& operator[](unsigned int i) const;

    Parity GetParity() const;
    int GetTwoM() const;

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

protected:
    std::vector<ElectronInfo> config;
};

typedef std::list<Projection> ProjectionList;

#endif

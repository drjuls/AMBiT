#ifndef PROJECTION_H
#define PROJECTION_H

#include "ElectronInfo.h"
#include "HartreeFock/Configuration.h"
#include <vector>
#include <list>

namespace Ambit
{
class RelativisticConfiguration;

/** A projection is kind of like a configuration, however there
    is no occupancy: either a state has one electron or it isn't
    part of the configuration.
    Like std::vector, the elements (ElectronInfo) are guaranteed to be stored contiguously,
    so a pointer retrieved via, e.g. data(), can be offset to access other elements.
 */
class Projection
{
public:
    Projection(const RelativisticConfiguration& relconfig, const std::vector<int>& twoMs);
    virtual ~Projection() = default;

    typedef std::vector<ElectronInfo>::iterator iterator;
    typedef std::vector<ElectronInfo>::const_iterator const_iterator;

    iterator begin() { return config.begin(); }
    const_iterator begin() const { return config.begin(); }
    iterator end() { return config.end(); }
    const_iterator end() const { return config.end(); }

    /** Get pointer to front(). */
    const ElectronInfo* data() const { return config.data(); }

    /** Get pointer just past the last element. */
    const ElectronInfo* data_end() const { return config.data() + config.size(); }

    unsigned int size() const { return config.size(); }
    ElectronInfo& operator[](unsigned int i);
    const ElectronInfo& operator[](unsigned int i) const;

    Parity GetParity() const;
    int GetTwoM() const;

    std::string Name() const;

protected:
    std::vector<ElectronInfo> config;
};

typedef std::vector<Projection> ProjectionList;

}
#endif

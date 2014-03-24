#ifndef ANGULAR_DATA_H
#define ANGULAR_DATA_H

#include "HartreeFock/OrbitalInfo.h"
#include <list>
#include <boost/shared_ptr.hpp>
#include <unordered_map>
#include <boost/functional/hash.hpp>

class RelativisticConfiguration;

class AngularData
{
    /** Store projections and configuration state functions (CSF) corresponding to
        the angular part of a RelativisticConfiguration.
        Access is only provided by const iterators.
     */
    friend class AngularDataCollection;

public:
    AngularData(): CSFs(nullptr) {}
    AngularData(const RelativisticConfiguration& config, int two_m); //!< Generate projections but not CSFs.
    AngularData(const RelativisticConfiguration& config, int two_m, int two_j); //!< Generate projections and CSFs.
    ~AngularData();

//    const_projection_iterator;
//    const_CSF_iterator;

protected:
    bool GenerateProjections(const RelativisticConfiguration& config, int two_m);
    bool GenerateCSFs(int two_j);

    static bool ProjectionCompare(const std::vector<int>& first, const std::vector<int>& second);

    std::list< std::vector<int> > projections;
    double* CSFs;
};

typedef boost::shared_ptr<AngularData> pAngularData;
typedef boost::shared_ptr<const AngularData> pAngularDataConst;

class AngularDataLibrary
{
    /** Collection of AngularData elements, indexed by RelativisticConfiguration.
        The collection is stored on disk in
            AMBiT/AngularData/particle_number.two_m.two_j.angular
     */
public:
    AngularDataLibrary(int particle_number, int two_m, int two_j);
    ~AngularDataLibrary() {}

    /** PRE: config is sorted correctly. */
    pAngularData operator[](const RelativisticConfiguration& config);

    void Read();
    void Write() const;

protected:
    typedef std::vector<int> KeyType;
    KeyType GenerateKey(const RelativisticConfiguration& config) const;

    std::string filename;
    std::unordered_map<KeyType, pAngularData, boost::hash<KeyType> > library;
};

typedef boost::shared_ptr<AngularDataLibrary> pAngularDataLibrary;

#endif

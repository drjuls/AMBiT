#ifndef ELECTRON_INFO_H
#define ELECTRON_INFO_H

#include "RelativisticInfo.h"

class ElectronInfo : public RelativisticInfo
{
    /** Stores a projection of angular momentum, M.
        Since M is half-integer, use TwoM which is always integer.
     */
public:
    ElectronInfo(unsigned int principal_qn, int kappa, int two_m):
        RelativisticInfo(principal_qn, kappa)
    {
        this->two_m = two_m;
    }
    ElectronInfo(const ElectronInfo& other):
        RelativisticInfo(other), two_m(other.two_m)
    {}
    virtual ~ElectronInfo(void) {}

    virtual bool operator<(const ElectronInfo& other) const;
    virtual bool operator==(const ElectronInfo& other) const;
    virtual bool operator!=(const ElectronInfo& other) const;

    inline void SetTwoM(int two_m) { this->two_m = two_m; }
    inline int TwoM() const { return two_m; }
    inline double M() const { return double(two_m)/2.; }

    virtual std::string Name() const;

protected:
    int two_m;
};

#endif

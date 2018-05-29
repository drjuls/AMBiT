#ifndef ELECTRON_INFO_H
#define ELECTRON_INFO_H

#include "HartreeFock/OrbitalInfo.h"

class ElectronInfo : public OrbitalInfo
{
    /** Stores a projection of angular momentum, M.
        Since M is half-integer, use TwoM which is always integer.
     */
public:
    ElectronInfo(): OrbitalInfo() {};
    ElectronInfo(unsigned int principal_qn, int kappa, int two_m, bool is_hole = false):
        OrbitalInfo(principal_qn, kappa), two_m(two_m), is_hole(is_hole)
    {}
    ElectronInfo(const ElectronInfo& other) = default;
    virtual ~ElectronInfo() = default;

    virtual bool operator<(const ElectronInfo& other) const;
    virtual bool operator==(const ElectronInfo& other) const;
    virtual bool operator!=(const ElectronInfo& other) const;

    inline void SetTwoM(int two_m) { this->two_m = two_m; }
    inline int TwoM() const { return two_m; }
    inline double M() const { return double(two_m)/2.; }

    inline bool IsHole() const { return is_hole; }
    inline bool IsParticle() const { return !is_hole; }
    virtual std::string Name() const;

protected:
    int two_m = {0};
    bool is_hole = {false};
};

#endif

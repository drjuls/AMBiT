#ifndef DEBUG_H
#define DEBUG_H

class Debug
{
public:
    Debug();
    ~Debug() {}

public:
    inline bool DebugFirstBuild() const { return bFirstBuild; }
    inline void DebugFirstBuild(bool debugon) { bFirstBuild = debugon; }
    inline bool DebugHFIterations() const { return bHFIterations; }
    inline void DebugHFIterations(bool debugon) { bHFIterations = debugon; }
    inline bool DebugHFExcited() const { return bHFExcited; }
    inline void DebugHFExcited(bool debugon) { bHFExcited = debugon; }
    inline bool DebugHFContinuum() const { return bHFContinuum; }
    inline void DebugHFContinuum(bool debugon) { bHFContinuum = debugon; }
    inline bool DebugVolumeShift() const { return bVShift; }
    inline void DebugVolumeShift(bool debugon) { bVShift = debugon; }
    
public:
    /** Generic Options */
    inline void HartreeEnergyUnits(bool turnon) { bHartreeUnits = turnon; }
    inline bool HartreeEnergyUnits() const { return bHartreeUnits; }
    inline void InvCmEnergyUnits(bool turnon) { bInvCmUnits = turnon; }
    inline bool InvCmEnergyUnits() const { return bInvCmUnits; }

private:
    bool bFirstBuild;
    bool bHFIterations;
    bool bHFExcited;
    bool bHFContinuum;
    bool bVShift;
    
    bool bHartreeUnits;
    bool bInvCmUnits;
};

inline Debug::Debug()
{
    bFirstBuild = false;
    bHFIterations = false;
    bHFExcited = false;
    bHFContinuum = false;
    bVShift = false;

    bHartreeUnits = false;
    bInvCmUnits = false;
}

#endif
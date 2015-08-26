#ifndef DEBUG_H
#define DEBUG_H

#include <chrono>

class Debug
{
    /** Set of options for debugging and logging.
        We also hijack this class for directing output from outside of class Atom.
     */
    friend class Atom;
public:
    Debug();
    ~Debug() {}

public:
    /** Hartree-Fock */
    inline bool LogFirstBuild() const { return bFirstBuild; }
    inline bool LogHFInnerLoop() const { return bHFInnerLoop; }
    inline bool LogHFIterations() const { return bHFIterations; }
    inline bool OutputHFExcited() const { return bHFExcited; }
    inline bool LogHFContinuum() const { return bHFContinuum; }

    /** MBPT */
    inline bool LogMBPT() const { return bMBPT; }

    /** ScaLAPACK */
    inline bool LogScalapack() const { return bScalapack; }

    /** Rate calculator */
    inline bool LogAugerRate() const { return bAugerRate; }

    /** Generic Options */
    inline bool HartreeEnergyUnits() const { return bHartreeUnits; }
    inline bool InvCmEnergyUnits() const { return bInvCmUnits; }

public:
    inline void LogFirstBuild(bool debugon) { bFirstBuild = debugon; }
    inline void LogHFInnerLoop(bool debugon) { bHFInnerLoop = debugon; }
    inline void LogHFIterations(bool debugon) { bHFIterations = debugon; }
    inline void OutputHFExcited(bool debugon) { bHFExcited = debugon; }
    inline void LogHFContinuum(bool debugon) { bHFContinuum = debugon; }

    inline void LogMBPT(bool debugon) { bMBPT = debugon; }
    inline void LogScalapack(bool debugon) { bScalapack = debugon; }
    inline void LogAugerRate(bool debugon) { bAugerRate = debugon; }

    inline void HartreeEnergyUnits(bool turnon) { bHartreeUnits = turnon; }
    inline void InvCmEnergyUnits(bool turnon) { bInvCmUnits = turnon; }

public:
    /** Set marker for later use by GetInterval(). */
    inline void MarkTime() { mark_time = std::chrono::steady_clock::now(); }

    /** Get duration since last MarkTime() call. */
    inline std::chrono::steady_clock::duration GetInterval() const
    {   return std::chrono::steady_clock::now() - mark_time;
    }

    /** Get number of seconds since last MarkTime() call. */
    inline double GetIntervalInSeconds() const
    {   std::chrono::steady_clock::duration d = GetInterval();
        return std::chrono::duration_cast<std::chrono::duration<double, std::ratio<1,1>>>(d).count();
    }

private:
    bool bFirstBuild;
    bool bHFInnerLoop;
    bool bHFIterations;
    bool bHFExcited;
    bool bHFContinuum;

    bool bMBPT;
    bool bScalapack;
    bool bAugerRate;

    bool bHartreeUnits;
    bool bInvCmUnits;

    std::chrono::steady_clock::time_point mark_time;
};

inline Debug::Debug()
{
    bFirstBuild = false;
    bHFIterations = false;
    bHFExcited = false;
    bHFContinuum = false;

    bMBPT = false;
    bScalapack = false;

    bHartreeUnits = false;
    bInvCmUnits = false;

    MarkTime();
}

#endif

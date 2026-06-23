#ifndef AMBIT_TIMER_H
#define AMBIT_TIMER_H

#include <absl/time/clock.h>
#include <absl/time/time.h>

namespace Ambit {
/* Utilities for some built-in, coarse-grained profiling. Uses Abseil's functionality
 * under the hood, because it has nicer defaults and automatically converts between units
 * when printing.
 */

using duration = absl::Duration;
using time = absl::Time;

/* Struct to hold the start, end, and elapsed times for a given timer, as well as to
 * update the timestamps.
 *
 * Usage: we don't directly set the start or end times, but call the start() and stop()
 * methods, respectively, which get the current wall time and update the corresponding
 * variable. This currently uses Abseil under the hood but, hopefully, this should keep
 * the choice of timer backend somewhat abstracted so that the main body of the code
 * doesn't need to worry about how this works - just call the corresponding methods and
 * it should all Just Work (TM).
 */
class Timer
{
public:
    void start() { start_time = absl::Now(); };
    void stop()  { end_time = absl::Now(); };

    // Get the difference between the start and end wall times.
    // Pre: both start() and stop() **must** have been called, or else this will return
    // garbage results
    //
    // TODO EVK: Do we want to maybe make this a std::optional or something to have
    // cleaner failure modes?
    duration elapsed()
    {
        // Currently this just does a simple subtraction, but we might want to do more
        // clever things later on
        return(end_time - start_time);
    };

public:
    time start_time;
    time end_time;
};

/* struct to hold a set of wall time durations, each one corresponding to a different
 * section of the calculation. The members are of "Timer" type, which is defined above.
 *
 * This is a singleton object for now (as it needs to be available and unique for 
 * the entire program's lifetime), but we might eventually want to do something
 * smarter so that we can aggregate and analyse times across MPI processes.
 */

class Profiler {
public:
    static Profiler* Instance();

public:
    Timer total; // Total calculation
    Timer hf; // Hartree-Fock
    Timer slater; // Two-electron Slater integrals
    Timer one_body_mbpt; // One-Body MBPT
    Timer two_body_mbpt; // Two-Body MBPT
    Timer angular_momentum; // Angular momentum calculations
    Timer ci; // CI
    Timer transitions; // Transition matrix elements
protected:
    Profiler() {};
};

}; // namespace ambit

#endif // AMBIT_TIMER_H

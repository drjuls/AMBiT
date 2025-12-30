#include "Profiler.h"
#include <absl/time/clock.h>
#include <absl/time/time.h>

namespace Ambit
{
    Profiler* Profiler::Instance()
    {
        static Profiler instance;
        return &instance;
    };
};

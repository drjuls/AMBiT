#ifndef AMBIT_CLASS_H
#define AMBIT_CLASS_H

#include "Include.h"
#include "Atom/Atom.h"

namespace Ambit
{
struct AmbitInterface
{
    AmbitInterface(int Z, std::string id): atom(Z, id) {}

    void EnergyCalculations();
    void TransitionCalculations();
    void Recombination();
    void InternalConversion();

    static void PrintHelp(const std::string& ApplicationName);

    Atom atom;
};
}

#ifdef UNIX
// Signal handler: emits a stack-trace on SIGTERM or SIGSEGV. Only does something useful on UNIX systems
static void signal_handler(int signo);
#endif

#endif

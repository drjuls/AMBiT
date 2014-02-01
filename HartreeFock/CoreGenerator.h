#ifndef CORE_GENERATOR_H
#define CORE_GENERATOR_H

#include "Orbital.h"
#include "HFOperator.h"
#include "Atom/MultirunOptions.h"
#include <list>

/** Create HartreeFock operator (with relevant decorators) and
    core orbitals based on user input.
 */
class CoreGenerator
{
public:
    CoreGenerator(pLattice lat);
    virtual ~CoreGenerator();

    virtual void GenerateCore(MultirunOptions& userInput, Core* open_shell_core, Core* closed_shell_core);

    virtual pHFOperatorConst GetHFOperator() const { return hf; }

protected:
    pLattice lattice;
    pHFOperator hf;
};

#endif

#include "CIIntegrals.h"
#include "gtest/gtest.h"
#include "HartreeFock/StateManager.h"
#include "HartreeFock/ConfigurationParser.h"
#include "HartreeFock/HartreeFocker.h"
#include "Basis/BasisGenerator.h"
#include "Include.h"

TEST(CIIntegralsTester, StandardCoulombSlow)
{
    DebugOptions.LogFirstBuild(true);
    DebugOptions.LogHFIterations(true);

    pLattice lattice(new Lattice(1000, 1.e-6, 50.));
    MultirunOptions userInput("template.input", "//", "\n", ",");
    userInput.SetRun(1);

    BasisGenerator generator(lattice, userInput);
    pCore core = generator.GenerateHFCore();
    pStateManager valence = generator.GenerateBasis();

    pHFOperatorConst hf = generator.GetHFOperator();
    pCoulombOperator coulomb(new CoulombOperator(lattice));
    pHartreeY hartreeY(new HartreeY(hf->GetOPIntegrator(), coulomb));

    CIIntegrals integrals(hf, hartreeY, valence, "", true);

    unsigned int size = integrals.GetStorageSize();
    integrals.Update();

    EXPECT_EQ(size, integrals.size());
}

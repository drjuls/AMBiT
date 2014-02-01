#include "CoreGenerator.h"
#include "gtest/gtest.h"
#include "Include.h"

TEST(CoreGeneratorTester, StartCore)
{
    pLattice lattice(new Lattice(1000, 1.e-6, 50.));

    DebugOptions.LogFirstBuild(true);
    DebugOptions.LogHFIterations(true);

    MultirunOptions userInput("template.input", "//", "\n", ",");
    userInput.SetRun(1);

    CoreGenerator core_generator(lattice);
    Core open_core(lattice), closed_core(lattice);
    core_generator.GenerateCore(userInput, &open_core, &closed_core);
    pHFOperatorConst hf = core_generator.GetHFOperator();

    // Check that < 4s | h | 4s > = h.GetEnergy()
    pOrbitalConst core_4s = open_core.GetState(OrbitalInfo(4, -1));
    EXPECT_NEAR(1.0, core_4s->Norm(lattice), 1.e-8);
    EXPECT_NEAR(-1.333798, open_core.GetState(OrbitalInfo(3, -2))->GetEnergy(), 1.e-6 * 1.33380);
    EXPECT_NEAR(core_4s->GetEnergy(), hf->GetMatrixElement(*core_4s, *core_4s), 1.e-6 * fabs(core_4s->GetEnergy()));
}
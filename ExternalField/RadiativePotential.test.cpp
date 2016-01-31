#include "gtest/gtest.h"
#include "Include.h"
#include "RadiativePotential.h"
#include "Atom/MultirunOptions.h"
#include "Basis/BasisGenerator.h"

TEST(RadiativePotentialTester, Uehling)
{
    pLattice lattice(new Lattice(1000, 1.e-6, 50.));

    // Na - comparison with Ginges & Berengut
    std::string user_input_string = std::string() +
        "NuclearRadius = 2.9374\n" +
        "NuclearThickness = 2.3\n" +
        "Z = 11\n" +
        "[HF]\n" +
        "N = 10\n" +
        "Configuration = '1s2 2s2 2p6'\n" +
        "[Basis]\n" +
        "--hf-basis\n" +
        "ValenceBasis = 3sp\n";

    std::stringstream user_input_stream(user_input_string);
    MultirunOptions userInput(user_input_stream, "//", "\n", ",");

    // Get core and excited basis
    BasisGenerator basis_generator(lattice, userInput);
    pCore core = basis_generator.GenerateHFCore();
    pOrbitalManagerConst orbitals = basis_generator.GenerateBasis();
    pPhysicalConstant constants = basis_generator.GetPhysicalConstant();

    pHFOperator hf = basis_generator.GetClosedHFOperator();
    pHFOperator ueh = pHFOperator(new UehlingDecorator(hf, basis_generator.GetNuclearRMSRadius()));

    pOrbitalConst s = orbitals->valence->GetState(OrbitalInfo(3, -1));
    double step_ueh = ueh->GetMatrixElement(*s, *s) - hf->GetMatrixElement(*s, *s);

    EXPECT_NEAR(-5.559e-7, step_ueh, 1.e-10);

    // Check density function
    pNucleusDecorator nuc = basis_generator.GetNucleusDecorator();
    ueh = pHFOperator(new UehlingDecorator(hf, nuc->GetNuclearDensity()));
    EXPECT_NEAR(step_ueh, ueh->GetMatrixElement(*s, *s) - hf->GetMatrixElement(*s, *s), 1.e-10);

    // Check electron density
//    RadialFunction density = nuc->GetNuclearDensity();
//    for(const auto& orb: *core)
//    {
//        density -= orb.second->GetDensity() * core->GetOccupancy(orb.first);
//    }
//
//    ueh.reset(new UehlingDecorator(hf, density));
//    *logstream << "Electron Uehling: " << ueh->GetMatrixElement(*s, *s) << std::endl;
}

TEST(RadiativePotentialTester, UehlingCs)
{
    pLattice lattice(new Lattice(2000, 1.e-6, 50.));

    // Cs
    std::string user_input_string = std::string() +
        "NuclearRadius = 5.67073\n" +
        "NuclearThickness = 2.3\n" +
        "Z = 55\n" +
        "[HF]\n" +
        "--local-exchange\n" +
        "Xalpha = 0.0\n" +
        "N = 54\n" +
        "Configuration = '1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d10 5s2 5p6'\n" +
        "[Basis]\n" +
        "--hf-basis\n" +
        "ValenceBasis = 6s\n";

    std::stringstream user_input_stream(user_input_string);
    MultirunOptions userInput(user_input_stream, "//", "\n", ",");

    // Get core and excited basis
    BasisGenerator basis_generator(lattice, userInput);
    pCore core = basis_generator.GenerateHFCore();
    pOrbitalManagerConst orbitals = basis_generator.GenerateBasis();
    pPhysicalConstant constants = basis_generator.GetPhysicalConstant();

    pHFOperator hf = basis_generator.GetClosedHFOperator();
    pHFOperator ueh = std::make_shared<UehlingDecorator>(hf, basis_generator.GetNuclearRMSRadius());

    pOrbitalConst s = orbitals->valence->GetState(OrbitalInfo(6, -1));
    double E0 = hf->GetMatrixElement(*s, *s);
    double step_ueh = ueh->GetMatrixElement(*s, *s) - E0;

    EXPECT_NEAR(-1.1312e-5, step_ueh, 1.e-9);

    // Check density function
    pNucleusDecorator nuc = basis_generator.GetNucleusDecorator();
    ueh = std::make_shared<UehlingDecorator>(hf, nuc->GetNuclearDensity());
    EXPECT_NEAR(step_ueh, ueh->GetMatrixElement(*s, *s) - E0, 1.e-10);

    // Check electron density
//    RadialFunction density;
//    for(const auto& orb: *core)
//    {
//        density -= orb.second->GetDensity() * core->GetOccupancy(orb.first);
//    }
//
//    ueh = std::make_shared<UehlingDecorator>(hf, density);
//    double el_ueh = ueh->GetMatrixElement(*s, *s) - E0;
//    EXPECT_NEAR(0.78175e-7, el_ueh, 1.e-10);
//
//    density += nuc->GetNuclearDensity();
//    ueh = std::make_shared<UehlingDecorator>(hf, density);
//    EXPECT_NEAR(el_ueh, ueh->GetMatrixElement(*s, *s) - (E0 + step_ueh), 1.e-10);
}

TEST(RadiativePotentialTester, Magnetic)
{
    pLattice lattice(new Lattice(1000, 1.e-6, 50.));

    // Cs - comparison with Ginges & Berengut
    std::string user_input_string = std::string() +
        "NuclearRadius = 5.67073\n" +
        "NuclearThickness = 2.3\n" +
        "Z = 55\n" +
        "[HF]\n" +
        "--local-exchange\n" +
        "Xalpha = 0.666666666667\n" +
        "N = 55\n" +
        "Configuration = '1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d10 5s2 5p6: 6s1'\n" +
        "[Basis]\n" +
        "--hf-basis\n" +
        "ValenceBasis = 6s\n";

    std::stringstream user_input_stream(user_input_string);
    MultirunOptions userInput(user_input_stream, "//", "\n", ",");

    // Get core and excited basis
    BasisGenerator basis_generator(lattice, userInput);
    pCore core = basis_generator.GenerateHFCore();
    pOrbitalManagerConst orbitals = basis_generator.GenerateBasis();
    pPhysicalConstant constants = basis_generator.GetPhysicalConstant();

    pHFOperator hf = basis_generator.GetOpenHFOperator();
    pHFOperator mag = pHFOperator(new MagneticSelfEnergyDecorator(hf, basis_generator.GetNuclearRMSRadius()));

    pOrbitalConst s = orbitals->valence->GetState(OrbitalInfo(6, -1));
    double step_mag = mag->GetMatrixElement(*s, *s) - hf->GetMatrixElement(*s, *s);

    // Check step function
    EXPECT_NEAR(1.390e-5, step_mag, 1.e-8);

    // Check density function
    pNucleusDecorator nuc = basis_generator.GetNucleusDecorator();
    mag = pHFOperator(new MagneticSelfEnergyDecorator(hf, nuc->GetNuclearDensity()));
    EXPECT_NEAR(step_mag, mag->GetMatrixElement(*s, *s) - hf->GetMatrixElement(*s, *s), 1.e-10);

    // Add electron density
    RadialFunction density = nuc->GetNuclearDensity();
    for(const auto& orb: *core)
    {
        density -= orb.second->GetDensity() * core->GetOccupancy(orb.first);
    }

    mag = pHFOperator(new MagneticSelfEnergyDecorator(hf, density));
    EXPECT_NEAR(1.377e-5, mag->GetMatrixElement(*s, *s) - hf->GetMatrixElement(*s, *s), 1.e-8);
}

TEST(RadiativePotentialTester, ElectricSlow)
{
    pLattice lattice(new Lattice(1000, 1.e-6, 50.));

    // Cs - comparison with Ginges & Berengut
    std::string user_input_string = std::string() +
        "NuclearRadius = 5.67073\n" +
        "NuclearThickness = 2.3\n" +
        "Z = 55\n" +
        "[HF]\n" +
        "--local-exchange\n" +
        "Xalpha = 0.666666666667\n" +
        "N = 55\n" +
        "Configuration = '1s2 2s2 2p6 3s2 3p6 3d10 4s2 4p6 4d10 5s2 5p6: 6s1'\n" +
        "[Basis]\n" +
        "--hf-basis\n" +
        "ValenceBasis = 6s\n";

    std::stringstream user_input_stream(user_input_string);
    MultirunOptions userInput(user_input_stream, "//", "\n", ",");

    // Get core and excited basis
    BasisGenerator basis_generator(lattice, userInput);
    pCore core = basis_generator.GenerateHFCore();
    pOrbitalManagerConst orbitals = basis_generator.GenerateBasis();
    pPhysicalConstant constants = basis_generator.GetPhysicalConstant();

    pHFOperator hf = basis_generator.GetOpenHFOperator();
    pOrbitalConst s = orbitals->valence->GetState(OrbitalInfo(6, -1));

    // Check step function: external offmass term
    pHFOperator electric = pHFOperator(new ElectricSelfEnergyDecorator(hf, basis_generator.GetNuclearRMSRadius(), false));
    double step_electric_ext = electric->GetMatrixElement(*s, *s) - hf->GetMatrixElement(*s, *s);
    EXPECT_NEAR(7.361e-5, step_electric_ext, 2.e-8);

    // Check step function: integrated offmass term
    electric = pHFOperator(new ElectricSelfEnergyDecorator(hf, basis_generator.GetNuclearRMSRadius(), true));
    double step_electric_int = electric->GetMatrixElement(*s, *s) - hf->GetMatrixElement(*s, *s);
    EXPECT_NEAR(7.346e-5, step_electric_int, 1.e-7);    // Higher tolerance due to sensitivity to Rn

    // Check density function
    pNucleusDecorator nuc = basis_generator.GetNucleusDecorator();
    electric = pHFOperator(new ElectricSelfEnergyDecorator(hf, nuc->GetNuclearDensity()));
    EXPECT_NEAR(7.346e-5, electric->GetMatrixElement(*s, *s) - hf->GetMatrixElement(*s, *s), 1.e-8);
}

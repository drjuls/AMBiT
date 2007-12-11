#ifdef _MPI
#include <mpi.h>
#endif
#include "Include.h"
#include "Atom.h"
#include "Universal/Constant.h"
#include "HartreeFock/NonRelInfo.h"
#include "Configuration/ConfigGenerator.h"
#include "Configuration/HamiltonianMatrix.h"
#include "Configuration/MPIHamiltonianMatrix.h"
#include "Basis/BSplineBasis.h"

void Atom::RunMultipleOpen()
{
    DebugOptions.LogFirstBuild(false);
    DebugOptions.LogHFIterations(false);
    DebugOptions.OutputHFExcited(true);
    DebugOptions.HartreeEnergyUnits(true);
    DebugOptions.LogMBPT(false);

    //CreateCustomBasis();
    //CreateRBasis();
    CreateBSplineBasis();

    DebugOptions.OutputHFExcited(false);

    multiple_SMS = true;
    multiple_alpha = false;
    multiple_volume = false;
    multiple_radius = false;

    unsigned int i;

    if(multiple_SMS)
    {   // Inverse nuclear mass
        double parameter = -0.002;

        for(i = 0; i < 5; i++)
        {
            std::stringstream ss;
            multiple_parameters.push_back(parameter);
            ss << (int)(parameter * 1000);
            multiple_ids.push_back(identifier + '_' + ss.str());

            parameter += 0.001;
        }
    }
    else if(multiple_alpha)
    {   // x = (alpha/alpha0)^2 - 1
        alpha0 = Constant::Alpha;
        double parameter = -0.25;

        for(i = 0; i < 5; i++)
        {
            std::stringstream ss;
            multiple_parameters.push_back(parameter);
            ss << (int)(parameter * 8);
            multiple_ids.push_back(identifier + '_' + ss.str());

            parameter += 0.125;
        }
    }
    else if(multiple_volume)
    {   // Field shift multiplier
        double parameter = -100.;

        for(i = 0; i < 5; i++)
        {
            std::stringstream ss;
            multiple_parameters.push_back(parameter);
            ss << (int)(parameter);
            multiple_ids.push_back(identifier + '_' + ss.str());

            parameter += 50.;
        }

        double deltaR = 0.05;
        core->CalculateVolumeShiftPotential(deltaR/Constant::AtomicToFermi);
        *outstream << "\nThickness = " << std::setprecision(3) << core->GetNuclearThickness()*Constant::AtomicToFermi;
        *outstream << "\nRadius    = " << std::setprecision(3) << core->GetNuclearRadius()*Constant::AtomicToFermi;
        *outstream << "\nd(Radius) = " << std::setprecision(3) << deltaR;
    }
    else if(multiple_radius)
    {   // Nuclear radius in fermi
        double parameter = 0.;
        
        for(i = 0; i < 5; i++)
        {
            std::stringstream ss;
            multiple_parameters.push_back(parameter);
            ss << (int)(parameter);
            multiple_ids.push_back(identifier + '_' + ss.str());

            parameter += 3.;
        }

        *outstream << "\nThickness = " << std::setprecision(3) << core->GetNuclearThickness()*Constant::AtomicToFermi;
    }

    // delta = E_CI - E_HF
    //double[] mbpt_delta = {0.0, 0.0, 0.0, 0.0, 0.0};

    //GenerateMultipleIntegralsMBPT(true, false);
    //CollateMultipleIntegralsMBPT(32);
    GenerateMultipleIntegrals(true);
    ChooseSymmetries();

    // Uncomment to include sigma3.
    //sigma3 = new Sigma3Calculator(lattice, core, excited);
    //sigma3->SetEnergyShift(mbpt_delta/Constant::HartreeEnergy_cm);

    //CheckMatrixSizes();
    CalculateMultipleEnergies();
}

void Atom::GenerateMultipleIntegralsMBPT(bool CoreMBPT, bool ValenceMBPT, double* delta)
{
    integrals = new CIIntegralsMBPT(*excited);
    integralsMBPT = dynamic_cast<CIIntegralsMBPT*>(integrals);

    // Create excited state basis. Should be a superset of the CI basis.
    excited_mbpt = new BSplineBasis(lattice, core);
    excited_mbpt->SetIdentifier(&identifier);
    dynamic_cast<BSplineBasis*>(excited_mbpt)->SetParameters(40, 7, 45.);
    std::vector<unsigned int> num_states_per_l;
    num_states_per_l.push_back(29);
    num_states_per_l.push_back(29);
    num_states_per_l.push_back(30);
    num_states_per_l.push_back(29);
    num_states_per_l.push_back(28);
    excited_mbpt->CreateExcitedStates(num_states_per_l);

    if(CoreMBPT)
    {   mbpt = new MBPTCalculator(lattice, core, excited_mbpt);
        integralsMBPT->IncludeMBPT1(true, mbpt);
        integralsMBPT->IncludeMBPT2(true, mbpt);
        integralsMBPT->IncludeExtraBoxDiagrams(true);
    }
    else
    {   if(mbpt)
            delete mbpt;
        mbpt = NULL;
    }

    if(ValenceMBPT)
    {   valence_mbpt = new ValenceCalculator(lattice, core, excited_mbpt);
        integralsMBPT->IncludeValenceMBPT1(true, valence_mbpt);
        integralsMBPT->IncludeValenceMBPT2(true, valence_mbpt);
        integralsMBPT->IncludeValenceExtraBoxDiagrams(true);
    }
    else
    {   if(valence_mbpt)
            delete valence_mbpt;
        valence_mbpt = NULL;
    }

    integralsMBPT->SetTwoElectronStorageLimits(4, 4);

    // Affects both core and valence MBPT if extra box diagrams are included.
    // To include box diagrams in Hamiltonian, uncomment the #defines at the top of HamiltonianMatrix.cpp.
    integralsMBPT->SetExtraBoxDiagramLimits(4, 4);

    //integrals->GetStorageSize();
    //return;

    for(unsigned int i = 0; i < multiple_parameters.size(); i++)
    {
        if(mbpt && delta)
            mbpt->SetEnergyShift(delta[i]/Constant::HartreeEnergy_cm);
        if(valence_mbpt && delta)
            valence_mbpt->SetEnergyShift(delta[i]/Constant::HartreeEnergy_cm);

        SetMultipleIntegralsAndCore(i);
    }
}

void Atom::CollateMultipleIntegralsMBPT(unsigned int num_processors)
{
    if(!integrals)
    {   integrals = new CIIntegralsMBPT(*excited);
        integralsMBPT = dynamic_cast<CIIntegralsMBPT*>(integrals);
    }
    else if(integralsMBPT)
    {   // Stop doing MBPT calculations
        integralsMBPT->IncludeMBPT1(false);
        integralsMBPT->IncludeMBPT2(false);
        integralsMBPT->IncludeExtraBoxDiagrams(false);
        integralsMBPT->IncludeValenceMBPT1(false);
        integralsMBPT->IncludeValenceMBPT2(false);
        integralsMBPT->IncludeValenceExtraBoxDiagrams(false);
    }

    for(unsigned int i = 0; i < multiple_ids.size(); i++)
    {
        integrals->SetIdentifier(multiple_ids[i]);
        integrals->Clear();

        if(ProcessorRank == 0)
        {   integralsMBPT->ReadMultipleOneElectronIntegrals(multiple_ids[i], num_processors);
            integralsMBPT->ReadMultipleTwoElectronIntegrals(multiple_ids[i], num_processors);
            integrals->WriteOneElectronIntegrals();
            integrals->WriteTwoElectronIntegrals();
        }
    }
}

void Atom::GenerateMultipleIntegrals(bool MBPT_CI)
{
    if(MBPT_CI)
        integrals = new CIIntegralsMBPT(*excited);
    else
        integrals = new CIIntegrals(*excited);

    integralsMBPT = dynamic_cast<CIIntegralsMBPT*>(integrals);
}

void Atom::SetMultipleIntegralsAndCore(unsigned int index)
{
    if(index >= multiple_ids.size())
    {   *errstream << "SetMultipleIntegralsAndCore: index " << index << " out of bounds" << std::endl;
        exit(1);
    }

    std::string original_identifier = identifier;

    core->ToggleOpenShellCore();
    identifier = multiple_ids[index];   // For updating core

    if(multiple_SMS)
    {   core->SetNuclearInverseMass(multiple_parameters[index]);
        *outstream << "\nNuclearInverseMass = " << multiple_parameters[index] << std::endl;
        integrals->IncludeValenceSMS(true);
    }
    else if(multiple_alpha)
    {   Constant::Alpha = alpha0 * sqrt(multiple_parameters[index]+1.);
        Constant::AlphaSquared = alpha0 * alpha0 * (multiple_parameters[index]+1.);
        *outstream << "\nx = " << multiple_parameters[index] << std::endl;
        integrals->IncludeValenceSMS(false);
    }
    else if(multiple_volume)
    {   core->SetVolumeShiftParameter(multiple_parameters[index]);
        *outstream << "\nVolumeShiftParameter = " << std::setprecision(3) << multiple_parameters[index] << std::endl;
        integrals->IncludeValenceSMS(false);
    }
    else if(multiple_radius)
    {   core->SetNuclearRadius(multiple_parameters[index]/Constant::AtomicToFermi);
        *outstream << "\nRadius = " << std::setprecision(3) << multiple_parameters[index] << std::endl;
        integrals->IncludeValenceSMS(false);
    }

    integrals->SetIdentifier(identifier);

    core->Update();
    excited->Update();
    if(excited_mbpt)
        excited_mbpt->Update();
    Read();
    //Write();
    core->ToggleClosedShellCore();

    integrals->Clear();
    integrals->Update();

    identifier = original_identifier;
}

void Atom::CalculateMultipleEnergies()
{
    SymmetryEigenstatesMap::iterator it = SymEigenstates.begin();

    while(it != SymEigenstates.end())
    {
        ConfigGenerator* conf_gen = GenerateConfigurations(it->first);

        for(unsigned int i = 0; i < multiple_ids.size(); i++)
        {
            Eigenstates* E = new Eigenstates(multiple_ids[i], conf_gen);
            SetMultipleIntegralsAndCore(i);

            if(!E->Read())
            {
                HamiltonianMatrix* H;

                #ifdef _MPI
                    H = new MPIHamiltonianMatrix(*integrals, conf_gen);
                #else
                    H = new HamiltonianMatrix(*integrals, conf_gen);
                #endif

                if(sigma3)
                    H->IncludeSigma3(sigma3);

                H->GenerateMatrix();
                //H->PollMatrix();

                #ifdef _SCALAPACK
                    H->WriteToFile("temp.matrix");
                    MPIHamiltonianMatrix* MpiH = dynamic_cast<MPIHamiltonianMatrix*>(H);
                    if(multiple_parameters[i] == 0.0)  // Get g-factors
                        MpiH->SolveScalapack("temp.matrix", -0.45, *E, true);
                    else
                        MpiH->SolveScalapack("temp.matrix", -0.45, *E);
                #else
                    if(multiple_parameters[i] == 0.0)  // Get g-factors
                        H->SolveMatrix(NumSolutions, *E, true);
                    else
                        H->SolveMatrix(NumSolutions, *E);
                #endif

                delete H;
            }
        }

        it++;
    }
}

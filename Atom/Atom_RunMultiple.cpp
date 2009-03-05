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

void Atom::RunMultiple(bool include_mbpt, bool closed_shell)
{
    generate_mbpt_integrals = false;

    DebugOptions.LogFirstBuild(false);
    DebugOptions.LogHFIterations(false);
    DebugOptions.OutputHFExcited(true);
    DebugOptions.HartreeEnergyUnits(true);
    DebugOptions.LogMBPT(false);

    //CreateCustomBasis(include_mbpt);
    //CreateRBasis(include_mbpt);
    CreateBSplineBasis(include_mbpt);
    //CreateHartreeFockBasis(include_mbpt);

    DebugOptions.OutputHFExcited(false);

    multiple_SMS = false;
    multiple_alpha = false;
    multiple_volume = false;
    multiple_radius = true;

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
        double parameter = -1.;

        for(i = 0; i < 5; i++)
        {
            std::stringstream ss;
            multiple_parameters.push_back(parameter);
            ss << (int)(parameter);
            multiple_ids.push_back(identifier + '_' + ss.str());

            parameter += .5;
        }

        double deltaR = 0.05;
        core->CalculateVolumeShiftPotential(deltaR/Constant::AtomicToFermi);
        *outstream << "\nThickness = " << std::setprecision(3) << core->GetNuclearThickness()*Constant::AtomicToFermi;
        *outstream << "\nRadius    = " << std::setprecision(3) << core->GetNuclearRadius()*Constant::AtomicToFermi;
        *outstream << "\nd(Radius) = " << std::setprecision(3) << deltaR << std::endl;
    }
    else if(multiple_radius)
    {   // Nuclear radius in fermi
        double parameter = 0.;
        
        multiple_parameters.push_back(6.80502);
        multiple_parameters.push_back(6.86598);
        multiple_parameters.push_back(6.9264);
        multiple_parameters.push_back(6.9863);
        multiple_parameters.push_back(7.04569);

        for(i = 0; i < 5; i++)
        {
            std::stringstream ss;
            //multiple_parameters.push_back(parameter);
            ss << (int)(parameter);
            multiple_ids.push_back(identifier + '_' + ss.str());

            parameter += 1.;
        }

        *outstream << "\nThickness = " << std::setprecision(3)
                   << core->GetNuclearThickness()*Constant::AtomicToFermi << std::endl;
    }

    original_id = identifier;

    if(closed_shell)
        CalculateMultipleClosedShell(include_mbpt);
    else
    {   // delta = E_CI - E_HF
        double delta[] = {-153032.343, -153171.873, -153311.989, -153452.693, -153593.985};
        for(unsigned int i = 0; i < 5; i++)
            mbpt_delta.push_back(delta[i]);

        // check_size_only should be off
        check_size_only = false;

        if(include_mbpt)
        {
            // Uncomment to include sigma3.
            //sigma3 = new Sigma3Calculator(lattice, core, excited);

            for(unsigned int i = 0; i < multiple_ids.size(); i++)
            {
                // Save orbitals
                SetMultipleIntegralsAndCore(i);

                if(generate_mbpt_integrals)
                {   // See RunOpen() for explanation of the order of these
                    CollateIntegralsMBPT(NumProcessors);
                    GenerateIntegralsMBPT(true, false, mbpt_delta[i]);
                    CollateIntegralsMBPT(NumProcessors);
                }
            }
        }

        GenerateMultipleIntegrals(include_mbpt);

        // Save integrals, if we haven't just generated mbpt integrals.
        if(!generate_mbpt_integrals)
        {
            for(unsigned int i = 0; i < multiple_ids.size(); i++)
            {
                SetMultipleIntegralsAndCore(i);

                if(ProcessorRank == 0)
                {   integrals->WriteOneElectronIntegrals(true);
                    integrals->WriteTwoElectronIntegrals(true);
                }
                #ifdef _MPI
                    // Wait for root node to finish writing
                    MPI::COMM_WORLD.Barrier();
                #endif
            }
        }

        ChooseSymmetries();
        CalculateMultipleEnergies();
    }
}

void Atom::GenerateMultipleIntegrals(bool MBPT_CI)
{
    if(integrals)
        delete integrals;

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

    core->ToggleOpenShellCore();
    identifier = multiple_ids[index];   // For updating core

    if(multiple_SMS)
    {   core->SetNuclearInverseMass(multiple_parameters[index]);
        *outstream << "\nNuclearInverseMass = " << multiple_parameters[index] << std::endl;
        if(integrals)
            integrals->IncludeValenceSMS(true);
    }
    else if(multiple_alpha)
    {   Constant::Alpha = alpha0 * sqrt(multiple_parameters[index]+1.);
        Constant::AlphaSquared = alpha0 * alpha0 * (multiple_parameters[index]+1.);
        *outstream << "\nx = " << multiple_parameters[index] << std::endl;
        if(integrals)
            integrals->IncludeValenceSMS(false);
    }
    else if(multiple_volume)
    {   core->SetVolumeShiftParameter(multiple_parameters[index]);
        *outstream << "\nVolumeShiftParameter = " << std::setprecision(3) << multiple_parameters[index] << std::endl;
        if(integrals)
            integrals->IncludeValenceSMS(false);
    }
    else if(multiple_radius)
    {   core->SetNuclearRadius(multiple_parameters[index]/Constant::AtomicToFermi);
        *outstream << "\nRadius = " << std::setprecision(3) << multiple_parameters[index] << std::endl;
        if(integrals)
            integrals->IncludeValenceSMS(false);
    }

    if(integrals)
        integrals->SetIdentifier(identifier);

    core->Update();
    excited->Update();
    if(excited_mbpt)
        excited_mbpt->Update();
    ReadOrWriteBasis();
    core->ToggleClosedShellCore();

    if(integrals)
    {   integrals->Clear();
        integrals->Update();
    }

    if(sigma3)
    {   sigma3->UpdateIntegrals(excited);
        sigma3->SetEnergyShift(mbpt_delta[index]/Constant::HartreeEnergy_cm);
    }
}

void Atom::CalculateMultipleEnergies()
{
    SymmetryEigenstatesMap::iterator it = SymEigenstates.begin();

    while(it != SymEigenstates.end())
    {
        identifier = original_id;
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
                    if(multiple_parameters[i])// == 0.0)  // Get g-factors
                        MpiH->SolveScalapack("temp.matrix", -1.43, *E, true);
                    else
                        MpiH->SolveScalapack("temp.matrix", -1.43, *E);
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

void Atom::CalculateMultipleClosedShell(bool include_mbpt)
{
    bool brueckner = true;

    // Report ordering of excited states
    ConstStateIterator it = excited->GetConstStateIterator();
    while(!it.AtEnd())
    {   *outstream << it.GetStateInfo().Name() << std::endl;
        it.Next();
    }

    if(include_mbpt)
         mbpt = new CoreMBPTCalculator(lattice, core, excited_mbpt);

    for(unsigned int i = 0; i < multiple_ids.size(); i++)
    {
        SetMultipleIntegralsAndCore(i);

        if(include_mbpt)
        {   if(brueckner)
            {   excited->SetIdentifier(multiple_ids[i]);
                excited->ClearSigmas();
            }
            else
                mbpt->UpdateIntegrals(excited);
        }

        // Calculate for all excited states
        it = excited->GetConstStateIterator();
        while(!it.AtEnd())
        {
            StateInfo si = it.GetStateInfo();
            if(include_mbpt && brueckner)
            {   excited->CreateSecondOrderSigma(si, *mbpt);
                DiscreteState ds = excited->GetStateWithSigma(si);
                *outstream << std::setprecision(15) << ds.Energy() * Constant::HartreeEnergy_cm << std::endl;
            }
            else
            {   const DiscreteState* ds = it.GetState();
                double dE = 0.;
                if(include_mbpt)
                    dE = mbpt->GetOneElectronDiagrams(si, si);
                *outstream << std::setprecision(15) << (ds->Energy() + dE)*Constant::HartreeEnergy_cm << std::endl;
            }

            it.Next();
        }
    }
}

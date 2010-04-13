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

void Atom::InitialiseRunIndex()
{
    // Set up any multiple run options that may be used
    multiple_length = 0;
    unsigned int multiple_options_used = 0;  // Number of multiple options being used.
    unsigned int length;
    unsigned int i;

    // Typical values for NuclearInverseMass: -0.002 -> 0.002
    length = userInput_.vector_variable_size("NuclearInverseMass");
    for(i = 0; i < length; i++)
    {   multiple_SMS.push_back(userInput_("NuclearInverseMass", 0.0, i));
    }
    multiple_length = mmax(multiple_length, length);
    if(length > 1)
    {   multiple_options_used++;
        multiple_parameters = &multiple_SMS;
    }

    // Typical values for AlphaSquaredVariation: -0.25 -> 0.25
    length = userInput_.vector_variable_size("AlphaSquaredVariation");
    for(i = 0; i < length; i++)
    {   multiple_alpha.push_back(userInput_("AlphaSquaredVariation", 0.0, i));
    }
    alpha0 = Constant::Alpha;
    multiple_length = mmax(multiple_length, length);
    if(length > 1)
    {   multiple_options_used++;
        multiple_parameters = &multiple_alpha;
    }

    // Typical values for NuclearVolumeVariation: -1. -> 1.
    // with R ~ 4 fm and dR = 0.05 fm
    length = userInput_.vector_variable_size("NuclearVolumeVariation");
    for(i = 0; i < length; i++)
    {   multiple_volume.push_back(userInput_("NuclearVolumeVariation", 0.0, i));
    }
    multiple_length = mmax(multiple_length, length);
    if(length > 1)
    {   multiple_options_used++;
        multiple_parameters = &multiple_volume;
    }

    // Typical values for NuclearRadius is a few fm.
    length = userInput_.vector_variable_size("NuclearRadius");
    for(i = 0; i < length; i++)
    {   multiple_radius.push_back(userInput_("NuclearRadius", 0.0, i));
    }
    multiple_length = mmax(multiple_length, length);
    if(length > 1)
    {   multiple_options_used++;
        multiple_parameters = &multiple_radius;
    }

    if(multiple_options_used > 1)
    {   *errstream << "Error: Used more than one multiple run option." << std::endl;
        exit(1);
    }

    original_id = identifier;

    // Choose which of the multiple run options are being used in the current calculation
    current_run_selection.clear();
    unsigned int total_run_selections = userInput_.vector_variable_size("-r");

    if((multiple_length <= 1) || userInput_.search("--check-sizes"))
    {   // Only do a single run if --check-sizes option is set (doesn't matter which one).
        current_run_selection.push_back(0);
    }
    else if(total_run_selections)
    {
        for(i = 0; i < total_run_selections; i++)
        {
            int selection = userInput_("-r", 0, i);

            if(selection == 0)
            {   // "-r=0" is shorthand for "just do the zero variation variety"
                if(total_run_selections > 1)
                {   *outstream << "USAGE: \"-r=0\" can only be used by itself." << std::endl;
                    exit(1);
                }

                for(unsigned int j = 0; j < multiple_length; j++)
                    if((*multiple_parameters)[j] == 0.0)
                    {   current_run_selection.push_back(j);
                        break;
                    }

                if(current_run_selection.size() == 0)
                {   *outstream << "ERROR: Multiple run option couldn't find zero-variation option." << std::endl;
                    exit(1);
                }
            }
            else if((selection < 0) || (selection > multiple_length))
            {   *outstream << "USAGE: Option \"-r\" included index outside of multiple run range." << std::endl;
                exit(1);
            }
            else
            {   // Subtract one from the user's selection to make it start at zero
                current_run_selection.push_back(selection-1);
            }
        }
    }
    else
    {   // If "-r" is not found, do all of the multiple run options
        for(i = 0; i < multiple_length; i++)
            current_run_selection.push_back(i);
    }
}

unsigned int Atom::NumberRunsSelected()
{
    return current_run_selection.size();
}

void Atom::RunIndexBegin(bool print)
{
    current_run_index = 0;
    if(multiple_length > 1)
        identifier = original_id + "_" + itoa(current_run_selection[current_run_index]);
    SetRunParameters(print);
}

void Atom::RunIndexNext(bool print)
{
    current_run_index++;
    if(!RunIndexAtEnd())
    {   identifier = original_id + "_" + itoa(current_run_selection[current_run_index]);
        SetRunParameters(print);
    }
}

bool Atom::RunIndexAtEnd()
{
    return (current_run_index >= NumberRunsSelected());
}

void Atom::InitialiseParameters()
{
    *outstream << "Physical Parameters:" << std::endl;

    // Nuclear Radius
    if(multiple_radius.size() == 1)
        core->SetNuclearRadius(multiple_radius[0]/Constant::AtomicToFermi);
    if(multiple_radius.size() <= 1)
        *outstream << "NuclearRadius = " << std::setprecision(3)
                   << core->GetNuclearRadius()*Constant::AtomicToFermi << std::endl;

    // Nuclear Thickness
    double nuclear_thickness = userInput_("NuclearThickness", -1.);
    if(nuclear_thickness < 0.0) // Not found
    {   if(core->GetNuclearRadius() < 1.5)
            nuclear_thickness = 0.0;
        else
            nuclear_thickness = 2.3;    // Standard value
    }
    core->SetNuclearThickness(nuclear_thickness/Constant::AtomicToFermi);
    *outstream << "NuclearThickness = " << std::setprecision(3) << nuclear_thickness << std::endl;

    if(multiple_volume.size())
    {   // Field shift radius
        double deltaR = userInput_("DeltaNuclearRadius", 0.05);
        core->CalculateVolumeShiftPotential(deltaR/Constant::AtomicToFermi);
        *outstream << "delta(NuclearRadius) = " << std::setprecision(3) << deltaR << std::endl;
    }

    if(multiple_SMS.size() == 1)
    {   core->SetNuclearInverseMass(multiple_SMS[0]);
        *outstream << "NuclearInverseMass = " << core->GetNuclearInverseMass() << std::endl;
    }

    if(multiple_alpha.size() == 1)
    {   alpha0 = Constant::Alpha;
        Constant::Alpha = alpha0 * sqrt(1.0 + multiple_alpha[0]);
        Constant::AlphaSquared = alpha0 * alpha0 * (1.0 + multiple_alpha[0]);
        *outstream << "AlphaSquaredVariation (x) = " << multiple_alpha[0] << std::endl;
   }
}

void Atom::SetRunParameters(bool print)
{
    unsigned int index = current_run_selection[current_run_index];

    if(multiple_SMS.size() > 1)
    {   core->SetNuclearInverseMass(multiple_SMS[index]);
        if(print)
            *outstream << "\nNuclearInverseMass = " << multiple_SMS[index] << std::endl;
        if(integrals)
        {   if(multiple_SMS[index])
                integrals->IncludeValenceSMS(true);
            else
                integrals->IncludeValenceSMS(false);
        }
    }
    if(multiple_alpha.size() > 1)
    {   Constant::Alpha = alpha0 * sqrt(multiple_alpha[index]+1.);
        Constant::AlphaSquared = alpha0 * alpha0 * (multiple_alpha[index]+1.);
        if(print)
            *outstream << "\nAlphaSquaredVariation (x) = " << multiple_alpha[index] << std::endl;
    }
    if(multiple_volume.size() > 1)
    {   core->SetVolumeShiftParameter(multiple_volume[index]);
        if(print)
            *outstream << "\nVolumeShiftParameter = " << std::setprecision(3) << multiple_volume[index] << std::endl;
    }
    if(multiple_radius.size() > 1)
    {   core->SetNuclearRadius(multiple_radius[index]/Constant::AtomicToFermi);
        if(print)
            *outstream << "\nNuclearRadius = " << std::setprecision(3) << multiple_radius[index] << std::endl;
    }
}

void Atom::SetRunCore(bool force)
{
    static int previous_index = -1;
    unsigned int index = current_run_selection[current_run_index];

    if(force || (previous_index != index))
    {
        core->ToggleOpenShellCore();

        // Try read, if it fails then update
        if(!useRead || !Read())
        {
            if(!useRead || !ReadCore())
                core->Update();
            if(excited)
                excited->Update();
            if(excited_mbpt)
                excited_mbpt->Update();
        }

        core->ToggleClosedShellCore();
        previous_index = index;
    }
}

void Atom::SetRunIntegrals(bool force)
{
    static int previous_index = -1;
    unsigned int index = current_run_selection[current_run_index];

    if(force || (previous_index != index))
    {
        unsigned int index = current_run_selection[current_run_index];

        double delta;
        if(mbpt || valence_mbpt || sigma3)
        {   if(mbpt_delta.size() == 0)
                delta = 0.0;
            else if(mbpt_delta.size() == 1)
                delta = mbpt_delta[0];
            else
                delta = mbpt_delta[index];
        }

        if(mbpt)
            mbpt->SetEnergyShift(delta);
        if(valence_mbpt)
            valence_mbpt->SetEnergyShift(delta);

        if(integrals)
        {   integrals->SetIdentifier(identifier);
            integrals->Clear();
            integrals->Update();
        }

        if(sigma3)
        {   sigma3->UpdateIntegrals(excited);
            sigma3->SetEnergyShift(delta);
        }
        previous_index = index;
    }
}
/*
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
    SymmetryEigenstatesMap::iterator it = symEigenstates.begin();

    while(it != symEigenstates.end())
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
*/
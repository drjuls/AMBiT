#include "Include.h"
#include "CustomBasis.h"
#include "BasisGenerator.h"
#include "HartreeFock/HartreeFocker.h"
#include "HartreeFock/ConfigurationParser.h"
#include "Universal/MathConstant.h"
#include "boost/algorithm/string/split.hpp"
#include "boost/algorithm/string/classification.hpp"
#include "boost/algorithm/string/trim.hpp"

namespace Ambit
{
// This file contains the basis creation function from BasisGenerator as well as CustomBasis
pOrbitalMap BasisGenerator::GenerateXRExcited(const std::vector<int>& max_pqn)
{
    pOrbitalMap excited(new OrbitalMap(lattice));

    if(!max_pqn.size())
        return excited;

    bool debug = DebugOptions.OutputHFExcited();

    // Create Hartree-Fock solver; define integrators.
    pIntegrator integrator(new SimpsonsIntegrator(lattice));
    pODESolver ode_solver(new AdamsSolver(integrator));
    HartreeFocker HF_Solver(ode_solver);

    // Do HF orbitals first. This ensures the lattice is large enough and that the starting levels exist.
    unsigned int numorbitals = user_input.vector_variable_size("Basis/CustomOrbitals");
    for(unsigned int i = 0; i < numorbitals; i++)
    {
        // Split each instruction on space
        std::vector<std::string> instruction;
        boost::split(instruction, boost::algorithm::trim_copy(user_input("Basis/CustomOrbitals", "", i)),
                     boost::is_any_of(" -"), boost::token_compress_on);

        // Check orbital is required for current max_pqn
        NonRelInfo nrcurrent = ConfigurationParser::ParseOrbital(instruction[0]);
        if((max_pqn.size() <= nrcurrent.L()) || (max_pqn[nrcurrent.L()] < nrcurrent.PQN()))
            continue;

        if(instruction.size() == 1)
        {
            // Create orbitals using Hartree Fock
            for(OrbitalInfo current: nrcurrent.GetRelativisticInfos())
            {
                pOrbitalConst s = open_core->GetState(current);
                if(s)
                {   // Already exists in core
                    if(!closed_core->GetOccupancy(current))
                    {
                        pOrbital s_copy = std::make_shared<Orbital>(s);
                        double nu = s->Nu();
                        *s_copy = *s;
                        excited->AddState(s_copy);

                        if(debug)
                            *logstream << "  " << s_copy->Name() << " en:   " << std::setprecision(12) << s_copy->Energy() << "  size:  " << s_copy->size() << std::endl;
                    }
                }
                else
                {
                    pOrbital ds = std::make_shared<Orbital>(current);

                    // Look for orbital with (pqn-1, kappa)
                    pOrbitalConst onelower = open_core->GetState(OrbitalInfo(current.PQN()-1, current.Kappa()));
                    if(!onelower)
                        onelower = excited->GetState(OrbitalInfo(current.PQN()-1, current.Kappa()));

                    if(onelower)
                        ds->SetNu(onelower->Nu() + 1./hf->GetCharge());

                    HF_Solver.CalculateExcitedState(ds, hf);
                    double nu = ds->Nu();
                    excited->AddState(ds);

                    if(debug)
                        *logstream << "  " << ds->Name() << " en:   " << std::setprecision(12) << ds->Energy() << "  size:  " << ds->size() << std::endl;
                }
            }
        }
    }

    // Do times R orbitals
    CustomBasis timesr(hf);

    for(unsigned int i = 0; i < numorbitals; i++)
    {
        // Split each instruction on space
        std::vector<std::string> instruction;
        boost::split(instruction, boost::algorithm::trim_copy(user_input("Basis/CustomOrbitals", "", i)),
                     boost::is_any_of(" -"), boost::token_compress_on);

        // Check orbital is required for current max_pqn
        NonRelInfo nrcurrent = ConfigurationParser::ParseOrbital(instruction[0]);
        if((max_pqn.size() <= nrcurrent.L()) || (max_pqn[nrcurrent.L()] < nrcurrent.PQN()))
            continue;

        if (instruction.size() == 3)
        {
            NonRelInfo nrprev = ConfigurationParser::ParseOrbital(instruction[2]);

            OrbitalInfo current = nrcurrent.GetFirstRelativisticInfo();
            OrbitalInfo prev = nrprev.GetFirstRelativisticInfo();

            int pass = 0;
            while(pass < 2)
            {
                pass++;
                pOrbital ds = open_core->GetState(current);
                if(ds == NULL)
                {   ds = std::make_shared<Orbital>(current);
                    excited->AddState(ds);
                }

                pOrbitalConst previous = excited->GetState(prev);
                if(previous == NULL)
                    previous = open_core->GetState(prev);
                if(previous == NULL)
                {   *errstream << "Input error: " << prev.Name() << " undefined." << std::endl;
                    exit(1);
                }

                if(instruction[1] == std::string("R"))
                    timesr.MultiplyByR(previous, ds);
                else if(instruction[1] == std::string("S"))
                    timesr.MultiplyBySinR(previous, ds);
                else if(instruction[1] == std::string("RS"))
                    timesr.MultiplyByRSinR(previous, ds);
                else
                {   *errstream << "Input error: " << user_input("Basis/CustomOrbitals", "", i) << " unknown." << std::endl;
                    exit(1);
                }

                // Orthogonalise and set energy of new state
                Orthogonalise(ds);
                Orthogonalise(ds, excited);

                if(DebugOptions.OutputHFExcited())
                {   unsigned int count, r;
                    double fmax = 0.;
                    for(count = 0; count < ds->size(); count++)
                        if(fabs(ds->f[count]) > fmax)
                        {   fmax = fabs(ds->f[count]);
                            r = count;
                        }

                    *outstream << "  " << ds->Name() << "  en: " << std::setprecision(8) << ds->Energy()
                               << "  size: " << ds->size() << "  Rmax: " << r << "  NZ: "
                               << int(ds->PQN()) - int(ds->NumNodes()) - int(ds->L()) - 1 << std::endl;
                }

                if(current.L() != 0)
                {
                    current = nrcurrent.GetSecondRelativisticInfo();
                    prev = nrprev.GetSecondRelativisticInfo();
                }
                else
                    pass++;
            }
        }
    }

    return excited;
}

CustomBasis::CustomBasis(pHFOperator hf_core)
{
    pIntegrator integrator = std::make_shared<SimpsonsIntegrator>(hf_core->GetLattice());
    pODESolver ode_solver = std::make_shared<AdamsSolver>(integrator);
    pCoulombOperator coulomb = std::make_shared<CoulombOperator>(hf_core->GetLattice(), ode_solver);

    if(hf_core->IncludeExchange())
    {
        // Replace exchange with local approximation to get effective direct potential
        pHFOperator localexch = std::make_shared<LocalExchangeApproximation>(hf_core, coulomb, 1.0);
        localexch->SetCore(hf_core->GetCore());
        potential = localexch->GetDirectPotential();
    }
    else
        potential = hf_core->GetDirectPotential();

    hf = hf_core;
}

void CustomBasis::MultiplyByR(pOrbitalConst previous, pOrbital current) const
{
    current->resize(previous->size());

    const double* R = hf->GetLattice()->R();
    double alpha = hf->GetPhysicalConstant()->GetAlpha();
    double AlphaSquared = hf->GetPhysicalConstant()->GetAlphaSquared(); // Set to zero to remove effect of core potential

    double kappa_c = current->Kappa();
    double kappa_p = previous->Kappa();

    for(unsigned int i=0; i<previous->size(); i++)
    {
        current->f[i] = previous->f[i] * R[i];
        current->g[i] = alpha/(2. + AlphaSquared*potential.f[i])
                * (R[i]*previous->dfdr[i] + (1. + kappa_c)*previous->f[i]);
        current->dfdr[i] = previous->dfdr[i] * R[i] + previous->f[i];
        current->dgdr[i] = previous->dgdr[i] * R[i]
                + alpha/(2. + AlphaSquared*potential.f[i])
                  * (kappa_p/R[i] * previous->f[i] + (2. + (kappa_c - kappa_p)/R[i])*previous->dfdr[i]
                     + alpha * potential.dfdr[i] * (previous->g[i]*R[i] - current->g[i]));
    }

    current->ReNormalise(hf->GetIntegrator());
}

void CustomBasis::MultiplyBySinR(pOrbitalConst previous, pOrbital current) const
{
    current->resize(previous->size());

    const double* R = hf->GetLattice()->R();
    double alpha = hf->GetPhysicalConstant()->GetAlpha();
    double AlphaSquared = hf->GetPhysicalConstant()->GetAlphaSquared(); // Set to zero to remove effect of core potential

    double kappa_c = current->Kappa();
    double kappa_p = previous->Kappa();
    double k = MathConstant::Instance()->Pi()/hf->GetLattice()->R(previous->size());

    for(unsigned int i=0; i<previous->size(); i++)
    {
        double sinp = sin(k*R[i]);
        double kcosp = k * cos(k*R[i]);

        current->f[i] = previous->f[i] * sinp;
        current->g[i] = alpha/(2. + AlphaSquared*potential.f[i])
                * (sinp*previous->dfdr[i] + (kcosp + kappa_c/R[i]*sinp)*previous->f[i]);
        current->dfdr[i] = previous->dfdr[i] * sinp + previous->f[i] * kcosp;
        current->dgdr[i] = previous->dgdr[i] * sinp
                + alpha/(2. + AlphaSquared*potential.f[i])
                  * ((kappa_c/R[i]*kcosp - (k*k + (kappa_c-kappa_p)/(R[i]*R[i])) * sinp) * previous->f[i]
                     + (2.* kcosp + (kappa_c-kappa_p)/R[i] * sinp) * previous->dfdr[i]
                     + alpha * potential.dfdr[i] * (previous->g[i]*sinp - current->g[i]));
    }

    current->ReNormalise(hf->GetIntegrator());
}

void CustomBasis::MultiplyByRSinR(pOrbitalConst previous, pOrbital current) const
{
    current->resize(previous->size());

    const double* R = hf->GetLattice()->R();
    double alpha = hf->GetPhysicalConstant()->GetAlpha();
    double AlphaSquared = hf->GetPhysicalConstant()->GetAlphaSquared(); // Set to zero to remove effect of core potential

    double kappa_c = current->Kappa();
    double kappa_p = previous->Kappa();
    double k = MathConstant::Instance()->Pi()/hf->GetLattice()->R(previous->size());

    for(unsigned int i=0; i<previous->size(); i++)
    {
        double sinp = sin(k*R[i]);
        double kcosp = k * cos(k*R[i]);

        current->f[i] = previous->f[i] * R[i] * sinp;
        current->g[i] = alpha/(2. + AlphaSquared*potential.f[i])
                * (R[i] * sinp * previous->dfdr[i]
                        + (R[i] * kcosp + (1.+ kappa_c) * sinp)*previous->f[i]);
        current->dfdr[i] = previous->dfdr[i] * R[i] * sinp
                + previous->f[i] * (R[i] * kcosp + sinp);
        current->dgdr[i] = previous->dgdr[i] * R[i] * sinp
                + alpha/(2. + AlphaSquared*potential.f[i])
                  * (((kappa_p/R[i] - R[i]*k*k) * sinp + (2.+ kappa_c) * kcosp) * previous->f[i]
                     + (2.*R[i]*kcosp + (2. + kappa_c - kappa_p) * sinp) * previous->dfdr[i]
                     + alpha * potential.dfdr[i] * (previous->g[i]*R[i]*sinp - current->g[i]));
    }

    current->ReNormalise(hf->GetIntegrator());
}
}

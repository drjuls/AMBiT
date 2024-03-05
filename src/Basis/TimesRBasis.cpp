#include "Include.h"
#include "TimesRBasis.h"
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

    // Create list of custom orbitals
    typedef std::tuple<NonRelInfo, std::string, NonRelInfo> XRInstruction;
    auto XRInstructionComparator = [](const XRInstruction& a, const XRInstruction& b)
    {   return std::get<0>(a) < std::get<0>(b);
    };
    auto XRInstructionFromNonRelInfo = [](const NonRelInfo& nrinfo)
    {   return std::make_tuple(nrinfo, "", nrinfo);
    };

    std::vector<XRInstruction> xRinstructions;
    std::vector<NonRelInfo> HFinstructions;

    unsigned int numcustomorbitals = user_input.vector_variable_size("Basis/CustomOrbitals");
    for(unsigned int i = 0; i < numcustomorbitals; i++)
    {
        std::vector<std::string> line;
        std::string trimmed = boost::algorithm::trim_copy(user_input("Basis/CustomOrbitals", "", i));
        boost::algorithm::split(line, trimmed, boost::is_any_of(" -"), boost::token_compress_on);

        NonRelInfo nrcurrent = ConfigurationParser::ParseOrbital(line[0]);
        if(line.size() == 1)
            HFinstructions.push_back(nrcurrent);
        else if (line.size() == 3)
        {   NonRelInfo nrprev = ConfigurationParser::ParseOrbital(line[2]);
            xRinstructions.push_back(std::make_tuple(nrcurrent, line[1], nrprev));
        }
        else
        {   *errstream << "Input error: CustomOrbitals: " << user_input("Basis/CustomOrbitals", "", i) << '\n';
        }
    }

    // Sort instruction lists for searching on target NonRelInfo
    std::sort(HFinstructions.begin(), HFinstructions.end());
    std::sort(xRinstructions.begin(), xRinstructions.end(), XRInstructionComparator);
    auto xrend = std::unique(xRinstructions.begin(), xRinstructions.end());
    xRinstructions.resize(std::distance(xRinstructions.begin(), xrend));

    // Default first state in each kappa is HF, so add to HF list if it is not in CustomOrbitals
    for(int l = 0; l < max_pqn.size(); l++)
    {
        if(!max_pqn[l])
            continue;

        for(int kappa = - (l+1); kappa <= l; kappa += 2*l + 1)
        {
            if(kappa == 0)
                break;

            unsigned int pqn = l + 1;
            while(pqn <= max_pqn[l] && closed_core->GetState(OrbitalInfo(pqn, kappa)))
                pqn++;

            if(pqn <= max_pqn[l])
            {
                if(!std::binary_search(xRinstructions.begin(), xRinstructions.end(),
                                       XRInstructionFromNonRelInfo(NonRelInfo(pqn, kappa)), XRInstructionComparator))
                    HFinstructions.push_back(OrbitalInfo(pqn, kappa));
            }
        }
    }

    std::sort(HFinstructions.begin(), HFinstructions.end());
    auto hfend = std::unique(HFinstructions.begin(), HFinstructions.end());
    HFinstructions.resize(std::distance(HFinstructions.begin(), hfend));

    // Create Hartree-Fock solver; define integrators.
    pIntegrator integrator(new SimpsonsIntegrator(lattice));
    pODESolver ode_solver(new AdamsSolver(integrator));
    HartreeFocker HF_Solver(ode_solver);

    // Do HF orbitals first. This ensures the lattice is large enough and that the starting levels exist.
    std::string hf_valence_states = user_input("Basis/HFOrbitals", "");
    if(!hf_valence_states.empty())
    {
        pOrbitalMap hf_valence = GenerateHFExcited(ConfigurationParser::ParseBasisSize(hf_valence_states));
        for(auto porb: *hf_valence)
            excited->AddState(porb.second);
    }

    for(auto nrcurrent: HFinstructions)
    {
        // Create orbitals using Hartree Fock
        for(OrbitalInfo current: nrcurrent.GetRelativisticInfos())
        {
            if(excited->GetState(current))
                continue;

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

    // Do times R orbitals
    TimesRBasis timesr(hf);

    // Default first state in each kappa is HF, so add to HF list if it is not in CustomOrbitals
    for(int l = 0; l < max_pqn.size(); l++)
    {
        if(!max_pqn[l])
            continue;

        for(int kappa = - (l+1); kappa <= l; kappa += 2*l + 1)
        {
            if(kappa == 0)
                break;

            unsigned int pqn = l + 1;
            bool use_sinkr = false;
            while(pqn <= max_pqn[l])
            {
                OrbitalInfo current = OrbitalInfo(pqn, kappa);
                OrbitalInfo prev;
                std::string method;
                pqn++;

                if(closed_core->GetState(current) || excited->GetState(current))
                {   use_sinkr = false;
                    continue;
                }
                // Get method and previous orbitalinfo
                else if(std::binary_search(xRinstructions.begin(), xRinstructions.end(),
                                           XRInstructionFromNonRelInfo(current), XRInstructionComparator))
                {
                    auto it = std::lower_bound(xRinstructions.begin(), xRinstructions.end(),
                                               XRInstructionFromNonRelInfo(current), XRInstructionComparator);

                    // Use custom instruction
                    method = std::get<1>(*it);
                    NonRelInfo nrprev = std::get<2>(*it);
                    prev = (kappa < 0? nrprev.GetFirstRelativisticInfo(): nrprev.GetSecondRelativisticInfo());
                    use_sinkr = false;
                }
                else
                {   prev = OrbitalInfo(current.PQN()-1, kappa);
                    method = (use_sinkr? "S": "R");
                    use_sinkr = !use_sinkr;
                }

                // Make new orbital
                pOrbital ds = std::make_shared<Orbital>(current);
                excited->AddState(ds);

                // Find previous orbital
                pOrbitalConst previous = excited->GetState(prev);
                if(previous == NULL)
                    previous = open_core->GetState(prev);
                if(previous == NULL)
                {   *errstream << "Input error: " << prev.Name() << " not found." << std::endl;
                    exit(1);
                }

                if(method == std::string("R"))
                    timesr.MultiplyByR(previous, ds);
                else if(method == std::string("S"))
                    timesr.MultiplyBySinR(previous, ds);
                else if(method == std::string("RS"))
                    timesr.MultiplyByRSinR(previous, ds);
                else
                {   *errstream << "Input error: method " << method << " unknown." << std::endl;
                    exit(1);
                }

                // Orthogonalise and set energy of new state
                Orthogonalise(ds, open_core);
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
            }
        }
    }

    return excited;
}

TimesRBasis::TimesRBasis(pHFOperator hf_core)
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

void TimesRBasis::MultiplyByR(pOrbitalConst previous, pOrbital current) const
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
                  * (kappa_p/R[i] * previous->f[i] + (2. + kappa_c - kappa_p)*previous->dfdr[i]
                     + alpha * potential.dfdr[i] * (previous->g[i]*R[i] - current->g[i]));
    }

    current->ReNormalise(hf->GetIntegrator());
}

void TimesRBasis::MultiplyBySinR(pOrbitalConst previous, pOrbital current) const
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

void TimesRBasis::MultiplyByRSinR(pOrbitalConst previous, pOrbital current) const
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

#include "LorentzInvarianceT2.h"

namespace Ambit
{
SpinorFunction LorentzInvarianceT2Operator::ReducedApplyTo(const SpinorFunction& a, int kappa_b) const
{
    SpinorFunction ret(kappa_b);
    MathConstant* math = MathConstant::Instance();

    // prefactor = A(j', j) or B(j', j) or C(j', j)
    // from Hohensee et al. PRL 111 050401 (2013), Supplemental Material
    double prefactor = 0.;
    int two_jb = ret.TwoJ();
    int two_jmin = mmin(two_jb, a.TwoJ());

    if(abs(kappa_b + a.Kappa()) == 1)
    {   // A(j', j)
        prefactor = 6. * (two_jmin + 3) * (two_jmin + 1)/
                    double((two_jmin + 2) * two_jmin * (two_jmin + 4));
        prefactor = std::sqrt(prefactor) * math->minus_one_to_the_power((two_jmin - two_jb)/2 + 1);
    }
    else if(abs(kappa_b - a.Kappa()) == 2)
    {   // B(j', j)
        prefactor = 3. * (two_jmin + 5) * (two_jmin + 3) * (two_jmin + 1)/
                    double(2 * (two_jmin + 4) * (two_jmin + 2));
        prefactor = std::sqrt(prefactor) * math->minus_one_to_the_power((two_jmin - two_jb)/2 + 1);
    }
    else if(kappa_b == a.Kappa())
    {   // C(j', j)
        prefactor = (two_jmin + 3) * (two_jmin + 1) * (two_jmin - 1)/
                    double(two_jmin * (two_jmin + 2));
        prefactor = std::sqrt(prefactor);
    }

    if(prefactor)
    {
        ret.resize(a.size());

        // All T2 operators have the form
        //     ret.f = S a.dgdr + T a.g/r
        //     ret.g = U a.dfdr + V a.f/r
        // where S, T, U, V are constants depending on kappas.
        int S = 0, T = 0, U = 0, V = 0;
        const double kap = a.Kappa();

        if(kappa_b == -a.Kappa() - 1)
        {   // I1
            S = 2 * kap + 3;
            T = - kap * (2 * kap + 3);
            U = 2 * kap - 1;
            V = - (2 * kap - 1) * (kap + 1);
            prefactor *= -0.5;
        }
        else if (kappa_b == -a.Kappa() + 1)
        {   // I2
            S = 2 * kap + 1;
            T = (2 * kap + 1) * (kap - 1);
            U = 2 * kap - 3;
            V = kap * (2 * kap - 3);
            prefactor *= -0.5;
        }
        else if (kappa_b == a.Kappa() - 2)
        {   // I3
            S = 1;
            T = kap - 1;
            prefactor *= -2.;
        }
        else if (kappa_b == a.Kappa() + 2)
        {   // I4
            U = 1;
            V = -(kap + 1);
            prefactor *= 2.;
        }
        else if (kappa_b == a.Kappa())
        {   // I5
            S = -1;
            T = kap;
            U = 1;
            V = kap;
        }

        const double* R = lattice->R();
        const double* R2 = lattice->Rpower(2);

        if(S)
        {   std::vector<double> second_derivative_g(a.size());
            differentiator->GetDerivative(a.dgdr, second_derivative_g);

            for(unsigned int i = 0; i < a.size(); i++)
            {
                ret.f[i] = S * a.dgdr[i] + T * a.g[i]/R[i];
                ret.dfdr[i] = S * second_derivative_g[i] + T * (a.dgdr[i]/R[i] - a.g[i]/R2[i]);
            }
        }

        if(U)
        {   std::vector<double> second_derivative_f(a.size());
            differentiator->GetDerivative(a.dfdr, second_derivative_f);

            for(unsigned int i = 0; i < a.size(); i++)
            {
                ret.g[i] = U * a.dfdr[i] + V * a.f[i]/R[i];
                ret.dgdr[i] = U * second_derivative_f[i] + V * (a.dfdr[i]/R[i] - a.f[i]/R2[i]);
            }
        }

        ret *= prefactor * math->SpeedOfLightAU();
    }

    return ret;
}

LorentzInvarianceT2Calculator::LorentzInvarianceT2Calculator(MultirunOptions& user_input, Atom& atom):
    TransitionCalculator(user_input, atom.GetBasis(), atom.GetLevels())
{
    auto hf = atom.GetHFOperator();
    op = std::make_shared<LorentzInvarianceT2Operator>(hf->GetIntegrator());

    if(user_input.search("--rpa"))
    {   op = MakeRPA(std::static_pointer_cast<LorentzInvarianceT2Operator>(op), hf, atom.GetHartreeY());
        double scale = user_input("Scale", 1.);
        std::static_pointer_cast<RPAOperator>(op)->SetScale(scale);
    }
}

void LorentzInvarianceT2Calculator::PrintHeader() const
{
    if(user_input.search("--reduced-elements"))
        *outstream << "Lorentz reduced matrix elements (a.u.): " << std::endl;
    else
        *outstream << "Lorentz matrix elements (stretched states) in a.u.: " << std::endl;
}

void LorentzInvarianceT2Calculator::PrintTransition(const LevelID& left, const LevelID& right, double matrix_element) const
{
    MathConstant* math = MathConstant::Instance();
    double value = matrix_element;
    if(user_input.search("--reduced-elements"))
    {   int twoj1 = left.first->GetTwoJ();
        int twoj2 = right.first->GetTwoJ();
        value = value/math->Electron3j(twoj2, twoj1, op->GetK(), twoj2, -twoj1);
    }

    *outstream << "  " << Name(left) << " -> " << Name(right)
               << " = " << std::setprecision(6) << value << std::endl;
}

}

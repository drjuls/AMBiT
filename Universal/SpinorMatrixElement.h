#ifndef SPINOR_MATRIX_ELEMENT_H
#define SPINOR_MATRIX_ELEMENT_H

#include "Universal/SpinorFunction.h"
#include "HartreeFock/Orbital.h"
#include "HartreeFock/OrbitalInfo.h"
#include "Universal/MathConstant.h"
#include "Include.h"

/** SpinorMatrixElement is a base class for calculating radial matrix elements with an orbital basis. */
class SpinorMatrixElement : public std::enable_shared_from_this<SpinorMatrixElement>
{
public:
    SpinorMatrixElement(int K, Parity P, pIntegrator integration_strategy = nullptr):
        K(K), P(P), integrator(integration_strategy)
    {}
    /** Initialise with P = (-1)^K */
    SpinorMatrixElement(int K, pIntegrator integration_strategy = nullptr):
        integrator(integration_strategy), K(K)
    {   P = (K%2? Parity::odd: Parity::even);
    }

    /** Matrix element < b | t(K) | a > for our operator t(K) assuming stretched states
        for a and b (i.e. m_a = j_a and m_b = j_b) and q set to facilitate this.
        Default behaviour: take GetReducedMatrixElement(b, a) and scale.
     */
    virtual double GetMatrixElement(const Orbital& b, const Orbital& a) const
    {
        double reduced_matrix_element = GetReducedMatrixElement(b, a);
        return reduced_matrix_element * MathConstant::Instance()->Electron3j(a.TwoJ(), b.TwoJ(), K, a.TwoJ(), -b.TwoJ());
    }

    /** Reduced matrix element < b || t(K) || a > for our operator t(K).
        Usually the operation t|a> makes sense for any Orbital, but in order to have a
        matrix element one of the wavefunctions must be bounded, and hence an "Orbital".
     */
    virtual double GetReducedMatrixElement(const Orbital& b, const Orbital& a) const
    {   return 0.;
    }

    /** Get multipolarity K for this operator. */
    virtual int GetK() const
    {   return K;
    }

    /** Get parity of this operator. */
    virtual Parity GetParity() const
    {   return P;
    }

    virtual pIntegrator GetIntegrator() const
    {   return integrator;
    }

    /** Return false if matrix element is zero by symmetry (angular) properties only. */
    bool IsNonZero(const Orbital& b, const Orbital& a) const
    {
        if(K == 0)
        {   if(P == Parity::even)
                return a.Kappa() == b.Kappa();
            else
                return a.Kappa() == -b.Kappa();
        }
        else
            return ((abs(a.TwoJ() - b.TwoJ()) <= 2 * K) && (a.GetParity() * b.GetParity() == P));
    }

    /** Return false if matrix element is zero by symmetry (angular) properties only. */
    bool IsNonZero(const OrbitalInfo& b, const OrbitalInfo& a) const
    {
        if(K == 0)
        {   if(P == Parity::even)
                return a.Kappa() == b.Kappa();
            else
                return a.Kappa() == -b.Kappa();
        }
        else
            return ((abs(a.TwoJ() - b.TwoJ()) <= 2 * K) && (a.GetParity() * b.GetParity() == P));
    }

protected:
    pIntegrator integrator;
    int K;
    Parity P;
};

typedef std::shared_ptr<SpinorMatrixElement> pSpinorMatrixElement;
typedef std::shared_ptr<const SpinorMatrixElement> pSpinorMatrixElementConst;

#endif

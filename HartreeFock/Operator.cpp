#include "Operator.h"

double OneBodyOperator::GetMatrixElement(const SpinorFunction& b, const SpinorFunction& a) const
{
    SpinorFunction ta = ApplyTo(a, b.Kappa());
    if(!integrator)
        throw "OneBodyOperator::GetMatrixElement(): no integrator found.";
    return integrator->GetInnerProduct(ta, b);
}

SpinorFunction OneBodyOperator::ApplyTo(const SpinorFunction& a, int kappa_b) const
{
    if(kappa_b == a.Kappa())
        return ApplyTo(a);
    else
        return SpinorFunction(kappa_b);
}


double ZeroOperator::GetMatrixElement(const SpinorFunction& a, const SpinorFunction& b) const
{
    return 0.0;
}

SpinorFunction ZeroOperator::ApplyTo(const SpinorFunction& a) const
{
    return SpinorFunction(a.Kappa());
}

SpinorFunction ZeroOperator::ApplyTo(const SpinorFunction& a, int kappa_b) const
{
    return SpinorFunction(kappa_b);
}

#include "Operator.h"

double OneBodyOperator::GetMatrixElement(const SpinorFunction& b, const SpinorFunction& a) const
{
    SpinorFunction ta = ApplyTo(a);
    if(!integrator)
        throw "OneBodyOperator::GetMatrixElement(): no integrator found.";
    return integrator->GetInnerProduct(ta, b);
}

double ZeroOperator::GetMatrixElement(const SpinorFunction& a, const SpinorFunction& b) const
{
    return 0.0;
}

SpinorFunction ZeroOperator::ApplyTo(const SpinorFunction& a) const
{
    return SpinorFunction(a.Kappa(), a.Size());
}
